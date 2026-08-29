"""
gunw_reader.py
--------------
Read NISAR L2 GUNW products into line-of-sight displacement.

    unwrapped phase  ->  LOS displacement (mm)
    coherence        ->  quality gate
    connectedComponents -> unwrapping-reliability gate
    AOI polygon      ->  spatial clip
    reference region ->  removes the arbitrary phase offset

Why paths are discovered, not hardcoded
=======================================
NISAR L2 went public in 2026 and the group layout has moved between product
versions. Rather than hardcode /science/LSAR/GUNW/grids/... and fail on a
version bump, this walks the file and finds datasets by name. Run --inspect
on your first file to see exactly what was found and confirm it is sane.

Sign convention
===============
    d_los = -(lambda / 4*pi) * unwrapped_phase

Positive = motion AWAY from the satellite (range increase). For a descending
pass on a west-facing slope that usually reads as downslope motion, but ALWAYS
confirm against a known signal before interpreting. Use --flip-sign if your
convention is the opposite.

Unwrapped phase is RELATIVE. It carries an arbitrary constant per connected
component, so an absolute displacement number is meaningless until you
reference it to something you believe is stable. --ref-lat/--ref-lon does that.
Without a reference the script warns and reports relative values only.

Usage
-----
    python gunw_reader.py --inspect FILE.h5
    python gunw_reader.py --read FILE.h5 --ref-lat 28.21 --ref-lon 28.55
    python gunw_reader.py --read FILE.h5 --geotiff out.tif --quicklook out.png
    python gunw_reader.py --batch ./nisar_l2 --csv timeseries.csv

Requires: h5py, numpy.  Optional: rasterio + pyproj (GeoTIFF), matplotlib (PNG).
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

try:
    import h5py
except ImportError:
    print("h5py is required:  pip install h5py", file=sys.stderr)
    raise SystemExit(1)

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("gunw")

SPEED_OF_LIGHT = 299_792_458.0
NISAR_L_BAND_HZ = 1_257_500_000.0          # fallback if not in the file
DEFAULT_COHERENCE = 0.3

# Langtang AOI. Replace with LHENDE_RING for the source-zone box.
AOI_RING = [
    (85.46683434336315, 28.324709534140283),
    (85.45910958140026, 28.277704299412660),
    (85.47267083017955, 28.253664665263376),
    (85.50734379401987, 28.245345485689977),
    (85.53583958259410, 28.244740597906280),
    (85.55884220710581, 28.249882034713380),
    (85.56485035529917, 28.272864249930198),
    (85.56193211189097, 28.289493024692200),
    (85.55918552985972, 28.307177143980763),
    (85.54854252448862, 28.316849264798610),
    (85.52468159309214, 28.328333765691790),
    (85.50614216438120, 28.329693690212682),
]


# ---------------------------------------------------------------------------
# HDF5 discovery
# ---------------------------------------------------------------------------
def walk(h5file) -> dict[str, tuple]:
    """Map every dataset path -> (shape, dtype)."""
    found: dict[str, tuple] = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            found[name] = (obj.shape, obj.dtype)

    h5file.visititems(visit)
    return found


def find(datasets: dict, *patterns: str, ndim: int | None = None) -> str | None:
    """
    First dataset whose path matches all patterns (case-insensitive substrings).
    Prefers frequencyA and HH/VV polarisations when several match.
    """
    hits = []
    for path, (shape, _dtype) in datasets.items():
        low = path.lower()
        if all(p.lower() in low for p in patterns):
            if ndim is not None and len(shape) != ndim:
                continue
            hits.append(path)
    if not hits:
        return None
    hits.sort(key=lambda p: (
        0 if "frequencya" in p.lower() else 1,
        0 if re.search(r"/(hh|vv)/", p.lower()) else 1,
        len(p),
    ))
    return hits[0]


def read_wavelength(h5file, datasets: dict) -> float:
    """Radar wavelength in metres, from the file if possible."""
    for key in ("centerfrequency", "center_frequency", "radarcenterfrequency"):
        path = find(datasets, key)
        if path:
            try:
                freq = float(np.ravel(h5file[path][()])[0])
                if 1e8 < freq < 1e11:
                    lam = SPEED_OF_LIGHT / freq
                    logger.info("Wavelength from file: %.4f m (%.1f MHz)", lam, freq / 1e6)
                    return lam
            except Exception:
                pass
    lam = SPEED_OF_LIGHT / NISAR_L_BAND_HZ
    logger.warning("Center frequency not found; assuming NISAR L-band %.4f m", lam)
    return lam


def read_pair_dates(path: Path) -> tuple[str, str]:
    stamps = re.findall(r"_(\d{8})T\d{6}", path.name)
    if len(stamps) >= 4:
        return stamps[0], stamps[2]
    return "unknown", "unknown"


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
def point_in_ring(x: float, y: float, ring: list[tuple[float, float]]) -> bool:
    """Ray-casting point-in-polygon. No shapely dependency."""
    inside = False
    n = len(ring)
    for i in range(n):
        x1, y1 = ring[i]
        x2, y2 = ring[(i + 1) % n]
        if (y1 > y) != (y2 > y):
            xin = (x2 - x1) * (y - y1) / (y2 - y1) + x1
            if x < xin:
                inside = not inside
    return inside


def build_aoi_mask(xs: np.ndarray, ys: np.ndarray, ring, epsg: int | None) -> np.ndarray | None:
    """
    Boolean mask of grid cells inside the AOI ring.

    Works directly when the grid is geographic. For a projected grid it needs
    pyproj to convert the ring into grid coordinates.
    """
    ring_xy = ring
    if epsg and epsg != 4326:
        try:
            from pyproj import Transformer
            tf = Transformer.from_crs(4326, epsg, always_xy=True)
            ring_xy = [tf.transform(lon, lat) for lon, lat in ring]
        except ImportError:
            logger.warning("Grid is EPSG:%s and pyproj is missing - skipping AOI clip.", epsg)
            return None

    xmin = min(p[0] for p in ring_xy); xmax = max(p[0] for p in ring_xy)
    ymin = min(p[1] for p in ring_xy); ymax = max(p[1] for p in ring_xy)

    gx, gy = np.meshgrid(xs, ys)
    mask = (gx >= xmin) & (gx <= xmax) & (gy >= ymin) & (gy <= ymax)
    idx = np.argwhere(mask)
    for i, j in idx:
        if not point_in_ring(float(gx[i, j]), float(gy[i, j]), ring_xy):
            mask[i, j] = False
    return mask


# ---------------------------------------------------------------------------
# Core read
# ---------------------------------------------------------------------------
def read_gunw(
    path: Path,
    coherence_threshold: float = DEFAULT_COHERENCE,
    ref_lat: float | None = None,
    ref_lon: float | None = None,
    ref_radius_px: int = 5,
    clip_aoi: bool = True,
    flip_sign: bool = False,
) -> dict:
    with h5py.File(path, "r") as f:
        datasets = walk(f)

        p_phase = find(datasets, "unwrappedphase", ndim=2)
        p_coh = find(datasets, "coherencemagnitude", ndim=2) or find(datasets, "coherence", ndim=2)
        p_conn = find(datasets, "connectedcomponents", ndim=2)
        p_mask = find(datasets, "layovershadowmask", ndim=2) or find(datasets, "mask", ndim=2)
        p_x = find(datasets, "xcoordinates", ndim=1)
        p_y = find(datasets, "ycoordinates", ndim=1)

        if not p_phase:
            raise RuntimeError(
                "No 'unwrappedPhase' dataset found. Run --inspect and check the layout."
            )

        logger.info("phase      : %s %s", p_phase, datasets[p_phase][0])
        logger.info("coherence  : %s", p_coh or "NOT FOUND")
        logger.info("components : %s", p_conn or "not found")

        phase = np.asarray(f[p_phase][()], dtype=np.float64)
        coherence = np.asarray(f[p_coh][()], dtype=np.float64) if p_coh else None
        components = np.asarray(f[p_conn][()]) if p_conn else None
        losmask = np.asarray(f[p_mask][()]) if p_mask else None
        xs = np.asarray(f[p_x][()], dtype=np.float64) if p_x else None
        ys = np.asarray(f[p_y][()], dtype=np.float64) if p_y else None

        epsg = None
        p_proj = find(datasets, "projection")
        if p_proj:
            try:
                val = f[p_proj]
                epsg = int(val.attrs.get("epsg_code", np.ravel(val[()])[0]))
            except Exception:
                pass

        wavelength = read_wavelength(f, datasets)

    # phase -> LOS displacement, millimetres
    scale = -(wavelength / (4.0 * math.pi)) * 1000.0
    if flip_sign:
        scale = -scale
    disp = phase * scale

    valid = np.isfinite(disp) & (phase != 0)
    gates = {"finite": int(valid.sum())}

    if coherence is not None:
        valid &= coherence >= coherence_threshold
        gates[f"coherence>={coherence_threshold}"] = int(valid.sum())
    if components is not None:
        valid &= components > 0          # component 0 = not reliably unwrapped
        gates["connected_component>0"] = int(valid.sum())
    if losmask is not None and losmask.shape == valid.shape:
        try:
            valid &= losmask == 0        # 0 = good in the usual convention
            gates["layover/shadow clear"] = int(valid.sum())
        except Exception:
            pass

    # Quality mask WITHOUT the AOI clip. The reference point is usually chosen
    # on stable ground outside the area of interest, so it must not be clipped
    # away before we can use it.
    quality_valid = valid.copy()

    aoi_mask = None
    if clip_aoi and xs is not None and ys is not None:
        aoi_mask = build_aoi_mask(xs, ys, AOI_RING, epsg)
        if aoi_mask is not None:
            valid &= aoi_mask
            gates["inside AOI"] = int(valid.sum())

    # reference correction - unwrapped phase is relative
    ref_value = None
    if ref_lat is not None and ref_lon is not None and xs is not None and ys is not None:
        rx, ry = ref_lon, ref_lat
        if epsg and epsg != 4326:
            try:
                from pyproj import Transformer
                rx, ry = Transformer.from_crs(4326, epsg, always_xy=True).transform(ref_lon, ref_lat)
            except ImportError:
                logger.warning("pyproj missing - cannot place reference point on a projected grid.")
                rx = ry = None
        if rx is not None:
            j = int(np.argmin(np.abs(xs - rx)))
            i = int(np.argmin(np.abs(ys - ry)))
            r = ref_radius_px
            window = disp[max(0, i - r): i + r + 1, max(0, j - r): j + r + 1]
            wvalid = quality_valid[max(0, i - r): i + r + 1, max(0, j - r): j + r + 1]
            if wvalid.sum() > 0:
                ref_value = float(np.median(window[wvalid]))
                disp = disp - ref_value
                logger.info("Referenced to %.4f N %.4f E: subtracted %.2f mm (%d px)",
                            ref_lat, ref_lon, ref_value, int(wvalid.sum()))
            else:
                logger.warning(
                    "Reference window at %.4f N %.4f E has no usable pixels - NOT referenced. "
                    "Pick a point with good coherence, or widen --ref-radius.",
                    ref_lat, ref_lon)
    else:
        logger.warning(
            "No reference point given. Values are RELATIVE - differences within "
            "the scene are meaningful, absolute magnitudes are not."
        )

    return {
        "file": path.name,
        "reference_date": read_pair_dates(path)[0],
        "secondary_date": read_pair_dates(path)[1],
        "displacement_mm": disp,
        "coherence": coherence,
        "valid": valid,
        "xs": xs, "ys": ys, "epsg": epsg,
        "wavelength_m": wavelength,
        "ref_value_mm": ref_value,
        "gates": gates,
    }


# ---------------------------------------------------------------------------
# Reporting / output
# ---------------------------------------------------------------------------
def report(result: dict) -> dict:
    disp, valid = result["displacement_mm"], result["valid"]
    total = disp.size
    print(f"\n=== {result['file']}")
    print(f"  pair            {result['reference_date']} -> {result['secondary_date']}")
    print(f"  grid            {disp.shape[0]} x {disp.shape[1]} = {total:,} px")
    print(f"  wavelength      {result['wavelength_m']:.4f} m")
    if result["epsg"]:
        print(f"  projection      EPSG:{result['epsg']}")

    print("\n  quality gates (cumulative):")
    for name, n in result["gates"].items():
        print(f"    {name:<28} {n:>10,}  {100*n/total:>5.1f}%")

    n = int(valid.sum())
    stats = {"file": result["file"], "reference": result["reference_date"],
             "secondary": result["secondary_date"], "valid_px": n,
             "valid_pct": round(100 * n / total, 2)}

    if n == 0:
        print("\n  NO VALID PIXELS. Lower --coh-threshold or check the AOI overlaps the frame.")
        return stats

    vals = disp[valid]
    qs = np.percentile(vals, [1, 5, 25, 50, 75, 95, 99])
    print(f"\n  LOS displacement over {n:,} valid px (mm, +ve = away from satellite):")
    for label, q in zip(["p1", "p5", "p25", "median", "p75", "p95", "p99"], qs):
        print(f"    {label:<8}{q:>10.2f}")
    print(f"    {'min':<8}{vals.min():>10.2f}")
    print(f"    {'max':<8}{vals.max():>10.2f}")
    print(f"    {'std':<8}{vals.std():>10.2f}")

    if result["coherence"] is not None:
        print(f"\n  mean coherence (valid px): {result['coherence'][valid].mean():.3f}")

    if stats["valid_pct"] < 10:
        print("\n  WARNING: under 10% of the scene survived gating. Treat with suspicion -")
        print("  likely monsoon decorrelation or a 24-day span. Consider GOFF instead.")

    stats.update({k: round(float(v), 3) for k, v in
                  zip(["p1", "p5", "p25", "median", "p75", "p95", "p99"], qs)})
    stats["min"] = round(float(vals.min()), 3)
    stats["max"] = round(float(vals.max()), 3)
    stats["std"] = round(float(vals.std()), 3)
    stats["mean_coherence"] = (round(float(result["coherence"][valid].mean()), 4)
                               if result["coherence"] is not None else None)
    return stats


def write_geotiff(result: dict, out: Path) -> None:
    try:
        import rasterio
        from rasterio.transform import from_origin
    except ImportError:
        logger.error("rasterio is required for GeoTIFF output:  pip install rasterio")
        return
    xs, ys = result["xs"], result["ys"]
    if xs is None or ys is None:
        logger.error("No coordinate arrays - cannot georeference.")
        return
    data = np.where(result["valid"], result["displacement_mm"], np.nan).astype("float32")
    resx = abs(float(xs[1] - xs[0])); resy = abs(float(ys[1] - ys[0]))
    transform = from_origin(float(xs.min()) - resx / 2, float(ys.max()) + resy / 2, resx, resy)
    if ys[0] < ys[-1]:
        data = np.flipud(data)
    with rasterio.open(
        out, "w", driver="GTiff", height=data.shape[0], width=data.shape[1],
        count=1, dtype="float32", crs=f"EPSG:{result['epsg'] or 4326}",
        transform=transform, nodata=np.nan, compress="deflate",
    ) as dst:
        dst.write(data, 1)
        dst.update_tags(units="mm", pair=f"{result['reference_date']}_{result['secondary_date']}")
    logger.info("Wrote %s", out)


def write_quicklook(result: dict, out: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.error("matplotlib is required for quicklooks:  pip install matplotlib")
        return
    data = np.where(result["valid"], result["displacement_mm"], np.nan)
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        logger.error("Nothing valid to plot.")
        return
    lim = float(np.percentile(np.abs(finite), 98)) or 1.0
    fig, ax = plt.subplots(figsize=(9, 7))
    im = ax.imshow(data, cmap="RdBu_r", vmin=-lim, vmax=lim)
    ax.set_title(f"{result['reference_date']} to {result['secondary_date']}  |  LOS mm")
    ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(im, ax=ax, label="LOS displacement (mm), +ve away from satellite")
    fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)
    logger.info("Wrote %s", out)


# ---------------------------------------------------------------------------
def inspect(path: Path) -> None:
    with h5py.File(path, "r") as f:
        datasets = walk(f)
    print(f"\n{path.name}: {len(datasets)} datasets\n")
    interesting = ("unwrapped", "coherence", "connected", "mask", "coordinate",
                   "projection", "frequency", "incidence", "offset", "ionosphere")
    print("--- likely relevant ---")
    for p, (shape, dt) in sorted(datasets.items()):
        if any(k in p.lower() for k in interesting):
            print(f"  {p}\n      shape={shape} dtype={dt}")
    print("\n--- 2-D grids ---")
    for p, (shape, dt) in sorted(datasets.items()):
        if len(shape) == 2 and shape[0] > 50 and shape[1] > 50:
            print(f"  {p}  {shape} {dt}")


def main() -> int:
    ap = argparse.ArgumentParser(description="NISAR GUNW -> LOS displacement")
    ap.add_argument("--inspect", metavar="FILE", help="dump the HDF5 layout and exit")
    ap.add_argument("--read", metavar="FILE", help="read one GUNW")
    ap.add_argument("--batch", metavar="DIR", help="read every *GUNW*.h5 in a directory")
    ap.add_argument("--coh-threshold", type=float, default=DEFAULT_COHERENCE)
    ap.add_argument("--ref-lat", type=float, help="reference point latitude (stable ground)")
    ap.add_argument("--ref-lon", type=float, help="reference point longitude")
    ap.add_argument("--ref-radius", type=int, default=5, help="reference window half-width in pixels")
    ap.add_argument("--no-clip", action="store_true", help="do not clip to the AOI ring")
    ap.add_argument("--flip-sign", action="store_true", help="invert the LOS sign convention")
    ap.add_argument("--geotiff", metavar="OUT.tif")
    ap.add_argument("--quicklook", metavar="OUT.png")
    ap.add_argument("--csv", metavar="OUT.csv", help="write per-pair statistics")
    args = ap.parse_args()

    if args.inspect:
        inspect(Path(args.inspect)); return 0

    files: list[Path] = []
    if args.read:
        files = [Path(args.read)]
    elif args.batch:
        files = sorted(p for p in Path(args.batch).glob("*.h5") if "GUNW" in p.name.upper())
        if not files:
            logger.error("No *GUNW*.h5 files in %s", args.batch); return 1
        logger.info("Found %d GUNW files", len(files))
    else:
        ap.print_help(); return 0

    rows = []
    for fp in files:
        try:
            result = read_gunw(
                fp,
                coherence_threshold=args.coh_threshold,
                ref_lat=args.ref_lat, ref_lon=args.ref_lon, ref_radius_px=args.ref_radius,
                clip_aoi=not args.no_clip, flip_sign=args.flip_sign,
            )
        except Exception as exc:
            logger.error("%s: %s", fp.name, exc)
            continue
        rows.append(report(result))
        if args.geotiff and len(files) == 1:
            write_geotiff(result, Path(args.geotiff))
        if args.quicklook and len(files) == 1:
            write_quicklook(result, Path(args.quicklook))

    if args.csv and rows:
        keys = sorted({k for r in rows for k in r})
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
        logger.info("Wrote %s (%d pairs)", args.csv, len(rows))

    return 0


if __name__ == "__main__":
    sys.exit(main())
