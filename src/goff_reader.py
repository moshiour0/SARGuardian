"""
goff_reader.py
--------------
Read NISAR L2 GOFF (geocoded pixel offsets) into ground displacement.

Why GOFF and not GUNW
=====================
Interferometric phase is ambiguous once displacement between passes exceeds
lambda/4 - about 61 mm per 12-day pair at NISAR's 0.2439 m wavelength, i.e.
5 mm/day. Blatten was moving at 650 mm/day six days before it failed. Phase
cannot measure a slope that is actually failing; offset tracking can, because
it cross-correlates image patches and never wraps.

The trade is precision. GUNW resolves millimetres, GOFF resolves a fraction of
a resolution cell. Between the phase ceiling and the offset noise floor there
is a band of velocities neither product sees - measure that floor on the winter
pairs with --noise-floor rather than assuming it.

What is in the product
======================
    slantRangeOffset        metres, along the line of sight
    alongTrackOffset        metres, along the flight direction (~N-S)
    slantRangeOffsetVariance, alongTrackOffsetVariance   metres^2
    correlationSurfacePeak  0-1, the matching quality
    snr                     offset signal-to-noise

Two layers at different correlation window sizes. layer1 is the finer window
(sharper, noisier); layer2 the coarser (smoother, more robust). Compare them:
a signal present in both is far more credible than one in either alone.

NASA labels these "raw (unculled, unfiltered)", so gating and outlier rejection
are the caller's job. This module does both and says how much it removed.

Usage
-----
    python src/goff_reader.py --inspect FILE.h5
    python src/goff_reader.py --read FILE.h5 --auto-ref --quicklook out.png
    python src/goff_reader.py --batch data/nisar_l2/GOFF --csv outputs/goff.csv
    python src/goff_reader.py --noise-floor data/nisar_l2/GOFF/2025-12_winter

Requires: h5py, numpy.  Optional: pyproj (AOI clip on projected grids),
matplotlib (quicklooks).
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from pathlib import Path

import numpy as np

try:
    import h5py
except ImportError:
    print("h5py is required:  pip install h5py", file=sys.stderr)
    raise SystemExit(1)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gunw_reader import AOI_RING, build_aoi_mask, walk, read_pair_dates  # noqa: E402

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("goff")

DEFAULT_CORRELATION = 0.3
DEFAULT_SNR = 3.0
MAD_SIGMA = 1.4826          # MAD -> Gaussian sigma


# ---------------------------------------------------------------------------
def layer_groups(datasets: dict) -> list[str]:
    """Every .../pixelOffsets/<pol>/<layerN> group present, in order."""
    groups = set()
    for path in datasets:
        m = re.match(r"(.*/pixelOffsets/[^/]+/layer\d+)/", path)
        if m:
            groups.add(m.group(1))
    return sorted(groups)


def robust_sigma(a: np.ndarray) -> float:
    """MAD-based scatter - immune to the gross outliers in unculled offsets."""
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan")
    return float(MAD_SIGMA * np.median(np.abs(a - np.median(a))))


def read_goff(
    path: Path,
    layer: str | int = "best",
    correlation_min: float = DEFAULT_CORRELATION,
    snr_min: float = DEFAULT_SNR,
    cull_sigma: float = 4.0,
    clip_aoi: bool = True,
    auto_ref: bool = True,
    deramp: bool = True,
) -> dict:
    with h5py.File(path, "r") as f:
        datasets = walk(f)
        groups = layer_groups(datasets)
        if not groups:
            raise RuntimeError("No pixelOffsets layer groups found. Run --inspect.")

        if layer == "best":
            chosen = groups              # read all, pick later
        elif isinstance(layer, int) or str(layer).isdigit():
            want = f"layer{int(layer)}"
            chosen = [g for g in groups if g.endswith("/" + want)] or groups[:1]
        else:
            chosen = [g for g in groups if g.endswith("/" + str(layer))] or groups[:1]

        results = {}
        for grp in chosen:
            def ds(name):
                key = f"{grp}/{name}"
                return np.asarray(f[key][()], dtype=np.float64) if key in datasets else None

            rng = ds("slantRangeOffset")
            azi = ds("alongTrackOffset")
            corr = ds("correlationSurfacePeak")
            snr = ds("snr")
            rng_var = ds("slantRangeOffsetVariance")
            if rng is None:
                continue

            xs = np.asarray(f[f"{grp}/xCoordinates"][()], dtype=np.float64)
            ys = np.asarray(f[f"{grp}/yCoordinates"][()], dtype=np.float64)
            epsg = None
            if f"{grp}/projection" in datasets:
                try:
                    v = f[f"{grp}/projection"]
                    epsg = int(v.attrs.get("epsg_code", np.ravel(v[()])[0]))
                except Exception:
                    pass
            # Key by polarisation AND layer. A quad-pol product carries
            # HH/layer1..3 and VV/layer1..3; keying on the layer name alone
            # lets VV silently overwrite HH and the read collapses to nothing.
            pol_name = grp.split("/")[-2]
            results[f"{pol_name}/{grp.rsplit('/', 1)[-1]}"] = dict(
                rng=rng, azi=azi, corr=corr, snr=snr, rng_var=rng_var,
                xs=xs, ys=ys, epsg=epsg, group=grp,
            )

    if not results:
        raise RuntimeError("No slantRangeOffset found in any layer.")

    out = {}
    for name, L in results.items():
        rng, azi = L["rng"], L["azi"]
        gates = {"finite": int(np.isfinite(rng).sum())}
        valid = np.isfinite(rng)
        if azi is not None:
            valid &= np.isfinite(azi)

        if L["corr"] is not None:
            valid &= L["corr"] >= correlation_min
            gates[f"correlation>={correlation_min}"] = int(valid.sum())
        if L["snr"] is not None:
            valid &= L["snr"] >= snr_min
            gates[f"snr>={snr_min}"] = int(valid.sum())

        quality_valid = valid.copy()

        aoi_mask = None
        if clip_aoi:
            aoi_mask = build_aoi_mask(L["xs"], L["ys"], AOI_RING, L["epsg"])
            if aoi_mask is not None:
                valid &= aoi_mask
                gates["inside AOI"] = int(valid.sum())

        # Reference: offsets carry a bias from orbit and coregistration error,
        # exactly as unwrapped phase carries an arbitrary constant.
        ref_rng = ref_azi = None
        if auto_ref:
            blk = 11
            ny, nx = rng.shape
            by, bx = ny // blk, nx // blk
            if by and bx:
                fullv = quality_valid[:by * blk, :bx * blk].reshape(by, blk, bx, blk).all(axis=(1, 3))
                if aoi_mask is not None:
                    ina = aoi_mask[:by * blk, :bx * blk].reshape(by, blk, bx, blk).any(axis=(1, 3))
                    fullv &= ~ina
                if L["corr"] is not None:
                    sc = L["corr"][:by * blk, :bx * blk].reshape(by, blk, bx, blk).mean(axis=(1, 3))
                    score = np.where(fullv, sc, -1.0)
                else:
                    score = fullv.astype(float) - 1.0 * (~fullv)
                if score.max() > 0:
                    bi, bj = np.unravel_index(int(np.argmax(score)), score.shape)
                    i, j = bi * blk + blk // 2, bj * blk + blk // 2
                    w = quality_valid[i - 5:i + 6, j - 5:j + 6]
                    if w.sum() > 0:
                        ref_rng = float(np.median(rng[i - 5:i + 6, j - 5:j + 6][w]))
                        rng = rng - ref_rng
                        if azi is not None:
                            ref_azi = float(np.median(azi[i - 5:i + 6, j - 5:j + 6][w]))
                            azi = azi - ref_azi
                        logger.info("[%s] auto reference: corr %.2f, subtracted "
                                    "range %+.3f m, azimuth %+.3f m",
                                    name, float(score[bi, bj]), ref_rng, ref_azi or 0.0)

        # Deramp. Residual orbit and coregistration error leaves a near-planar
        # tilt across an offset field, so subtracting a single constant is not
        # enough - a stable winter pair still showed a -646 mm median before
        # this. Fit a + bx + cy on good ground OUTSIDE the AOI and remove it,
        # so the plane is never fitted to the signal being measured.
        ramp_rms = None
        if deramp:
            fit_mask = quality_valid.copy()
            if aoi_mask is not None:
                fit_mask &= ~aoi_mask
            if fit_mask.sum() > 5000:
                gy, gx = np.nonzero(fit_mask)
                step = max(1, gy.size // 200_000)          # cap the design matrix
                gy, gx = gy[::step], gx[::step]
                zr = rng[gy, gx]
                keep = np.abs(zr - np.median(zr)) <= 4 * robust_sigma(zr)
                gy, gx, zr = gy[keep], gx[keep], zr[keep]
                A = np.column_stack([np.ones_like(gx, dtype=float), gx.astype(float),
                                     gy.astype(float)])
                coef, *_ = np.linalg.lstsq(A, zr, rcond=None)
                YY, XX = np.mgrid[0:rng.shape[0], 0:rng.shape[1]]
                plane = coef[0] + coef[1] * XX + coef[2] * YY
                ramp_rms = float(np.sqrt(np.mean((plane[fit_mask] * 1000.0) ** 2)))
                rng = rng - plane
                if azi is not None:
                    za = azi[gy, gx]
                    coefa, *_ = np.linalg.lstsq(A, za, rcond=None)
                    azi = azi - (coefa[0] + coefa[1] * XX + coefa[2] * YY)
                logger.info("[%s] deramped: removed plane of %.0f mm RMS "
                            "(fitted on %d px outside the AOI)",
                            name, ramp_rms, int(gy.size))
            else:
                logger.warning("[%s] too little stable ground to deramp", name)

        # Cull gross outliers. NASA ships these unculled and a raw offset field
        # always contains mismatches that are metres wrong.
        if valid.sum() > 10 and cull_sigma > 0:
            s = robust_sigma(rng[valid])
            med = float(np.median(rng[valid]))
            if np.isfinite(s) and s > 0:
                keep = np.abs(rng - med) <= cull_sigma * s
                valid &= keep
                gates[f"within {cull_sigma:g}x MAD"] = int(valid.sum())

        out[name] = {
            "layer": name,
            "range_m": rng,
            "azimuth_m": azi,
            "correlation": L["corr"],
            "snr": L["snr"],
            "range_var": L["rng_var"],
            "valid": valid,
            "xs": L["xs"], "ys": L["ys"], "epsg": L["epsg"],
            "gates": gates,
            "ramp_rms_mm": ramp_rms,
            "ref_range_m": ref_rng,
            "ref_azimuth_m": ref_azi,
        }

    ref, sec = read_pair_dates(path)
    return {"file": path.name, "reference_date": ref, "secondary_date": sec,
            "layers": out}


# ---------------------------------------------------------------------------
def report(res: dict, span_days: float | None = None) -> list[dict]:
    print(f"\n=== {res['file']}")
    print(f"  pair  {res['reference_date']} -> {res['secondary_date']}"
          + (f"   ({span_days:.0f} d)" if span_days else ""))
    rows = []
    for name, L in sorted(res["layers"].items()):
        v = L["valid"]
        total = L["range_m"].size
        n = int(v.sum())
        print(f"\n  --- {name} ---")
        for k, c in L["gates"].items():
            print(f"    {k:<26}{c:>12,}  {100*c/total:>5.1f}%")
        if n == 0:
            print("    NO VALID PIXELS")
            continue

        r = L["range_m"][v] * 1000.0            # mm
        sig_r = robust_sigma(r)
        row = {"file": res["file"], "layer": name,
               "reference": res["reference_date"], "secondary": res["secondary_date"],
               "valid_px": n, "valid_pct": round(100 * n / total, 2),
               "range_median_mm": round(float(np.median(r)), 2),
               "range_mad_sigma_mm": round(sig_r, 2),
               "range_p5_mm": round(float(np.percentile(r, 5)), 2),
               "range_p95_mm": round(float(np.percentile(r, 95)), 2),
               "range_max_abs_mm": round(float(np.max(np.abs(r))), 2)}
        print(f"    range offset  median {row['range_median_mm']:>10.1f} mm   "
              f"robust sigma {sig_r:>8.1f} mm")
        print(f"                  p5..p95 {row['range_p5_mm']:>9.1f} .. "
              f"{row['range_p95_mm']:.1f} mm   |max| {row['range_max_abs_mm']:.1f} mm")

        if L["azimuth_m"] is not None:
            a = L["azimuth_m"][v] * 1000.0
            row["azimuth_median_mm"] = round(float(np.median(a)), 2)
            row["azimuth_mad_sigma_mm"] = round(robust_sigma(a), 2)
            mag = np.hypot(r, a)
            row["magnitude_p95_mm"] = round(float(np.percentile(mag, 95)), 2)
            print(f"    azimuth       median {row['azimuth_median_mm']:>10.1f} mm   "
                  f"robust sigma {row['azimuth_mad_sigma_mm']:>8.1f} mm")
            print(f"    magnitude p95 {row['magnitude_p95_mm']:>10.1f} mm")

        if L["correlation"] is not None:
            row["mean_correlation"] = round(float(L["correlation"][v].mean()), 3)
            print(f"    mean correlation {row['mean_correlation']:.3f}")

        if span_days:
            row["span_days"] = span_days
            row["range_velocity_mm_day"] = round(row["range_median_mm"] / span_days, 3)
            row["detect_floor_mm_day"] = round(3 * sig_r / span_days, 2)
            print(f"    implied range velocity {row['range_velocity_mm_day']:+.2f} mm/day")
            print(f"    3-sigma detection floor {row['detect_floor_mm_day']:.1f} mm/day")
        rows.append(row)
    return rows


def span_of(name: str) -> float | None:
    from datetime import datetime
    s = re.findall(r"_(\d{8})T\d{6}", name)
    if len(s) < 4:
        return None
    a = datetime.strptime(s[0], "%Y%m%d")
    b = datetime.strptime(s[2], "%Y%m%d")
    return (b - a).days


def noise_floor(directory: Path, **kw) -> None:
    """
    Measure the offset-tracking noise floor on pairs you believe are stable.

    On a static winter glacier the scatter in the referenced offset field IS
    the detection threshold. This turns an assumed noise level into a measured
    one, and it is the number that defines the blind band between the GUNW
    phase ceiling and GOFF sensitivity.
    """
    files = sorted(p for p in directory.rglob("*.h5") if "GOFF" in p.name.upper())
    if not files:
        logger.error("No GOFF files under %s", directory)
        return
    logger.info("Measuring noise floor on %d stable pairs", len(files))

    all_rows = []
    for fp in files:
        try:
            res = read_goff(fp, **kw)
        except Exception as exc:
            logger.error("%s: %s", fp.name, exc)
            continue
        all_rows += report(res, span_days=span_of(fp.name))

    if not all_rows:
        return
    print("\n" + "=" * 78)
    print("OFFSET-TRACKING NOISE FLOOR")
    print("=" * 78)
    for lay in sorted({r["layer"] for r in all_rows}):
        sub = [r for r in all_rows if r["layer"] == lay]
        sig = np.array([r["range_mad_sigma_mm"] for r in sub])
        floors = np.array([r.get("detect_floor_mm_day", np.nan) for r in sub], dtype=float)
        print(f"\n  {lay}: {len(sub)} pairs")
        print(f"    range sigma      median {np.median(sig):>8.1f} mm   "
              f"range {sig.min():.1f} - {sig.max():.1f} mm")
        good = floors[np.isfinite(floors)]
        if good.size:
            print(f"    3-sigma floor    median {np.median(good):>8.1f} mm/day  "
                  f"range {good.min():.1f} - {good.max():.1f}")
            print(f"    -> GOFF cannot see motion slower than ~{np.median(good):.0f} mm/day "
                  f"on a 12-day pair")
    print("\n  GUNW phase ceiling for comparison: 61 mm per 12-day pair = 5.1 mm/day")
    print("  Anything between those two numbers is invisible to both products.")


def quicklook(res: dict, out: Path, layer: str | None = None) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.error("matplotlib required for --quicklook")
        return
    name = layer or sorted(res["layers"])[0]
    L = res["layers"][name]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, key, title in ((axes[0], "range_m", "Slant range offset"),
                           (axes[1], "azimuth_m", "Along-track offset")):
        if L[key] is None:
            ax.axis("off"); continue
        d = np.where(L["valid"], L[key] * 1000.0, np.nan)
        fin = d[np.isfinite(d)]
        lim = float(np.percentile(np.abs(fin), 98)) if fin.size else 1.0
        im = ax.imshow(d, cmap="RdBu_r", vmin=-lim, vmax=lim)
        ax.set_title(f"{title} (mm)"); ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=.046)
    fig.suptitle(f"{res['reference_date']} to {res['secondary_date']}  |  GOFF {name}")
    fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)
    logger.info("Wrote %s", out)


def inspect(path: Path) -> None:
    with h5py.File(path, "r") as f:
        datasets = walk(f)
        groups = layer_groups(datasets)
        print(f"\n{path.name}")
        print(f"  {len(datasets)} datasets, {len(groups)} offset layer(s)")
        for g in groups:
            print(f"\n  {g}")
            for k, (shp, dt) in sorted(datasets.items()):
                if k.startswith(g + "/"):
                    u = f[k].attrs.get("units")
                    u = u.decode() if isinstance(u, bytes) else u
                    print(f"    {k.rsplit('/',1)[-1]:<28}{str(shp):<16}{dt}  {u or ''}")


def main() -> int:
    ap = argparse.ArgumentParser(description="NISAR GOFF -> ground displacement")
    m = ap.add_mutually_exclusive_group(required=True)
    m.add_argument("--inspect", metavar="FILE")
    m.add_argument("--read", metavar="FILE")
    m.add_argument("--batch", metavar="DIR")
    m.add_argument("--noise-floor", metavar="DIR",
                   help="measure the detection floor on pairs you believe are stable")
    ap.add_argument("--layer", default="best", help="layer1, layer2, or best (both)")
    ap.add_argument("--correlation", type=float, default=DEFAULT_CORRELATION)
    ap.add_argument("--snr", type=float, default=DEFAULT_SNR)
    ap.add_argument("--cull-sigma", type=float, default=4.0)
    ap.add_argument("--no-clip", action="store_true")
    ap.add_argument("--no-auto-ref", action="store_true")
    ap.add_argument("--no-deramp", action="store_true",
                    help="keep the planar orbital ramp (default: remove it)")
    ap.add_argument("--quicklook", metavar="OUT.png")
    ap.add_argument("--csv", metavar="OUT.csv")
    args = ap.parse_args()

    kw = dict(layer=args.layer, correlation_min=args.correlation, snr_min=args.snr,
              cull_sigma=args.cull_sigma, clip_aoi=not args.no_clip,
              auto_ref=not args.no_auto_ref, deramp=not args.no_deramp)

    if args.inspect:
        inspect(Path(args.inspect)); return 0

    if args.noise_floor:
        noise_floor(Path(args.noise_floor), **kw); return 0

    files = [Path(args.read)] if args.read else sorted(
        p for p in Path(args.batch).rglob("*.h5") if "GOFF" in p.name.upper())
    if not files:
        logger.error("No GOFF files found"); return 1

    rows = []
    for fp in files:
        try:
            res = read_goff(fp, **kw)
        except Exception as exc:
            logger.error("%s: %s", fp.name, exc); continue
        rows += report(res, span_days=span_of(fp.name))
        if args.quicklook and len(files) == 1:
            quicklook(res, Path(args.quicklook))

    if args.csv and rows:
        keys = sorted({k for r in rows for k in r})
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
        logger.info("Wrote %s (%d rows)", args.csv, len(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
