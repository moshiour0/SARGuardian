"""
pixel_stack.py
--------------
Per-pixel displacement from the exported AOI rasters.

Why not just the AOI median
===========================
Every number in this project so far is a median over the whole AOI. That is
the right statistic for asking whether a valley is moving, and the wrong one
for finding a landslide inside it: a 200 m scar creeping at 3 mm/day inside a
2 km box is a few hundred pixels among twelve thousand, and the median never
notices. Only a per-pixel fit can see it.

Removing the atmosphere without a weather model
===============================================
Two separate AOIs in the same scene were measured independently and their
winter displacements correlate at +0.974 - 82% of the variance is shared.
Different ground kilometres apart cannot move identically, so that shared part
is the troposphere, and the residual after removing it is 5.3 mm against
12.6 mm before: a 2.4-fold improvement from arithmetic alone.

The same logic applies spatially. Subtracting the scene median epoch by epoch
removes whatever is common to the whole AOI and keeps whatever is local. It
cannot remove a turbulent screen with structure at landslide scale, and it
cannot see motion that fills the AOI uniformly - both limits are stated in the
output rather than hidden.

What "significant" means here
=============================
A per-pixel velocity is compared with the scatter of its own fit, not with a
global threshold, and then against the distribution of every other pixel. A
handful of significant pixels scattered at random is noise finding its tail;
a significant CLUSTER is a candidate. The tool reports both counts and refuses
to call either one a detection on its own.

Usage
-----
    python src/pixel_stack.py --dir outputs/export --season winter
    python src/pixel_stack.py --dir outputs/export --season winter --geotiff out.tif
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("stack")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import resolve  # noqa: E402

WINTER_MONTHS = (11, 12, 1, 2, 3)


def parse_export(path: Path):
    """
    Parse dates and the complete measurement-network identity from the source tag.

    Path/look alone is not a sufficient time-series key. Frame, subswath,
    bandwidth, polarization and processing CRID can change the effective
    observation geometry or processing lineage. Combining pair intervals from
    different keys can manufacture a temporal series from incompatible inputs.
    """
    m = re.match(r"GUNW_(\d{8})_(\d{8})(?:_([A-Z]{2}))?\.tif$", path.name)
    if not m:
        return None
    import rasterio

    try:
        with rasterio.open(path) as src:
            source_name = str(src.tags().get("source") or "")
    except Exception as exc:
        logger.warning("Cannot read source identity from %s: %s", path.name, type(exc).__name__)
        return None

    granule = re.search(
        r"_GUNW_(?P<version>\d+)_(?P<path>\d+)_(?P<look>[AD])_"
        r"(?P<frame>\d+)_(?P<subswath>\d+)_(?P<bandwidth>\d+)_"
        r"(?P<polarization>[A-Z]{2})_.*?_(?P<crid>[A-Z]\d{5})_[A-Z]_",
        source_name,
    )
    if not granule:
        logger.warning("Cannot derive a complete GUNW network key from %s", path.name)
        return None

    geom = f"{granule.group('look')}{int(granule.group('path')):03d}"
    network_id = (
        f"{geom}-F{int(granule.group('frame')):03d}"
        f"-SW{int(granule.group('subswath')):03d}"
        f"-BW{int(granule.group('bandwidth')):04d}"
        f"-{granule.group('polarization')}-{granule.group('crid')}"
    )
    return {
        "path": path,
        "ref": datetime.strptime(m.group(1), "%Y%m%d").date(),
        "sec": datetime.strptime(m.group(2), "%Y%m%d").date(),
        "proc": m.group(3) or "PR",
        "geom": geom,
        "network_id": network_id,
    }


def align(rasters: list[dict]) -> tuple[np.ndarray, np.ndarray, dict]:
    """
    Put every raster on one grid using its geotransform, not its shape.

    The exporter crops to each pair's own valid bounding box, so two rasters of
    the same AOI can differ by a column. Stacking them by array index silently
    shears the scene by one pixel per epoch - a displacement gradient that is
    pure bookkeeping. Aligning on coordinates makes that impossible.
    """
    import rasterio

    metas = []
    for r in rasters:
        with rasterio.open(r["path"]) as src:
            transform = src.transform
            if src.crs is None:
                raise SystemExit(f"Raster has no CRS: {r['path']}")
            if abs(transform.b) > 1e-12 or abs(transform.d) > 1e-12:
                raise SystemExit(f"Rotated/sheared raster is unsupported by this north-up stack: {r['path']}")
            disp = src.read(1, masked=True).filled(np.nan).astype("float32")
            coh = src.read(2, masked=True).filled(np.nan).astype("float32") if src.count > 1 else None
            metas.append({"t": transform, "h": src.height, "w": src.width,
                          "crs": src.crs, "disp": disp, "coh": coh})

    res = {(round(m["t"].a, 9), round(m["t"].e, 9)) for m in metas}
    if len(res) != 1:
        raise SystemExit(f"Rasters have different resolutions: {res}")
    if len({str(m["crs"]) for m in metas}) != 1:
        raise SystemExit("Rasters are in different CRSs")

    px, py = metas[0]["t"].a, metas[0]["t"].e
    if not (np.isfinite(px) and np.isfinite(py) and px > 0 and py < 0):
        raise SystemExit(f"Expected a finite north-up raster grid with positive x and negative y posting; got {(px, py)}")
    if any(
        not np.isclose(m["t"].a, px, rtol=0, atol=1e-9)
        or not np.isclose(m["t"].e, py, rtol=0, atol=1e-9)
        for m in metas
    ):
        raise SystemExit("Rasters have incompatible x/y posting.")
    x0 = min(m["t"].c for m in metas)
    y0 = max(m["t"].f for m in metas)
    x1 = max(m["t"].c + m["w"] * px for m in metas)
    y1 = min(m["t"].f + m["h"] * py for m in metas)

    # Do not round a half-pixel or off-lattice offset into false spatial
    # agreement. Every source origin and the union extent must fall on the
    # same lattice before any array assignment.
    tolerance_px = 1e-6
    for m in metas:
        raw_col = (m["t"].c - x0) / px
        raw_row = (m["t"].f - y0) / py
        if abs(raw_col - round(raw_col)) > tolerance_px or abs(raw_row - round(raw_row)) > tolerance_px:
            raise SystemExit(
                f"Raster origin is not aligned to the common pixel lattice: {m['path'] if 'path' in m else m['t']}"
            )
    width_float = (x1 - x0) / px
    height_float = (y1 - y0) / py
    if abs(width_float - round(width_float)) > tolerance_px or abs(height_float - round(height_float)) > tolerance_px:
        raise SystemExit("Raster union extent is not an integer number of common-grid pixels.")
    W = int(round(width_float))
    H = int(round(height_float))
    if W <= 0 or H <= 0:
        raise SystemExit(f"Invalid common-grid dimensions {(H, W)}.")

    cube = np.full((len(metas), H, W), np.nan, dtype="float32")
    cohs = np.full((len(metas), H, W), np.nan, dtype="float32")
    for k, m in enumerate(metas):
        j = int(round((m["t"].c - x0) / px))
        i = int(round((m["t"].f - y0) / py))
        cube[k, i:i + m["h"], j:j + m["w"]] = m["disp"]
        if m["coh"] is not None:
            cohs[k, i:i + m["h"], j:j + m["w"]] = m["coh"]

    import rasterio.transform as rt
    return cube, cohs, {"transform": rt.Affine(px, 0, x0, 0, py, y0),
                        "crs": metas[0]["crs"], "shape": (H, W)}


def chain(rasters: list[dict]) -> list[list[dict]]:
    """
    Contiguous runs within one exact measurement network.

    ASC/DESC, frame, subswath, bandwidth, polarization and CRID are kept
    separate. A common path/look label alone does not prove that two GUNW
    rasters observe the same LOS projection or processing lineage.
    """
    runs = []
    network_ids = sorted({r.get("network_id") for r in rasters if r.get("network_id")})
    for network_id in network_ids:
        group = sorted(
            (r for r in rasters if r.get("network_id") == network_id),
            key=lambda r: (r["ref"], r["sec"]),
        )
        cur = [group[0]]
        for a, b in zip(group, group[1:]):
            if a["sec"] == b["ref"]:
                cur.append(b)
            else:
                runs.append(cur)
                cur = [b]
        runs.append(cur)
    return runs


def stack(run: list[dict], remove_common_mode: bool = True):
    if len(run) < 2:
        raise ValueError(
            "Per-pixel velocity requires at least two consecutive pair intervals "
            "(three acquisition epochs); one pair cannot estimate residual-based uncertainty."
        )
    cube, cohs, geo = align(run)
    n_ep = len(run) + 1
    H, W = geo["shape"]

    common_mode_pair_medians_mm: list[float | None] = []
    for k in range(cube.shape[0]):
        finite_values = cube[k][np.isfinite(cube[k])]
        med = float(np.median(finite_values)) if finite_values.size else float("nan")
        common_mode_pair_medians_mm.append(med if np.isfinite(med) else None)
        if remove_common_mode and np.isfinite(med):
            # Removing each pair's spatial median can also remove true broad
            # AOI-wide motion. Its exact applied value is retained for audit,
            # and the pipeline also runs a no-removal sensitivity diagnostic.
            cube[k] -= med

    cum = np.full((n_ep, H, W), np.nan, dtype="float32")
    # The first epoch is a relative zero only where at least one actual
    # pairwise measurement exists; it is not synthetic support for nodata.
    observed_pair_count = np.isfinite(cube).sum(axis=0).astype("uint16")
    # A missing pair breaks the per-pixel temporal chain. Keep both the total
    # number of observed pairs and the length of the longest contiguous segment
    # so the spatial attrition is visible instead of disappearing into the
    # complete-case velocity mask.
    longest_contiguous_pair_count = np.zeros((H, W), dtype="uint16")
    current_contiguous_pair_count = np.zeros((H, W), dtype="uint16")
    for pair_index in range(cube.shape[0]):
        pair_is_valid = np.isfinite(cube[pair_index])
        current_contiguous_pair_count = np.where(
            pair_is_valid, current_contiguous_pair_count + 1, 0
        ).astype("uint16")
        longest_contiguous_pair_count = np.maximum(
            longest_contiguous_pair_count, current_contiguous_pair_count
        )
    has_any_observation = observed_pair_count > 0
    cum[0, has_any_observation] = 0.0
    for k in range(cube.shape[0]):
        previous = cum[k]
        observed = np.isfinite(cube[k]) & np.isfinite(previous)
        cum[k + 1, observed] = previous[observed] + cube[k, observed]

    complete = np.isfinite(cum).all(axis=0)
    partial_support = (observed_pair_count > 0) & ~complete
    days = np.array([0.0] + list(np.cumsum([(r["sec"] - r["ref"]).days for r in run])))
    valid_pair_count = int(complete.sum())

    vel = np.full((H, W), np.nan, dtype="float32")
    sig = np.full((H, W), np.nan, dtype="float32")
    dof = n_ep - 2
    if dof <= 0:
        # Defensive invariant even though the public API gate above requires
        # at least three epochs.
        raise ValueError("Insufficient temporal degrees of freedom for a residual-based velocity uncertainty.")
    if complete.any():
        y = cum[:, complete]                            # (n_ep, n_px)
        A = np.vstack([days, np.ones_like(days)]).T
        if np.linalg.matrix_rank(A) < 2 or not np.isfinite(A).all():
            raise ValueError("Temporal design matrix is singular or non-finite.")
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        resid = A @ coef - y
        s = np.sqrt((resid ** 2).sum(axis=0) / dof)
        se = s * np.sqrt(np.linalg.inv(A.T @ A)[0, 0])
        vel[complete] = coef[0]
        sig[complete] = se
    return {
        "velocity": vel, "sigma": sig, "cumulative": cum, "coherence": cohs,
        "valid": complete, "days": days, "geo": geo, "run": run,
        "observed_pair_count": observed_pair_count,
        "longest_contiguous_pair_count": longest_contiguous_pair_count,
        "partial_support": partial_support,
        "complete_series_pixels": valid_pair_count,
        "total_grid_pixels": int(complete.size),
        "minimum_pairs_for_velocity": 2,
        "uncertainty_status": "RESIDUAL_BASED_DIAGNOSTIC_NOT_INDEPENDENTLY_VALIDATED" if dof > 0 else "WITHHELD_NO_DOF",
        "common_mode_removed": bool(remove_common_mode),
        "common_mode_policy": "PAIRWISE_MEDIAN_REMOVED" if remove_common_mode else "PAIRWISE_MEDIAN_RETAINED",
        "common_mode_pair_medians_mm": common_mode_pair_medians_mm,
    }


def report(res: dict, label: str) -> bool | None:
    v, s, ok = res["velocity"], res["sigma"], res["valid"]
    run = res["run"]
    print(f"\n{'='*72}\n{label}   {run[0]['ref']} .. {run[-1]['sec']}  "
          f"({len(run)} pairs, {int(res['days'][-1])} d)\n{'='*72}")
    n = int(ok.sum())
    total = int(ok.size)
    partial = int(((res["observed_pair_count"] > 0) & ~ok).sum())
    never = int((res["observed_pair_count"] == 0).sum())
    print(f"  common-mode policy            : {res['common_mode_policy']}")
    print(f"  residual uncertainty status   : {res['uncertainty_status']}")
    print(f"  pixels with a complete series : {n:,} of {total:,} "
          f"({100*n/max(total, 1):.1f}%)")
    print(f"  partial temporal support      : {partial:,} px "
          "(some observed pairs, but not a complete connected chain)")
    print(f"  no valid pair observations    : {never:,} px")
    if n == 0:
        res["spatial_null_test_status"] = "WITHHELD_NO_COMPLETE_SERIES"
        print("  nothing to fit")
        print("\n  VERDICT: NOT TESTABLE. No complete pixel time series is available.")
        return None

    vv, ss = v[ok], s[ok]
    rob = 1.4826 * np.median(np.abs(vv - np.median(vv)))
    print(f"  velocity, robust scatter      : {rob:.3f} mm/day "
          f"({rob*365:.0f} mm/yr)")
    print(f"  median per-pixel fit sigma    : {np.median(ss):.3f} mm/day")

    # Significant against its own fit, and against the population.
    own = np.zeros_like(v, bool); own[ok] = np.abs(vv) > 3 * ss
    pop = np.zeros_like(v, bool); pop[ok] = np.abs(vv - np.median(vv)) > 3 * rob
    both = own & pop
    print(f"\n  3-sigma against its own fit   : {int(own.sum()):,} px "
          f"({100*own.sum()/max(n,1):.2f}%)")
    print(f"  3-sigma against the population: {int(pop.sum()):,} px "
          f"({100*pop.sum()/max(n,1):.2f}%)")
    print(f"  both                          : {int(both.sum()):,} px")
    print(f"  expected by chance at 3-sigma : {0.0027*n:.0f} px")

    lab, count = _largest_cluster(both)
    print(f"\n  largest connected cluster     : {count} px", end="")
    if count:
        ii, jj = np.nonzero(lab)
        print(f"  (rows {ii.min()}-{ii.max()}, cols {jj.min()}-{jj.max()})")
        print(f"  cluster mean velocity         : "
              f"{float(np.nanmean(v[lab])):+.3f} mm/day")
        pos = float(np.nanmean(v[lab] > 0))
        print(f"  cluster sign consistency      : {100*max(pos, 1-pos):.0f}% "
              f"one direction")
    else:
        print()

    # --- the null the earlier version of this tool did not have --------------
    # Comparing against a 0.27% chance rate assumes pixels are independent.
    # Tropospheric noise is correlated over kilometres, so it CLUSTERS on its
    # own, and cluster size proves nothing until compared with what this
    # scene's own correlation length produces from pure noise.
    corr_km = _correlation_length_km(v, ok, res["geo"])
    if corr_km is None:
        res["spatial_null_test_status"] = "WITHHELD_INSUFFICIENT_SPATIAL_NEIGHBOR_SUPPORT"
        res["spatial_correlation_length_km"] = None
        print("\n  measured correlation length   : NOT ESTIMABLE")
        print(f"  valid complete pixels          : {n} (requires at least "
              f"{MIN_SPATIAL_CORRELATION_PAIRS} pixels and neighbouring pairs)")
        print("  noise-only cluster test        : WITHHELD")
        print("\n  VERDICT: NOT TESTABLE. Too few nearby valid pixels are available "
              "to estimate the spatial background-noise level. This is not a "
              "no-detection result.")
        return None

    res["spatial_null_test_status"] = "COMPLETED"
    res["spatial_correlation_length_km"] = float(corr_km)
    null = _null_cluster_sizes(v.shape, corr_km, res["geo"], trials=20)
    p95 = int(np.percentile(null, 95))
    print(f"\n  measured correlation length   : {corr_km:.1f} km")
    print(f"  clusters from noise ALONE     : median {int(np.median(null))} px, "
          f"95th pct {p95} px, max {int(max(null))} px")

    dof = len(res["days"]) - 2
    if dof <= 2:
        print(f"\n  NOTE: {len(res['days'])} epochs leaves dof={dof}. Each pixel's own")
        print("  sigma is itself uncertain by roughly half, so the per-fit")
        print("  significance count above is inflated and means little.")

    # A cluster faster than lambda/4 per pair is not a fast measurement, it is
    # an unmeasurable one. Above that the phase has wrapped and the unwrapper
    # has guessed; the value is as likely to be its mistake as the ground's.
    from gunw_reader import NISAR_LAMBDA_M
    ceiling = NISAR_LAMBDA_M / 4 * 1000 / 12.0
    if count and abs(float(np.nanmean(v[lab]))) > ceiling:
        print(f"\n  WARNING: cluster mean exceeds lambda/4 per 12-day pair "
              f"({ceiling:.2f} mm/day).")
        print("  Phase cannot represent motion this fast. These pixels are at least")
        print("  as likely to be an unwrapping error as a measurement, and cannot")
        print("  support a detection on their own.")

    if count <= p95:
        res["spatial_null_test_status"] = "COMPLETED_NO_DETECTION"
        print("\n  VERDICT: NO DETECTION. The largest cluster is within what this")
        print("  scene's own estimated correlation length produces from noise with no")
        print("  signal present. Cluster size is not evidence here.")
        return False

    res["spatial_null_test_status"] = "COMPLETED_CLUSTER_EXCEEDS_NOISE_NULL"
    print("\n  VERDICT: cluster exceeds the noise-only 95th percentile. That is")
    print("  necessary, not sufficient. Before calling it anything: check it")
    print("  appears in the other geometry, in a different time window, and")
    print("  moves consistently in one direction.")
    return True


MIN_SPATIAL_CORRELATION_PAIRS = 30


def _correlation_length_km(v: np.ndarray, ok: np.ndarray, geo: dict) -> float | None:
    """
    Estimate the 1/e decorrelation distance, or withhold it when neighbouring
    valid pixels do not support an empirical spatial-correlation estimate.

    Pixel spacing is not a measured correlation length. Substituting it for a
    sparse/disconnected mask creates an unjustifiably permissive noise null.
    """
    values_array = np.asarray(v, dtype="float64")
    valid = np.asarray(ok, dtype=bool) & np.isfinite(values_array)
    if int(valid.sum()) < MIN_SPATIAL_CORRELATION_PAIRS:
        return None

    values = values_array[valid]
    mean_value = float(np.mean(values))
    ref = float(np.mean((values - mean_value) ** 2))
    if not np.isfinite(ref) or ref <= 0:
        return None

    a = np.full(values_array.shape, np.nan, dtype="float64")
    a[valid] = values - mean_value
    px_km = abs(float(geo["transform"].a)) / 1000.0
    if not np.isfinite(px_km) or px_km <= 0:
        return None

    # Measure both axes and keep the larger decorrelation length. Require
    # enough actual neighbour pairs at every tested lag.
    def decay(axis: int) -> float | None:
        n = a.shape[axis]
        for lag in range(1, min(60, n // 2)):
            if axis == 1:
                left, right = a[:, :-lag], a[:, lag:]
            else:
                left, right = a[:-lag, :], a[lag:, :]
            paired = np.isfinite(left) & np.isfinite(right)
            if int(paired.sum()) < MIN_SPATIAL_CORRELATION_PAIRS:
                return None
            correlation = float(np.mean(left[paired] * right[paired]) / ref)
            if not np.isfinite(correlation):
                return None
            if correlation < np.e ** -1:
                return lag * px_km
        return 60 * px_km

    x_length = decay(1)
    y_length = decay(0)
    if x_length is None or y_length is None:
        return None
    return max(x_length, y_length)


def _null_cluster_sizes(shape, corr_km: float, geo: dict, trials: int = 20):
    """Largest 3-sigma cluster from pure noise at this scene's own scale."""
    from numpy.fft import irfft2, rfft2
    px_km = abs(geo["transform"].a) / 1000.0
    k = max(3, int(round(corr_km / max(px_km, 1e-6))))
    k = min(k, min(shape) // 2 or 3)
    ker = np.ones((k, k)) / (k * k)
    rng = np.random.default_rng(0)
    out = []
    for _ in range(trials):
        f = rng.normal(0.0, 1.0, shape)
        sm = irfft2(rfft2(f) * rfft2(ker, shape), shape)
        s = sm.std()
        if s <= 0:
            out.append(0); continue
        _, n = _largest_cluster(np.abs(sm / s) > 3.0)
        out.append(n)
    return out


def _largest_cluster(mask: np.ndarray):
    """Largest 4-connected blob, by flood fill. No scipy for one call."""
    seen = np.zeros_like(mask, bool)
    best = np.zeros_like(mask, bool); best_n = 0
    H, W = mask.shape
    for si in range(H):
        for sj in range(W):
            if not mask[si, sj] or seen[si, sj]:
                continue
            stack_, cells = [(si, sj)], []
            seen[si, sj] = True
            while stack_:
                i, j = stack_.pop(); cells.append((i, j))
                for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    a, b = i + di, j + dj
                    if 0 <= a < H and 0 <= b < W and mask[a, b] and not seen[a, b]:
                        seen[a, b] = True; stack_.append((a, b))
            if len(cells) > best_n:
                best_n = len(cells)
                best = np.zeros_like(mask, bool)
                for i, j in cells:
                    best[i, j] = True
    return best, best_n


def cluster_footprint(res: dict):
    """Largest significant cluster, in longitude and latitude."""
    v, s, ok = res["velocity"], res["sigma"], res["valid"]
    if not ok.any():
        return None
    vv, ss = v[ok], s[ok]
    rob = 1.4826 * np.median(np.abs(vv - np.median(vv)))
    own = np.zeros_like(v, bool); own[ok] = np.abs(vv) > 3 * ss
    pop = np.zeros_like(v, bool); pop[ok] = np.abs(vv - np.median(vv)) > 3 * rob
    lab, n = _largest_cluster(own & pop)
    if not n:
        return None
    ii, jj = np.nonzero(lab)
    T = res["geo"]["transform"]
    xs = T.c + (jj + 0.5) * T.a
    ys = T.f + (ii + 0.5) * T.e
    try:
        from pyproj import Transformer
        lon, lat = Transformer.from_crs(res["geo"]["crs"], 4326,
                                        always_xy=True).transform(xs, ys)
    except ImportError:
        logger.warning("pyproj missing - cannot compare geometries on the ground.")
        return None
    return {"n": n, "lon": (float(lon.min()), float(lon.max())),
            "lat": (float(lat.min()), float(lat.max())),
            "velocity": float(np.nanmean(v[lab]))}


def cross_geometry_check(found: dict) -> None:
    """
    Do the geometries see it in the same PLACE?

    This is the test that decides. Ascending and descending fly at different
    times through different weather, so a residual atmospheric blob has no
    reason to land twice on the same ground - while a moving slope has no
    choice. Comparing indices would be meaningless because each geometry has
    its own grid; the comparison has to happen in longitude and latitude.
    """
    if len(found) < 2:
        print("\n  Only one geometry available - the strongest check on a candidate")
        print("  cannot be run. Treat any cluster as unconfirmed.")
        return

    print(f"\n{'=' * 72}")
    print("CROSS-GEOMETRY CHECK")
    print("=" * 72)
    for g, c in sorted(found.items()):
        print(f"  {g}  {c['n']:>4} px  {c['velocity']:+7.2f} mm/day   "
              f"lon {c['lon'][0]:.4f}..{c['lon'][1]:.4f}  "
              f"lat {c['lat'][0]:.4f}..{c['lat'][1]:.4f}")

    geoms = sorted(found)
    for a, b in [(geoms[i], geoms[j])
                 for i in range(len(geoms)) for j in range(i + 1, len(geoms))]:
        ca, cb = found[a], found[b]
        ov_lon = min(ca["lon"][1], cb["lon"][1]) - max(ca["lon"][0], cb["lon"][0])
        ov_lat = min(ca["lat"][1], cb["lat"][1]) - max(ca["lat"][0], cb["lat"][0])
        print(f"\n  {a} vs {b}: overlap {ov_lon:+.4f} deg lon, "
              f"{ov_lat:+.4f} deg lat")
        if ov_lon > 0 and ov_lat > 0:
            print("  SAME GROUND. Two geometries, two weather systems, one location -")
            print("  that is a real candidate. Check the sign against the slope")
            print("  aspect before interpreting direction.")
        else:
            km = abs(min(ov_lon, ov_lat)) * 111.0
            print(f"  DISJOINT - about {km:.1f} km apart, sharing no ground.")
            print("  Independent geometries throwing up unrelated blobs is exactly")
            print("  what residual atmosphere does. NOT a detection.")


def write_geotiff(res: dict, out: Path) -> None:
    import rasterio
    v, s = res["velocity"], res["sigma"]
    # Bands 1-2 hold only complete-chain velocity and its uncertainty.
    # Bands 3-4 preserve partial support explicitly; they are not displacement
    # measurements and must not be interpreted as zero velocity.
    out.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out, "w", driver="GTiff", height=v.shape[0], width=v.shape[1],
                       count=4, dtype="float32", crs=res["geo"]["crs"],
                       transform=res["geo"]["transform"], nodata=np.nan,
                       compress="deflate", predictor=3) as dst:
        dst.write(v, 1); dst.set_band_description(1, "LOS velocity (mm/day; complete connected chain only)")
        dst.write(s, 2); dst.set_band_description(2, "Residual-based OLS sigma (mm/day; diagnostic only)")
        dst.write(res["observed_pair_count"].astype("float32"), 3)
        dst.set_band_description(3, "Observed valid pair count (support diagnostic; not motion)")
        dst.write(res["longest_contiguous_pair_count"].astype("float32"), 4)
        dst.set_band_description(4, "Longest contiguous valid pair count (support diagnostic; not motion)")
        pair_records = [
            {
                "source": str(Path(p.get("path", "")).name),
                "reference_date": p["ref"].isoformat(),
                "secondary_date": p["sec"].isoformat(),
            }
            for p in res["run"]
        ]
        dst.update_tags(
            support_policy="velocity bands require complete connected pair chain",
            partial_support_pixels=str(int(res["partial_support"].sum())),
            complete_series_pixels=str(int(res["complete_series_pixels"])),
            no_pair_support_pixels=str(int((res["observed_pair_count"] == 0).sum())),
            degrees_of_freedom=str(max(len(res["days"]) - 2, 0)),
            uncertainty_status=str(res["uncertainty_status"]),
            common_mode_removed=str(res["common_mode_removed"]).lower(),
            common_mode_policy=str(res["common_mode_policy"]),
            common_mode_pair_medians_mm_json=json.dumps(res["common_mode_pair_medians_mm"], separators=(",", ":")),
            pair_records_json=json.dumps(pair_records, separators=(",", ":")),
            geometry=str(res["run"][0].get("geom", "unknown")),
            first_acquisition_date=str(res["run"][0]["ref"].isoformat()),
            last_acquisition_date=str(res["run"][-1]["sec"].isoformat()),
            source_pair_count=str(len(res["run"])),
            spatial_null_test_status=str(res.get("spatial_null_test_status", "NOT_REPORTED")),
            spatial_correlation_length_km=(
                "" if res.get("spatial_correlation_length_km") is None
                else str(res.get("spatial_correlation_length_km"))
            ),
            network_id=str(res["run"][0].get("network_id", "unknown")),
            scientific_claim_limit="Relative LOS trend diagnostic; not independent geophysical validation or a hazard finding.",
        )
    logger.info("Wrote %s", out)


def main() -> int:
    ap = argparse.ArgumentParser(description="Per-pixel velocity from exported AOI rasters")
    ap.add_argument("--dir", default="outputs/export")
    ap.add_argument("--season", choices=("winter", "summer", "all"), default="winter")
    ap.add_argument("--keep-common-mode", action="store_true",
                    help="do NOT remove the per-epoch scene median")
    ap.add_argument("--geotiff", metavar="OUT.tif")
    args = ap.parse_args()

    d = resolve(args.dir)
    rasters = [r for r in (parse_export(p) for p in sorted(d.glob("*.tif"))) if r]
    rasters = [r for r in rasters if r["proc"] == "PR"]
    if args.season != "all":
        want = args.season == "winter"
        rasters = [r for r in rasters if (r["ref"].month in WINTER_MONTHS) == want]
    if not rasters:
        logger.error("No %s rasters in %s", args.season, d)
        return 1

    logger.info("%d %s pair(s)", len(rasters), args.season)
    if args.keep_common_mode:
        logger.warning("Common mode kept - atmosphere will dominate every pixel.")

    exact_network_count = len({r["network_id"] for r in rasters})
    all_runs = chain(rasters)
    runs = [run for run in all_runs if len(run) >= 2]
    print(f"PIXEL_STACK_PAIRS_EXAMINED={len(rasters)}")
    print(f"PIXEL_STACK_EXACT_NETWORKS={exact_network_count}")
    print(f"PIXEL_STACK_CONTIGUOUS_RUNS={len(all_runs)}")
    print(f"PIXEL_STACK_NETWORKS_WITH_2_PLUS_PAIRS={len(runs)}")
    if not runs:
        print("PIXEL_STACK_STATUS=WITHHELD_NO_COMPATIBLE_REPEATED_NETWORK")
        print(
            "VERDICT: NOT TESTABLE. The GUNW pairs do not contain two consecutive "
            "intervals within the same frame, subswath, bandwidth, polarization "
            "and processing CRID. Combining different networks would mix LOS "
            "geometry or processing lineage."
        )
        return 0

    print("PIXEL_STACK_STATUS=COMPLETED_WITH_COMPATIBLE_NETWORKS")
    found = {}
    for k, run in enumerate(runs, 1):
        res = stack(run, remove_common_mode=not args.keep_common_mode)
        survived = report(res, f"{run[0]['geom']} {run[0]['network_id']} chain {k}")
        fp = cluster_footprint(res)
        # A withheld null test is not a no-detection result and cannot
        # contribute to cross-geometry corroboration.
        if fp and survived is True:
            found.setdefault(run[0]["geom"], fp)
        elif fp and survived is False:
            logger.info("%s cluster did not beat its completed noise null - excluded "
                        "from the cross-geometry check.", run[0]["network_id"])
        elif fp and survived is None:
            logger.warning("%s spatial noise null was not estimable; candidate "
                           "withheld from the cross-geometry check.", run[0]["network_id"])
        if args.geotiff:
            base = Path(args.geotiff)
            if len(runs) == 1:
                output_path = base
            else:
                network = re.sub(r"[^A-Za-z0-9_-]+", "_", str(run[0]["network_id"]))
                output_path = base.with_name(f"{base.stem}_{network}_chain{k:02d}{base.suffix or '.tif'}")
            write_geotiff(res, output_path)
            logger.info("Wrote target-free per-pixel stack diagnostic: %s", output_path)

    cross_geometry_check(found)
    return 0


if __name__ == "__main__":
    sys.exit(main())
