"""
common_ref.py
-------------
Put a stack of exported pairs on one datum, and measure what the lack of one
cost.

Why this exists
===============
SBAS needs every interferogram referenced to the SAME ground. The pipeline
never did that for GOFF: goff_reader picks its own reference block and fits
its own plane for every pair, so a differential atmosphere or ramp residual
between two different reference areas enters the time series as a step at
that epoch. find_common_reference() in gunw_reader existed, its docstring said
so, and timeseries.py never called it. That was listed as the largest unfixed
defect under the headline bound.

The products are gone from disk, but the exports are not, and every export
sits on the same fixed AOI lattice (gunw_reader.aoi_grid). That is enough to
fix it after the fact:

    1. keep the pixels valid in EVERY pair of the stack,
    2. drop those within --buffer-km of any target, so the datum is never
       set on the ground being measured,
    3. subtract, from each pair, its median over what is left.

Each pair's constant is then exactly the per-pair datum error: what the target
series carried as a spurious step. Dividing it by the span says, in the same
units as the floor, how much of any quoted velocity it could have been.

This re-references; it does not re-deramp. A residual tilt across the AOI that
differs between pairs is not removed by a constant. The stable set is spread
across the AOI, so the median is robust to it, but it is a limit.

Usage
-----
    python src/common_ref.py --dir outputs/export_goff_src --match layer2 \\
        --include 20260702 20260714 20260726 --exclude _UR_ \\
        --targets outputs/scar_candidates_source.csv --reading DETACHMENT-like \\
        --csv outputs/common_ref_summer.csv
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("common_ref")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from local_floor import (DATE, span_from_name, window, point_in_aoi,
                         aoi_tokens_from_dir)  # noqa: E402
from gunw_reader import AOIS  # noqa: E402
from paths import resolve  # noqa: E402

MIN_COMMON = 200          # fewer stable pixels than this and a median is noise


# ---------------------------------------------------------------------------
def common_mask(stack: list[np.ndarray]) -> np.ndarray:
    """Cells finite in every array of the stack. All arrays share one lattice."""
    m = np.ones(stack[0].shape, dtype=bool)
    for a in stack:
        if a.shape != m.shape:
            raise ValueError("stack is not on one lattice - re-run --export")
        m &= np.isfinite(a)
    return m


def buffer_mask(shape: tuple, centres: list[tuple[int, int]], radius_px: int) -> np.ndarray:
    """True within radius_px (Euclidean) of any centre."""
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    out = np.zeros(shape, dtype=bool)
    for r, c in centres:
        out |= (yy - r) ** 2 + (xx - c) ** 2 <= radius_px ** 2
    return out


def rereference(stack: list[np.ndarray], stable: np.ndarray) -> tuple[list[np.ndarray], list[float]]:
    """
    Subtract each pair's median over the common stable set.

    Returns the re-referenced arrays and the constant removed from each. The
    constants are the result as much as the arrays are: they are the datum
    error each pair carried under its own private reference.
    """
    if int(stable.sum()) < MIN_COMMON:
        raise ValueError(f"only {int(stable.sum())} common stable pixels; "
                         f"need {MIN_COMMON}")
    consts = [float(np.median(a[stable])) for a in stack]
    return [a - k for a, k in zip(stack, consts)], consts


def window_median(a: np.ndarray, rc: tuple[int, int], radius: int) -> tuple[float, int]:
    w = window(a, rc[0], rc[1], radius)
    v = w[np.isfinite(w)]
    return (float(np.median(v)) if v.size else float("nan"), int(v.size))


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description="Re-reference a stack of exports to one datum")
    ap.add_argument("--dir", required=True)
    ap.add_argument("--aoi", choices=sorted(AOIS), default=None,
                    help="exact AOI polygon represented by the exports; inferred only from an unambiguous directory name")
    ap.add_argument("--allow-event-affected-aoi", action="store_true",
                    help="explicitly allow LHENDE, which is an event-affected runout/impoundment corridor, for a non-control analysis")
    ap.add_argument("--match", default=None)
    ap.add_argument("--include", nargs="*", default=[],
                    help="reference dates to keep; one geometry and one "
                         "connected block per run, since blocks cannot share a zero")
    ap.add_argument("--exclude", nargs="*", default=[])
    ap.add_argument("--band", type=int, default=1)
    ap.add_argument("--targets", required=True, help="CSV with id, lat, lon")
    ap.add_argument("--reading", default=None)
    ap.add_argument("--radius", type=int, default=6, help="target window half-width, px")
    ap.add_argument("--buffer-km", type=float, default=2.0,
                    help="no datum pixel within this distance of any target")
    ap.add_argument("--csv", metavar="OUT.csv")
    args = ap.parse_args()

    import rasterio
    from pyproj import Transformer
    from rasterio.features import geometry_mask
    from rasterio.warp import transform_geom

    folder_aois = aoi_tokens_from_dir(args.dir)
    if len(folder_aois) > 1:
        ap.error(f"--dir name contains conflicting AOI markers {sorted(folder_aois)}; use a dedicated export directory.")
    inferred_aoi = next(iter(folder_aois)) if folder_aois else None
    if args.aoi is None and inferred_aoi is None:
        ap.error("could not infer AOI from --dir; pass --aoi source|langtang|lhende explicitly")
    if args.aoi is not None and inferred_aoi is not None and args.aoi != inferred_aoi:
        ap.error(f"--aoi {args.aoi!r} conflicts with the AOI marker in --dir ({inferred_aoi!r}); use matching exports.")
    aoi_name = args.aoi or inferred_aoi
    if aoi_name == "lhende" and not args.allow_event_affected_aoi:
        ap.error("LHENDE is event-affected and is not a stable-ground control; pass --allow-event-affected-aoi only for an explicitly justified analysis.")

    files = sorted(Path(resolve(args.dir)).glob("*.tif"))
    if args.match:
        files = [f for f in files if args.match in f.name]
    files = [f for f in files if not any(x in f.name for x in args.exclude)]
    if args.include:
        files = [f for f in files
                 if DATE.search(f.name) and DATE.search(f.name).group(1) in args.include]
    if len(files) < 2:
        logger.error("need at least two pairs, found %d", len(files))
        return 1

    if args.band < 1:
        ap.error("--band must be a positive 1-based raster band")
    if args.radius < 0 or args.buffer_km < 0:
        ap.error("--radius and --buffer-km must be non-negative")

    stack, names, spans, px = [], [], [], None
    grid_crs = grid_transform = grid_shape = None
    pair_names = []
    for f in files:
        pair_match = DATE.search(f.name)
        span = span_from_name(f.name)
        if pair_match is None or not np.isfinite(span) or span <= 0:
            logger.error("Invalid or non-positive acquisition span in %s", f.name)
            return 1
        pair_names.append(pair_match.group(0))
        with rasterio.open(f) as s:
            if s.crs is None or not s.crs.is_projected:
                logger.error("%s is not in a projected CRS; metric target buffers are unsafe.", f.name)
                return 1
            if s.crs.to_epsg() is None:
                logger.error("%s has no identifiable EPSG CRS; refusing grid alignment by assumption.", f.name)
                return 1
            if grid_crs is None:
                grid_crs, grid_transform, grid_shape = s.crs, s.transform, (s.height, s.width)
                px = abs(float(s.res[0]))
                if px <= 0 or not np.isfinite(px):
                    logger.error("Invalid raster posting in %s", f.name)
                    return 1
                tf = Transformer.from_crs(4326, s.crs, always_xy=True)
                idx = s.index
            else:
                same_transform = np.allclose(tuple(s.transform), tuple(grid_transform), rtol=0, atol=1e-9)
                if s.crs != grid_crs or (s.height, s.width) != grid_shape or not same_transform:
                    logger.error("Grid mismatch in %s; all common-reference inputs must share CRS, shape and affine transform.", f.name)
                    return 1
            a = s.read(args.band).astype(float)
            if s.nodata is not None:
                a[a == s.nodata] = np.nan
        stack.append(a)
        names.append(pair_match.group(0))
        spans.append(span)

    if len(pair_names) != len(set(pair_names)):
        logger.error("Multiple exports for the same acquisition pair were selected; constrain --match to one polarization/layer, e.g. HH-layer2.")
        return 1

    ring = [list(p) for p in AOIS[aoi_name]]
    ring.append(ring[0])
    polygon = transform_geom("EPSG:4326", grid_crs, {"type": "Polygon", "coordinates": [ring]})
    inside_aoi = geometry_mask([polygon], out_shape=grid_shape,
                               transform=grid_transform, all_touched=False,
                               invert=True)

    with open(resolve(args.targets), newline="") as fh:
        target_rows = [(r["id"], float(r["lat"]), float(r["lon"])) for r in csv.DictReader(fh)
                       if args.reading is None or r.get("reading") == args.reading]
    targets = []
    for tid, lat, lon in target_rows:
        if not point_in_aoi(lon, lat, aoi_name):
            logger.warning("Ignoring target %s outside AOI polygon %s (%.6f, %.6f).", tid, aoi_name, lat, lon)
            continue
        targets.append((tid, lat, lon))
    if not targets:
        logger.error("No selected candidate lies inside AOI polygon %s.", aoi_name)
        return 1

    centres = {}
    for tid, lat, lon in targets:
        row, col = idx(*tf.transform(lon, lat))
        if not (0 <= row < grid_shape[0] and 0 <= col < grid_shape[1]):
            logger.error("Target %s is inside the declared AOI but outside the export grid; re-export the same AOI.", tid)
            return 1
        centres[tid] = (row, col)

    common = common_mask(stack) & inside_aoi
    stable = common & ~buffer_mask(common.shape, list(centres.values()),
                                   int(round(args.buffer_km * 1000 / px)))
    fixed, consts = rereference(stack, stable)

    print(f"\n{'='*78}\nCOMMON DATUM over {len(stack)} pairs: {int(common.sum()):,} cells valid "
          f"in all, {int(stable.sum()):,} after the {args.buffer_km:g} km target buffer\n{'='*78}")
    print(f"  {'PAIR':<20}{'span':>5}{'datum error':>13}{'as velocity':>13}")
    for n, sp, k in zip(names, spans, consts):
        print(f"  {n:<20}{sp:>5.0f}{k:>10.1f} mm{k/sp:>9.2f} mm/d")
    print("  Each constant is what that pair's private reference left in the")
    print("  series. Under per-pair referencing it arrives as a step at that epoch.")

    rows = []
    for tid, (r0, c0) in centres.items():
        print(f"\n  target {tid}")
        print(f"  {'PAIR':<20}{'per-pair ref':>14}{'common ref':>12}{'px':>5}")
        for n, sp, raw, fx, k in zip(names, spans, stack, fixed, consts):
            m_raw, npx = window_median(raw, (r0, c0), args.radius)
            m_fix, _ = window_median(fx, (r0, c0), args.radius)
            print(f"  {n:<20}{m_raw:>11.1f} mm{m_fix:>9.1f} mm{npx:>5}")
            rows.append({"target_id": tid, "pair": n, "span_days": sp,
                         "datum_error_mm": round(k, 3),
                         "datum_error_mm_day": round(k / sp, 3),
                         "target_mm_per_pair_ref": round(m_raw, 3),
                         "target_mm_common_ref": round(m_fix, 3),
                         "target_mm_day_common_ref": round(m_fix / sp, 3),
                         "window_px": npx, "stable_px": int(stable.sum())})

    worst = max(abs(k / sp) for k, sp in zip(consts, spans))
    print(f"\n  Largest datum error in this stack: {worst:.2f} mm/day.")

    if args.csv:
        p = Path(resolve(args.csv))
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        logger.info("Wrote %s", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
