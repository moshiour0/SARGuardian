"""
local_floor.py
--------------
The detection floor at a point, not over the AOI.

Why this exists
===============
Every noise floor this project quotes is the MAD scatter of one product over
the whole area of interest. The source polygon is 82 km2. The failure was
localised. Those are not the same measurement, and until this module existed
nobody had checked how different they are.

They are different by a factor of 1.7, and the difference is worst exactly
where it matters. Ascending path 098, GOFF layer2, routine products, in a
13x13 pixel window (about 1 km) centred on the failure point at
28.28771 N 85.52809 E:

    pair                    span   AOI floor   point floor   valid px
    20251128_20251210        12      21.0         28.0       102/169
    20251210_20251222        12      21.4         42.8        93/169
    20251222_20260103        12      24.5         47.5       111/169
    20260103_20260115        12      23.8         32.5       103/169
    20260702_20260714        12      18.6         40.4        49/169
    20260714_20260726        12      13.4         19.2        65/169
    20260726_20260819        24       8.9         34.3        26/169   <-
    20260819_20260831        12      15.6         16.5        19/169
    ------------------------------------------------------------------
    median                           19.8         33.4

The marked row is the last ascending interval before the 26 August failure.
Over the AOI it has the LOWEST floor in the entire archive - 8.9 mm/day, the
number that made the published bound look strongest. At the failure point that
same pair has a floor of 34.3 mm/day and only 26 of 169 valid pixels.

The bound survives: nothing measured at the point exceeds even the local floor,
in windows of 250 m, 500 m and 1 km.

But do not quote 33.4 as the bound. It is the median across all eight
ascending pairs, five of them winter pairs from November to January, outside
the seven-week window being bounded. Over that window the three covering
intervals read 40.4, 19.2 and 34.3 at the point, and a bound that holds across
a window is set by its weakest interval. That point has since been shown not
to be a scar; the same measurement at every candidate (--targets) is what
bound.py quotes, and at the candidate at the published elevation the figure is
32.0 mm/day per pixel and 9.2 for a 1 km window median.

The AOI figure of 18.6 understates all of these.

The general lesson, which is the same one this project keeps relearning: a
statistic aggregated over terrain that did not fail describes terrain that did
not fail. An 82 km2 median is a fine description of an 82 km2 area and a poor
description of one hillside inside it.

What "floor" means here
=======================
3 * MAD-sigma of the range offsets in the window, divided by the span in days:
the constant velocity that would have to be present for the pair to see it at
three sigma. Identical definition to goff_reader's AOI-wide
`detect_floor_mm_day`, so the two are directly comparable - the only thing that
changes is which pixels go in.

Usage
-----
    python src/local_floor.py --dir outputs/export_goff_src --aoi source --match layer2 \\
        --lat 28.28771 --lon 85.52809 --radius 6

    python src/local_floor.py --dir outputs/export_goff_src --aoi source --match layer2 \\
        --lat 28.28771 --lon 85.52809 --sweep 3 6 12 --csv outputs/local_floor.csv
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import logging
import re
import sys
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("local_floor")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import resolve  # noqa: E402
from gunw_reader import AOIS  # noqa: E402

MIN_PX = 6                # below this a MAD is not a measurement
DATE = re.compile(r"(\d{8})_(\d{8})")


# ---------------------------------------------------------------------------
def robust_sigma(a: np.ndarray) -> float:
    """MAD-based scatter, the same definition goff_reader uses."""
    v = a[np.isfinite(a)]
    if v.size == 0:
        return float("nan")
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


def detection_floor(values_mm: np.ndarray, span_days: float) -> float:
    """
    The constant velocity this pair could see at three sigma, in mm/day.

    Same formula as goff_reader.report: 3 * sigma / span. Kept identical on
    purpose - the point of this module is that the AOI and the point give
    different answers to the SAME question, so the question must not change.
    """
    if not np.isfinite(span_days) or span_days <= 0:
        return float("nan")
    sigma = robust_sigma(values_mm)
    # A zero MAD is not evidence of zero instrument noise. It can result from
    # quantization, repeated fill/rounding values, or too little effective
    # information; never publish it as a perfect detection floor.
    if not np.isfinite(sigma) or sigma <= 0:
        return float("nan")
    return 3.0 * sigma / span_days


def window(arr: np.ndarray, row: int, col: int, radius: int) -> np.ndarray:
    """
    The square window of half-width `radius` about (row, col), clipped to the
    array. Clipping rather than padding, so a target near the edge yields
    fewer pixels instead of a window full of invented nodata.
    """
    r0, r1 = max(0, row - radius), min(arr.shape[0], row + radius + 1)
    c0, c1 = max(0, col - radius), min(arr.shape[1], col + radius + 1)
    if r0 >= r1 or c0 >= c1:
        return np.empty((0, 0), dtype=arr.dtype)
    return arr[r0:r1, c0:c1]


def floor_ci(values_mm: np.ndarray, span_days: float, n_boot: int = 2000,
             seed: int = 20260909) -> tuple:
    """
    Percentile-bootstrap 95% interval for the detection floor.

    Why a floor needs an interval
    ----------------------------
    The floor is 3*MAD-sigma/span, and a MAD from a few dozen pixels is an
    estimate with real sampling error - roughly 15% relative at n = 50, and
    far worse below that. The pre-event pairs that matter most here are the
    ones with the FEWEST valid pixels in the window (26 and 19 of 169), so the
    headline bound rests on exactly the samples where a point estimate is
    least trustworthy. Quoting "40.4 mm/day" to three figures from 49 pixels
    asserts a precision the data does not carry.

    One caveat this interval does NOT cover: GOFF layers come from overlapping
    correlation windows, so neighbouring offset estimates share input pixels
    and are not independent. The effective sample size is smaller than the
    pixel count, so this interval is a LOWER bound on the true uncertainty.
    The error neighbours share is measured separately, by
    block_median_floor(), rather than folded in here.
    """
    v = values_mm[np.isfinite(values_mm)]
    n = v.size
    if n < MIN_PX or not np.isfinite(span_days) or span_days <= 0:
        return (float("nan"), float("nan"))
    base_sigma = robust_sigma(v)
    if not np.isfinite(base_sigma) or base_sigma <= 0:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    samp = v[idx]
    med = np.median(samp, axis=1, keepdims=True)
    mad = np.median(np.abs(samp - med), axis=1)
    floors = 3.0 * 1.4826 * mad / span_days
    good = floors[np.isfinite(floors)]
    if good.size == 0:
        return (float("nan"), float("nan"))
    return (float(np.percentile(good, 2.5)), float(np.percentile(good, 97.5)))


def block_median_floor(field_mm: np.ndarray, size: int, span_days: float,
                       min_valid: int = 30) -> tuple[float, int]:
    """
    The 3-sigma floor on a WINDOW MEDIAN, measured empirically. mm/day.

    detection_floor() answers "what velocity would ONE PIXEL see at three
    sigma". A precursor measured by another instrument is an average over a
    slope, so the fair comparison is against the noise of an averaged
    estimate - and that cannot be had by dividing the per-pixel sigma by
    sqrt(n). GOFF offsets come from overlapping correlation windows, so
    neighbours within a kilometre share most of their error; a scatter
    computed inside one window cannot see the part the whole window shares.

    So measure it directly: tile the field into non-overlapping blocks the
    size of the window, take each block's median, and quote 3 MADs of those
    medians over the span. That includes every error the window median
    actually suffers, correlated or not, and any real motion elsewhere in the
    AOI - so it is an upper bound on the estimator noise, not an idealisation.
    It is AOI-wide by construction; the point-versus-area escalation this
    module exists for applies to it too, and is why it is quoted beside the
    per-pixel local floor rather than instead of it.

    Returns (floor, number of blocks). NaN when fewer than 6 blocks qualify.
    """
    if field_mm.ndim != 2 or size < 1 or span_days <= 0:
        return (float("nan"), 0)
    meds = []
    H, W = field_mm.shape
    for i in range(0, H - size + 1, size):
        for j in range(0, W - size + 1, size):
            v = field_mm[i:i + size, j:j + size]
            v = v[np.isfinite(v)]
            if v.size >= min(min_valid, size * size):
                meds.append(float(np.median(v)))
    if len(meds) < MIN_PX:
        return (float("nan"), len(meds))
    sigma = robust_sigma(np.array(meds))
    if not np.isfinite(sigma) or sigma <= 0:
        return (float("nan"), len(meds))
    return (3.0 * sigma / span_days, len(meds))


def compare(aoi_mm: np.ndarray, win_mm: np.ndarray, span_days: float) -> dict:
    """
    Both floors and the ratio between them, from one pair.

    Pure: takes arrays, returns numbers, no raster and no I/O, so the claim
    this module makes can be tested without a product. When the AOI field is
    2-D, the floor on a window MEDIAN of the same size is measured as well;
    see block_median_floor().
    """
    a = aoi_mm[np.isfinite(aoi_mm)]
    w = win_mm[np.isfinite(win_mm)]
    out = {"span_days": float(span_days),
           "aoi_px": int(a.size), "window_px": int(w.size),
           "window_total_px": int(win_mm.size),
           "aoi_floor_mm_day": float("nan"),
           "local_floor_mm_day": float("nan"),
           "local_floor_ci_lo": float("nan"),
           "local_floor_ci_hi": float("nan"),
           "block_median_floor_mm_day": float("nan"), "n_blocks": 0,
           "ratio": float("nan"), "usable": False,
           "status": "INSUFFICIENT_LOCAL_PIXELS"}
    if not np.isfinite(span_days) or span_days <= 0:
        out["status"] = "INVALID_DATE_SPAN"
        return out
    if a.size < MIN_PX or w.size < MIN_PX:
        return out

    fa = detection_floor(a, span_days)
    fl = detection_floor(w, span_days)
    if not np.isfinite(fl) or fl <= 0:
        out.update(aoi_floor_mm_day=fa, status="DEGENERATE_LOCAL_MAD")
        return out
    lo, hi = floor_ci(w, span_days)
    size = max(win_mm.shape) if win_mm.ndim == 2 else 0
    bm, nb = block_median_floor(aoi_mm, size, span_days)
    ratio = fl / fa if np.isfinite(fa) and fa > 0 else float("nan")
    status = "MEASURED" if np.isfinite(fa) and fa > 0 else "MEASURED_LOCAL_FLOOR_AOI_MAD_UNUSABLE"
    out.update(aoi_floor_mm_day=fa, local_floor_mm_day=fl,
               local_floor_ci_lo=lo, local_floor_ci_hi=hi,
               block_median_floor_mm_day=bm, n_blocks=nb,
               ratio=ratio, usable=True, status=status)
    return out


def span_from_name(name: str) -> float:
    """Days between the two YYYYMMDD stamps in an export filename."""
    m = DATE.search(name)
    if not m:
        return float("nan")
    a, b = (dt.date(int(s[:4]), int(s[4:6]), int(s[6:8])) for s in m.groups())
    return float((b - a).days)



def point_in_ring(lon: float, lat: float, ring) -> bool:
    """Return whether a WGS84 lon/lat point is inside or on an AOI ring.

    These AOIs cover a small Himalayan region, so a ray-casting test in
    geographic coordinates is sufficient for candidate eligibility. Raster
    georeferencing and metric distances still use the raster CRS / pyproj.
    A point exactly on a segment is treated as inside.
    """
    inside = False
    n = len(ring)
    if n < 3:
        return False
    eps = 1e-10
    for i in range(n):
        x1, y1 = ring[i]
        x2, y2 = ring[(i + 1) % n]
        cross = (lon - x1) * (y2 - y1) - (lat - y1) * (x2 - x1)
        if abs(cross) <= eps and min(x1, x2) - eps <= lon <= max(x1, x2) + eps and min(y1, y2) - eps <= lat <= max(y1, y2) + eps:
            return True
        if (y1 > lat) != (y2 > lat):
            x_intersect = x1 + (lat - y1) * (x2 - x1) / (y2 - y1)
            if lon < x_intersect:
                inside = not inside
    return inside


def point_in_aoi(lon: float, lat: float, aoi_name: str) -> bool:
    """Check target eligibility against the named polygon, not its bbox."""
    if aoi_name not in AOIS:
        raise ValueError(f"Unknown AOI {aoi_name!r}; choose from {sorted(AOIS)}")
    return point_in_ring(lon, lat, AOIS[aoi_name])


def aoi_tokens_from_dir(directory: str | Path) -> set[str]:
    """Return exact AOI markers in a folder basename; never infer by substring."""
    tokens = set(filter(None, re.split(r"[^a-z0-9]+", Path(directory).name.lower())))
    found = set()
    if "source" in tokens or "src" in tokens:
        found.add("source")
    if "langtang" in tokens:
        found.add("langtang")
    if "lhende" in tokens:
        found.add("lhende")
    return found


def infer_aoi_from_dir(directory: str | Path) -> str | None:
    """Infer AOI only when the directory names exactly one polygon."""
    found = aoi_tokens_from_dir(directory)
    return next(iter(found)) if len(found) == 1 else None


def validate_export_provenance(files, aoi_name: str, *,
                               allow_legacy_untagged: bool = False) -> None:
    """Refuse mixed/mismatched AOIs and unreferenced CRS before measuring floors.

    Legacy GeoTIFFs created before AOI provenance tags existed may be explicitly
    inspected with --allow-legacy-untagged. A tagged file is never allowed to
    claim a different AOI or say that the AOI mask was not applied.
    """
    import rasterio

    for path in files:
        with rasterio.open(path) as src:
            epsg = src.crs.to_epsg() if src.crs else None
            if epsg is None:
                raise ValueError(f"{path.name}: CRS/EPSG metadata is missing or not resolvable.")
            tags = src.tags()
            file_aoi = tags.get("aoi_name")
            mask_applied = tags.get("aoi_mask_applied")
            if not file_aoi:
                if allow_legacy_untagged:
                    logger.warning(
                        "%s has no AOI provenance tags; using caller-supplied AOI=%s as an explicit legacy assumption",
                        path.name, aoi_name
                    )
                    continue
                raise ValueError(
                    f"{path.name}: missing aoi_name provenance tag; re-export with the current reader "
                    "or pass --allow-legacy-untagged only after confirming the file's AOI."
                )
            if file_aoi != aoi_name:
                raise ValueError(
                    f"{path.name}: AOI tag is {file_aoi!r}, but --aoi is {aoi_name!r}; "
                    "refusing a floor from a mismatched AOI."
                )
            if str(mask_applied).lower() != "true":
                raise ValueError(
                    f"{path.name}: aoi_mask_applied is {mask_applied!r}; refusing an unclipped/mislabelled raster."
                )


def flag_is_true(value, *, default: bool = False) -> bool:
    """Parse boolean flags from in-memory results and CSV round-trips."""
    if value is None:
        return default
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def skipped_row(file: Path, target_id: str, lat: float, lon: float,
                radius: int, aoi_name: str, status: str) -> dict:
    """Machine-readable rejection row. A skipped point must never create a floor."""
    span = span_from_name(file.name)
    c = compare(np.array([], dtype=float), np.array([], dtype=float),
                span if np.isfinite(span) else 0.0)
    return {
        "file": file.name, "radius_px": radius, "target_id": target_id,
        "lat": lat, "lon": lon, "aoi_name": aoi_name,
        "aoi_eligible": status != "TARGET_OUTSIDE_AOI",
        **c, "status": status,
    }


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description="Detection floor at a point against the floor over the AOI")
    ap.add_argument("--dir", required=True, help="directory of exported GeoTIFFs")
    ap.add_argument("--match", metavar="SUBSTR",
                    help="only files whose name contains this, e.g. layer2")
    ap.add_argument("--lat", type=float)
    ap.add_argument("--lon", type=float)
    ap.add_argument("--aoi", choices=sorted(AOIS), default=None,
                    help="AOI polygon used for candidate eligibility. Inferred "
                         "only from an unambiguous export-directory suffix; "
                         "otherwise required to prevent bbox-only eligibility.")
    ap.add_argument("--targets", metavar="CANDIDATES.csv",
                    help="measure at every row of a CSV with id, lat, lon "
                         "columns (scar_map.py --candidates-csv writes one), "
                         "instead of one --lat/--lon. A floor measured at one "
                         "candidate is not the floor at the others")
    ap.add_argument("--reading", default=None,
                    help="with --targets, keep only rows whose 'reading' "
                         "column equals this, e.g. DETACHMENT-like")
    ap.add_argument("--allow-event-affected-aoi", action="store_true",
                    help="explicitly allow noise-floor calculation inside the LHENDE runout/impoundment corridor; this is not a stable-ground control.")
    ap.add_argument("--allow-legacy-untagged", action="store_true",
                    help="allow old GeoTIFFs that lack AOI provenance tags, after independently confirming their AOI; tagged mismatches are always refused.")
    ap.add_argument("--band", type=int, default=1)
    ap.add_argument("--radius", type=int, default=6,
                    help="window half-width in pixels; 6 is about 1 km at 80 m")
    ap.add_argument("--sweep", type=int, nargs="+", metavar="R",
                    help="report several radii, to show the answer is not an "
                         "artefact of one window size")
    ap.add_argument("--exclude", nargs="*", default=[],
                    help="substrings to skip, e.g. a co-event pair")
    ap.add_argument("--include", nargs="*", default=[],
                    help="keep only files matching one of these substrings. "
                         "Export filenames carry no track field, so this is "
                         "how you restrict a run to one geometry - and you "
                         "should, because pooling tracks whose floors differ "
                         "fivefold gives a median that describes neither")
    ap.add_argument("--csv", metavar="OUT.csv")
    ap.add_argument("--as-floors", metavar="OUT.csv",
                    help="also write the local floors in the layout "
                         "inverse_velocity.py --floors reads, so a series "
                         "measured at a point is gated on the floor AT that "
                         "point rather than the more lenient AOI floor")
    args = ap.parse_args()

    try:
        import rasterio
        from pyproj import Transformer
    except ImportError:
        logger.error("rasterio and pyproj required")
        return 1

    files = sorted(Path(resolve(args.dir)).glob("*.tif"))
    if args.match:
        files = [f for f in files if args.match in f.name]
    files = [f for f in files if not any(x in f.name for x in args.exclude)]
    if args.include:
        files = [f for f in files if any(x in f.name for x in args.include)]
    if not files:
        logger.error("No matching .tif under %s", args.dir)
        return 1

    if args.targets:
        with open(resolve(args.targets), newline="") as fh:
            targets = [(r["id"], float(r["lat"]), float(r["lon"]))
                       for r in csv.DictReader(fh)
                       if args.reading is None or r.get("reading") == args.reading]
        if not targets:
            logger.error("No targets in %s", args.targets)
            return 1
    elif args.lat is not None and args.lon is not None:
        targets = [("", args.lat, args.lon)]
    else:
        ap.error("give --lat and --lon, or --targets")

    target_ids = [str(t[0]) for t in targets]
    if len(target_ids) != len(set(target_ids)):
        ap.error("candidate IDs in --targets must be unique; duplicate IDs could mix floors across locations")
    for target_id, lat, lon in targets:
        if not (np.isfinite(lat) and np.isfinite(lon)
                and -90 <= lat <= 90 and -180 <= lon <= 180):
            ap.error(f"target {target_id or '(single)'} has invalid WGS84 coordinates ({lat}, {lon})")

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
        ap.error("LHENDE is an event-affected runout/impoundment corridor, not a stable-ground control; pass --allow-event-affected-aoi only for an explicitly justified analysis.")
    if aoi_name == "lhende":
        logger.warning("LHENDE IS AN EVENT-AFFECTED RUNOUT/IMPOUNDMENT CORRIDOR; do not interpret its scatter as a stable-ground control.")

    try:
        validate_export_provenance(
            files, aoi_name, allow_legacy_untagged=args.allow_legacy_untagged
        )
    except (ValueError, RuntimeError) as exc:
        logger.error("Export provenance validation failed: %s", exc)
        return 1

    radii = args.sweep or [args.radius]
    if any(radius < 1 for radius in radii):
        ap.error("--radius and every --sweep radius must be at least 1 pixel")
    rows = []
    for tid, lat, lon in targets:
        if not point_in_aoi(lon, lat, aoi_name):
            logger.error("TARGET_OUTSIDE_AOI id=%s lat=%.6f lon=%.6f AOI=%s; refusing to compute a local floor for this target",
                         tid or "(single target)", lat, lon, aoi_name)
            rows.extend(skipped_row(f, tid, lat, lon, radius, aoi_name,
                                    "TARGET_OUTSIDE_AOI")
                        for f in files for radius in radii)
            continue
        target_rows = measure_target(files, tid, lat, lon, radii, args.band,
                                     aoi_name=aoi_name)
        rows.extend(target_rows)

    measured_rows = [
        x for x in rows
        if flag_is_true(x.get("aoi_eligible", True), default=True)
        and flag_is_true(x.get("usable"))
    ]
    if args.csv:
        p = Path(resolve(args.csv))
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=sorted({k for x in rows for k in x}))
            w.writeheader()
            w.writerows(rows)
        logger.info("Wrote %s", p)

    if not measured_rows:
        print("  No eligible target has enough valid pixels for a measured floor.")
        return 1

    r = np.array([x["ratio"] for x in measured_rows if np.isfinite(x["ratio"])])
    if r.size == 0:
        print("  Floors exist, but no finite AOI/local ratio could be formed.")
    elif np.nanmedian(r) > 1.25:
        print("  The floor at the point is materially WORSE than the floor over")
        print("  the AOI. Any bound quoted at this location must use the local")
        print("  figure; the AOI figure overstates what the product could see.")
    elif r.size and np.nanmedian(r) < 0.8:
        print("  The point is quieter than the AOI. The AOI floor is")
        print("  conservative here, which is safe but not tight.")
    else:
        print("  The two agree. The AOI floor is representative at this point.")

    if args.as_floors:
        write_as_floors(measured_rows, Path(resolve(args.as_floors)), args.match or "")
    return 0


def write_as_floors(rows: list[dict], path: Path, layer: str) -> int:
    """
    Usable local floors as reference, secondary, layer, detect_floor_mm_day.

    One target only: a floors file keyed by pair cannot say which point it
    describes, so mixing targets would gate one slope on another's noise.
    """
    usable = [r for r in rows if r["usable"]]
    if len({str(r.get("target_id", "")) for r in usable}) > 1:
        raise ValueError("--as-floors needs a single target")
    coordinates = {
        (round(float(r["lat"]), 8), round(float(r["lon"]), 8))
        for r in usable if r.get("lat") is not None and r.get("lon") is not None
    }
    if len(coordinates) > 1:
        raise ValueError("--as-floors rows contain one target ID at multiple coordinates")

    validated = []
    for row in usable:
        try:
            px = int(float(row.get("window_px", "nan")))
            floor = float(row.get("local_floor_mm_day", "nan"))
        except (TypeError, ValueError, OverflowError):
            raise ValueError("Refusing --as-floors row with invalid pixel count or floor.")
        if px < MIN_PX or not np.isfinite(floor) or floor <= 0:
            raise ValueError(
                f"Refusing --as-floors row for {row.get('file', '(unknown file)')}: "
                "a floor requires at least MIN_PX finite pixels and a finite positive value."
            )
        status = row.get("status")
        if status not in (None, "MEASURED", "MEASURED_LOCAL_FLOOR_AOI_MAD_UNUSABLE"):
            raise ValueError(f"Refusing --as-floors row with non-measurement status {status!r}.")
        if not flag_is_true(row.get("aoi_eligible", True), default=True):
            raise ValueError("Refusing --as-floors row whose target is outside the chosen AOI.")
        try:
            lat, lon = float(row["lat"]), float(row["lon"])
            aoi_name = str(row["aoi_name"])
        except (KeyError, TypeError, ValueError, OverflowError):
            raise ValueError("Refusing --as-floors row without target coordinates and an explicit AOI name.")
        if not np.isfinite(lat) or not np.isfinite(lon) or not point_in_aoi(lon, lat, aoi_name):
            raise ValueError("Refusing --as-floors row whose point is outside its named AOI.")
        match = DATE.search(str(row.get("file", "")))
        if not match:
            raise ValueError(f"Refusing --as-floors row without a valid date pair: {row.get('file')!r}")
        a, b = match.groups()
        if span_from_name(row["file"]) <= 0:
            raise ValueError(f"Refusing --as-floors row with non-positive temporal span: {row['file']}")
        validated.append((a, b, floor))

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["reference", "secondary", "layer",
                                           "detect_floor_mm_day"])
        w.writeheader()
        for a, b, floor in validated:
            w.writerow({"reference": a, "secondary": b, "layer": layer,
                        "detect_floor_mm_day": round(floor, 3)})
    logger.info("Wrote %s (%d floors)", path, len(validated))
    return len(validated)


def measure_target(files, tid, lat, lon, radii, band, *, aoi_name: str) -> list[dict]:
    """Every file, every radius, at one target. Rows carry the target."""
    import rasterio
    from pyproj import Transformer
    rows = []
    print(f"\n{'='*84}")
    print(f"DETECTION FLOOR at {lat:.5f} N {lon:.5f} E"
          + (f"  (target {tid})" if tid else "") + ", against the whole AOI")
    print("=" * 84)
    print("A floor measured over terrain that did not fail describes terrain that")
    print("did not fail. Same formula both sides: 3 * MAD-sigma / span. The")
    print("'block' column is the floor on a window MEDIAN, measured from the")
    print("scatter of same-size block medians across the AOI.\n")

    for radius in radii:
        print(f"  --- window radius {radius} px "
              f"({2*radius+1}x{2*radius+1}) ---")
        print(f"  {'PAIR':<24}{'span':>5}{'AOI':>9}{'point':>9}"
              f"{'95% CI':>12}{'block':>8}{'ratio':>7}{'valid':>11}")
        print("  " + "-" * 78)
        block = []
        for f in files:
            span = span_from_name(f.name)
            if not np.isfinite(span) or span <= 0:
                logger.warning("no date pair in %s, skipped", f.name)
                rows.append(skipped_row(f, tid, lat, lon, radius, aoi_name,
                                         "INVALID_DATE_SPAN"))
                continue
            with rasterio.open(f) as s:
                arr = s.read(band).astype(float)
                if s.nodata is not None:
                    arr[arr == s.nodata] = np.nan
                tf = Transformer.from_crs(4326, s.crs.to_epsg(), always_xy=True)
                row, col = s.index(*tf.transform(lon, lat))
            # "too few valid pixels" and "the point is not in this raster" are
            # different failures and used to print the same line. A target
            # outside the footprint is a setup error the caller must fix; a
            # target inside a nodata hole is a measurement fact about the site.
            if not (0 <= row < arr.shape[0] and 0 <= col < arr.shape[1]):
                print(f"  {DATE.search(f.name).group(0):<24}{span:>5.0f}"
                      f"   TARGET OUTSIDE RASTER - no floor")
                for radius in radii:
                    rows.append(skipped_row(f, tid, lat, lon, radius, aoi_name,
                                             "TARGET_OUTSIDE_RASTER"))
                continue
            win = window(arr, row, col, radius)
            c = compare(arr, win, span)
            name = DATE.search(f.name).group(0)
            if not c["usable"]:
                if c["status"] == "DEGENERATE_LOCAL_MAD":
                    message = "zero/non-finite local MAD; no measurable noise floor"
                elif c["status"] == "INVALID_DATE_SPAN":
                    message = "invalid date span; no floor"
                else:
                    message = (f"too few valid pixels ({c['window_px']}/"
                               f"{c['window_total_px']})")
                print(f"  {name:<24}{span:>5.0f}   {message}")
                rows.append({"file": f.name, "radius_px": radius, "target_id": tid,
                             "lat": lat, "lon": lon, "aoi_name": aoi_name,
                             "aoi_eligible": True, **c})
                continue
            ci = (f"[{c['local_floor_ci_lo']:.0f}-{c['local_floor_ci_hi']:.0f}]"
                  if np.isfinite(c["local_floor_ci_lo"]) else "-")
            print(f"  {name:<24}{span:>5.0f}{c['aoi_floor_mm_day']:>9.1f}"
                  f"{c['local_floor_mm_day']:>9.1f}{ci:>12}"
                  f"{c['block_median_floor_mm_day']:>8.1f}"
                  f"{c['ratio']:>7.2f}{c['window_px']:>7}/{c['window_total_px']}")
            block.append(c)
            rows.append({"file": f.name, "radius_px": radius, "target_id": tid,
                         "lat": lat, "lon": lon, "aoi_name": aoi_name,
                         "aoi_eligible": True, "status": "MEASURED", **c})

        if not block:
            print("  nothing usable at this radius\n")
            continue
        fa = np.array([b["aoi_floor_mm_day"] for b in block])
        fl = np.array([b["local_floor_mm_day"] for b in block])
        frac = np.array([b["window_px"] / b["window_total_px"] for b in block])
        # Geometry dominates the floor here - ascending and descending differ
        # by about 5x over this AOI - so a median across both describes no
        # real observing configuration. This module used to pool them
        # silently, which is the same trap the GOFF season analysis documents
        # and then avoids. Detect the spread and say so.
        if fa.size >= 4 and np.nanmin(fa) > 0 and np.nanmax(fa) / np.nanmin(fa) > 3.0:
            print(f"  ! AOI floors span {np.nanmin(fa):.1f} to {np.nanmax(fa):.1f} "
                  f"mm/day ({np.nanmax(fa)/np.nanmin(fa):.1f}x).")
            print(f"    That is a geometry difference, not a seasonal one.")
            print(f"    A median across both tracks describes neither -"
                  f" re-run with --include to select one.")
        print("  " + "-" * 78)
        print(f"  median AOI floor   {np.median(fa):6.1f} mm/day")
        print(f"  median point floor {np.median(fl):6.1f} mm/day"
              f"   ratio {np.median(fl)/np.median(fa):.2f}x")
        print(f"  window valid fraction: median {100*np.median(frac):.0f}%, "
              f"worst {100*frac.min():.0f}%\n")
    return rows


if __name__ == "__main__":
    sys.exit(main())
