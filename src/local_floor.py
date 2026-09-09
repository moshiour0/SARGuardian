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
a window is set by its weakest interval, so **40.4 mm/day** is the figure - and
`40.4 [21-55]` once the bootstrap interval is attached, because it rests on 49
valid pixels. Expressed as downslope motion rather than line of sight it
weakens further, to as much as 143 mm/day, because the aspect at this point is
unconfirmed; see the README.

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
    python src/local_floor.py --dir outputs/export_goff_src --match layer2 \\
        --lat 28.28771 --lon 85.52809 --radius 6

    python src/local_floor.py --dir outputs/export_goff_src --match layer2 \\
        --lat 28.28771 --lon 85.52809 --sweep 3 6 12 --csv outputs/local_floor.csv
"""

from __future__ import annotations

import argparse
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
    if span_days <= 0:
        return float("nan")
    return 3.0 * robust_sigma(values_mm) / span_days


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
    Reported as such rather than corrected, because the correlation length of
    the offset field is not measured here.
    """
    v = values_mm[np.isfinite(values_mm)]
    n = v.size
    if n < MIN_PX or span_days <= 0:
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


def compare(aoi_mm: np.ndarray, win_mm: np.ndarray, span_days: float) -> dict:
    """
    Both floors and the ratio between them, from one pair.

    Pure: takes arrays, returns numbers, no raster and no I/O, so the claim
    this module makes can be tested without a product.
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
           "ratio": float("nan"), "usable": False}
    if a.size < MIN_PX or w.size < MIN_PX:
        return out
    fa = detection_floor(a, span_days)
    fl = detection_floor(w, span_days)
    lo, hi = floor_ci(w, span_days)
    out.update(aoi_floor_mm_day=fa, local_floor_mm_day=fl,
               local_floor_ci_lo=lo, local_floor_ci_hi=hi,
               ratio=(fl / fa) if fa > 0 else float("nan"), usable=True)
    return out


def span_from_name(name: str) -> float:
    """Days between the two YYYYMMDD stamps in an export filename."""
    m = DATE.search(name)
    if not m:
        return float("nan")
    a, b = (dt.date(int(s[:4]), int(s[4:6]), int(s[6:8])) for s in m.groups())
    return float((b - a).days)


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description="Detection floor at a point against the floor over the AOI")
    ap.add_argument("--dir", required=True, help="directory of exported GeoTIFFs")
    ap.add_argument("--match", metavar="SUBSTR",
                    help="only files whose name contains this, e.g. layer2")
    ap.add_argument("--lat", type=float, required=True)
    ap.add_argument("--lon", type=float, required=True)
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

    radii = args.sweep or [args.radius]
    rows = []
    print(f"\n{'='*84}")
    print(f"DETECTION FLOOR at {args.lat:.5f} N {args.lon:.5f} E, against the whole AOI")
    print("=" * 84)
    print("A floor measured over terrain that did not fail describes terrain that")
    print("did not fail. Same formula both sides: 3 * MAD-sigma / span.\n")

    for radius in radii:
        print(f"  --- window radius {radius} px "
              f"({2*radius+1}x{2*radius+1}) ---")
        print(f"  {'PAIR':<24}{'span':>5}{'AOI':>9}{'point':>9}"
              f"{'95% CI':>12}{'ratio':>7}{'valid':>11}")
        print("  " + "-" * 78)
        block = []
        for f in files:
            span = span_from_name(f.name)
            if not np.isfinite(span) or span <= 0:
                logger.warning("no date pair in %s, skipped", f.name)
                continue
            with rasterio.open(f) as s:
                arr = s.read(args.band).astype(float)
                if s.nodata is not None:
                    arr[arr == s.nodata] = np.nan
                tf = Transformer.from_crs(4326, s.crs.to_epsg(), always_xy=True)
                row, col = s.index(*tf.transform(args.lon, args.lat))
            # "too few valid pixels" and "the point is not in this raster" are
            # different failures and used to print the same line. A target
            # outside the footprint is a setup error the caller must fix; a
            # target inside a nodata hole is a measurement fact about the site.
            if not (0 <= row < arr.shape[0] and 0 <= col < arr.shape[1]):
                print(f"  {DATE.search(f.name).group(0):<24}{span:>5.0f}"
                      f"   TARGET OUTSIDE RASTER - check --lat/--lon")
                continue
            win = window(arr, row, col, radius)
            c = compare(arr, win, span)
            name = DATE.search(f.name).group(0)
            if not c["usable"]:
                print(f"  {name:<24}{span:>5.0f}   too few valid pixels "
                      f"({c['window_px']}/{c['window_total_px']})")
                continue
            ci = (f"[{c['local_floor_ci_lo']:.0f}-{c['local_floor_ci_hi']:.0f}]"
                  if np.isfinite(c["local_floor_ci_lo"]) else "-")
            print(f"  {name:<24}{span:>5.0f}{c['aoi_floor_mm_day']:>9.1f}"
                  f"{c['local_floor_mm_day']:>9.1f}{ci:>12}{c['ratio']:>7.2f}"
                  f"{c['window_px']:>7}/{c['window_total_px']}")
            block.append(c)
            rows.append({"file": f.name, "radius_px": radius, **c})

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

    if not rows:
        print("  Nothing measured.")
        return 1

    r = np.array([x["ratio"] for x in rows])
    if np.median(r) > 1.25:
        print("  The floor at the point is materially WORSE than the floor over")
        print("  the AOI. Any bound quoted at this location must use the local")
        print("  figure; the AOI figure overstates what the product could see.")
    elif np.median(r) < 0.8:
        print("  The point is quieter than the AOI. The AOI floor is")
        print("  conservative here, which is safe but not tight.")
    else:
        print("  The two agree. The AOI floor is representative at this point.")

    if args.csv:
        import csv as _csv
        p = Path(resolve(args.csv))
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "w", newline="") as fh:
            w = _csv.DictWriter(fh, fieldnames=sorted(rows[0]))
            w.writeheader()
            for x in rows:
                w.writerow(x)
        logger.info("Wrote %s", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
