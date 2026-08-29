"""
geometry_merge.py
-----------------
Combine ascending and descending line-of-sight time series into one denser
record of downslope motion.

The rule you cannot break
=========================
You can NEVER interfere an ascending acquisition with a descending one. The
look geometry differs, the scattering differs, coherence is zero. That is
physics.

But once each geometry has been inverted separately into its own LOS
displacement series - which is what timeseries.py produces - you are no longer
combining *phase*, you are combining *measurements*. Those can be merged, and
the merged record is much denser than any single track.

Over the Langtang AOI: five tracks (NISAR ASC 98 + DESC 48, Sentinel-1 ASC 85
+ DESC 19 + DESC 121) sample at a median of 4 days, against 12 days for any one
of them alone. Per the detectability sweep, that is the difference between
never detecting a 10-day precursor and marginally detecting one.

How the projection works
========================
InSAR measures only the component of motion along the line of sight:

    d_los = d . los_hat

With a single geometry you cannot separate subsidence from horizontal creep.
The standard fix on a hillslope is to assume motion is downslope, take that
direction from the DEM, and solve for its magnitude:

    d_downslope = d_los / (slope_hat . los_hat)

That denominator is the *sensitivity*. When it approaches zero the track is
effectively blind to the motion, and dividing by it amplifies noise without
limit. This tool computes it per track and refuses to use geometries below a
threshold - which is the whole point, because a slope facing across the look
direction can be moving fast and show nothing.

Satellite geometry
==================
Heading is derived from orbital inclination and latitude:

    sin(heading_ascending) = cos(inclination) / cos(latitude)

For a right-looking sensor the target-to-satellite unit vector in ENU is

    E = -sin(theta) cos(heading)
    N =  sin(theta) sin(heading)
    U =  cos(theta)

Check: descending Sentinel-1 at 28 N has heading ~189 deg, giving E > 0 - the
satellite is east of the target, which is correct for a west-looking descending
pass.

Usage
-----
    # Which tracks can even see motion on this slope? No data needed.
    python src/geometry_merge.py --sensitivity --lat 28.29 --lon 85.51

    # Merge the per-geometry series from timeseries.py
    python src/geometry_merge.py --merge --ts outputs/ts.csv \\
        --lat 28.29 --lon 85.51 --csv outputs/merged.csv --plot outputs/merged.png
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("geometry")


@dataclass(frozen=True)
class Track:
    name: str
    mission: str
    ascending: bool
    inclination_deg: float
    incidence_deg: float
    repeat_days: int


# Tracks covering the Langtang / Lhende AOIs. Incidence angles are nominal
# mid-swath values - read the incidenceAngle layer from the GUNW for the real
# per-pixel number when you have the products.
TRACKS = [
    Track("NISAR ASC 98", "NISAR", True, 98.4, 37.0, 12),
    Track("NISAR DESC 48", "NISAR", False, 98.4, 37.0, 12),
    Track("S1 ASC 85", "Sentinel-1", True, 98.18, 39.0, 12),
    Track("S1 DESC 19", "Sentinel-1", False, 98.18, 39.0, 12),
    Track("S1 DESC 121", "Sentinel-1", False, 98.18, 39.0, 12),
]


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
def heading_deg(inclination_deg: float, lat_deg: float, ascending: bool) -> float:
    """
    Satellite heading (azimuth of flight, clockwise from North).

    sin(heading_asc) = cos(inclination) / cos(latitude)
    """
    ratio = math.cos(math.radians(inclination_deg)) / math.cos(math.radians(lat_deg))
    ratio = max(-1.0, min(1.0, ratio))
    asc = math.degrees(math.asin(ratio))           # negative for a retrograde orbit
    return asc % 360.0 if ascending else (180.0 - asc) % 360.0


def los_unit(heading: float, incidence: float) -> np.ndarray:
    """Unit vector from the ground target toward the satellite, in ENU."""
    h, t = math.radians(heading), math.radians(incidence)
    return np.array([-math.sin(t) * math.cos(h),
                      math.sin(t) * math.sin(h),
                      math.cos(t)])


def downslope_unit(slope_deg: float, aspect_deg: float) -> np.ndarray:
    """
    Unit vector pointing down the slope, in ENU.

    slope 0 deg  -> horizontal, along aspect
    slope 90 deg -> straight down
    """
    b, a = math.radians(slope_deg), math.radians(aspect_deg)
    return np.array([math.cos(b) * math.sin(a),
                     math.cos(b) * math.cos(a),
                     -math.sin(b)])


# ---------------------------------------------------------------------------
# Terrain
# ---------------------------------------------------------------------------
def slope_aspect(lat: float, lon: float, spacing_m: float = 90.0) -> tuple[float, float, float]:
    """Slope, aspect and elevation at a point, from SRTM via OpenTopoData."""
    dlat = spacing_m / 110570.0
    dlon = spacing_m / (111320.0 * math.cos(math.radians(lat)))
    pts = [(lat + i * dlat, lon + j * dlon) for i in (1, 0, -1) for j in (-1, 0, 1)]
    loc = "|".join(f"{a:.6f},{b:.6f}" for a, b in pts)
    url = "https://api.opentopodata.org/v1/srtm30m?locations=" + urllib.parse.quote(loc)
    for attempt in range(4):
        try:
            d = json.loads(urllib.request.urlopen(url, timeout=90).read().decode())
            break
        except Exception as exc:
            if attempt == 3:
                raise
            logger.warning("retry (%s)", exc)
            time.sleep(3)
    z = np.array([r["elevation"] for r in d["results"]], dtype=float).reshape(3, 3)

    # Horn's method. Row 0 is north, so dz/dy uses (north - south).
    dzdx = ((z[0, 2] + 2 * z[1, 2] + z[2, 2]) - (z[0, 0] + 2 * z[1, 0] + z[2, 0])) / (8 * spacing_m)
    dzdy = ((z[0, 0] + 2 * z[0, 1] + z[0, 2]) - (z[2, 0] + 2 * z[2, 1] + z[2, 2])) / (8 * spacing_m)
    slope = math.degrees(math.atan(math.hypot(dzdx, dzdy)))
    aspect = (math.degrees(math.atan2(-dzdx, -dzdy))) % 360.0   # azimuth of steepest DESCENT
    return slope, aspect, float(z[1, 1])


# ---------------------------------------------------------------------------
def sensitivities(lat: float, lon: float, slope: float, aspect: float,
                  min_sensitivity: float) -> list[dict]:
    s_hat = downslope_unit(slope, aspect)
    rows = []
    for tr in TRACKS:
        h = heading_deg(tr.inclination_deg, lat, tr.ascending)
        l_hat = los_unit(h, tr.incidence_deg)
        sens = float(np.dot(s_hat, l_hat))
        rows.append({
            "track": tr.name, "ascending": tr.ascending,
            "heading_deg": h, "incidence_deg": tr.incidence_deg,
            "sensitivity": sens,
            "amplification": abs(1.0 / sens) if abs(sens) > 1e-6 else float("inf"),
            "usable": abs(sens) >= min_sensitivity,
            "repeat_days": tr.repeat_days,
        })
    return rows


def report_sensitivity(lat, lon, slope, aspect, elev, rows, min_sensitivity):
    print(f"\nPoint {lat:.4f} N, {lon:.4f} E   elevation {elev:.0f} m")
    print(f"Slope {slope:.1f} deg, aspect {aspect:.0f} deg "
          f"({compass(aspect)}-facing) - from SRTM at 90 m spacing")
    print(f"\nAssuming pure downslope motion. Sensitivity = slope_hat . los_hat;")
    print(f"a track needs |sensitivity| >= {min_sensitivity} to be usable.\n")
    print(f"{'TRACK':<16}{'HEADING':>9}{'INCID':>7}{'SENSITIVITY':>13}{'NOISE x':>9}  VERDICT")
    print("-" * 68)
    for r in sorted(rows, key=lambda x: -abs(x["sensitivity"])):
        amp = f"{r['amplification']:.1f}" if math.isfinite(r["amplification"]) else "inf"
        verdict = "usable" if r["usable"] else "BLIND - reject"
        print(f"{r['track']:<16}{r['heading_deg']:>8.1f}d{r['incidence_deg']:>6.0f}d"
              f"{r['sensitivity']:>13.3f}{amp:>9}  {verdict}")

    usable = [r for r in rows if r["usable"]]
    print(f"\n  {len(usable)} of {len(rows)} tracks usable on this slope.")
    if usable:
        eff = harmonic_revisit([r["repeat_days"] for r in usable])
        print(f"  Combined sampling from the usable tracks: ~{eff:.1f} days "
              f"(vs {min(r['repeat_days'] for r in usable)} d for any one alone).")
    signs = {r["sensitivity"] > 0 for r in usable}
    if len(signs) > 1:
        print("  Ascending and descending have OPPOSITE sign here - a genuine")
        print("  cross-check: real downslope motion must appear with opposite")
        print("  LOS sign in the two geometries. Atmosphere will not do that.")
    else:
        print("  All usable tracks share the same sign - no sign-based cross-check")
        print("  available at this point.")


def compass(azimuth: float) -> str:
    names = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
             "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]
    return names[int((azimuth % 360) / 22.5 + 0.5) % 16]


def harmonic_revisit(repeats: list[int]) -> float:
    """Mean gap when several independent tracks are interleaved."""
    return 1.0 / sum(1.0 / r for r in repeats)


# ---------------------------------------------------------------------------
# Merge
# ---------------------------------------------------------------------------
def load_timeseries(path: Path) -> dict[str, list[dict]]:
    series: dict[str, list[dict]] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            key = f"{row['geometry']}|{row['component']}"
            series.setdefault(key, []).append({
                "epoch": datetime.strptime(row["epoch"], "%Y-%m-%d").date(),
                "cumulative_mm": float(row["cumulative_mm"]),
                "error_mm": float(row["error_mm"]) if row.get("error_mm") else None,
            })
    for v in series.values():
        v.sort(key=lambda r: r["epoch"])
    return series


def match_track(geometry_label: str) -> Track | None:
    """Map a timeseries.py label like 'ASC path 98' onto a Track."""
    lab = geometry_label.upper().replace("PATH", "").split()
    asc = lab[0].startswith("ASC")
    num = next((t for t in lab if t.isdigit()), None)
    for tr in TRACKS:
        if tr.ascending == asc and num and num in tr.name:
            return tr
    return None


def merge(ts_path: Path, lat: float, lon: float, min_sensitivity: float):
    slope, aspect, elev = slope_aspect(lat, lon)
    s_hat = downslope_unit(slope, aspect)
    series = load_timeseries(ts_path)

    print(f"\nProjecting onto downslope at {lat:.4f} N {lon:.4f} E "
          f"(slope {slope:.1f} deg, aspect {aspect:.0f} deg / {compass(aspect)})\n")

    merged: list[dict] = []
    for key, rows in sorted(series.items()):
        geom, comp = key.split("|")
        tr = match_track(geom)
        if tr is None:
            logger.warning("No geometry known for '%s' - skipped", geom)
            continue
        h = heading_deg(tr.inclination_deg, lat, tr.ascending)
        sens = float(np.dot(s_hat, los_unit(h, tr.incidence_deg)))
        if abs(sens) < min_sensitivity:
            print(f"  {geom} block {comp}: sensitivity {sens:+.3f} - BLIND, rejected")
            continue
        print(f"  {geom} block {comp}: sensitivity {sens:+.3f}, "
              f"noise amplified x{abs(1/sens):.1f}, {len(rows)} epochs")
        for r in rows:
            merged.append({
                "epoch": r["epoch"], "geometry": geom, "component": int(comp),
                "los_mm": r["cumulative_mm"],
                "downslope_mm": r["cumulative_mm"] / sens,
                "sensitivity": sens,
                "weight": abs(sens),
                "error_mm": (abs(r["error_mm"] / sens) if r["error_mm"] else None),
            })

    merged.sort(key=lambda r: r["epoch"])
    if not merged:
        print("\nNothing usable - every geometry was blind to motion on this slope.")
        return []

    epochs = [r["epoch"] for r in merged]
    gaps = [(epochs[i + 1] - epochs[i]).days for i in range(len(epochs) - 1)]
    gaps = [g for g in gaps if g > 0]

    print(f"\n{'EPOCH':<12}{'GEOMETRY':<14}{'BLK':>4}{'LOS mm':>10}{'DOWNSLOPE mm':>14}{'SENS':>8}")
    print("-" * 64)
    for r in merged:
        print(f"{str(r['epoch']):<12}{r['geometry']:<14}{r['component']:>4}"
              f"{r['los_mm']:>10.2f}{r['downslope_mm']:>14.2f}{r['sensitivity']:>8.3f}")

    if gaps:
        print(f"\n  merged sampling: {len(merged)} epochs, "
              f"median gap {float(np.median(gaps)):.1f} d, min {min(gaps)} d, max {max(gaps)} d")
    print("\n  NOTE: blocks keep their own zero. Points from different components")
    print("  or different geometries are NOT on a common datum - the merged record")
    print("  is denser in TIME, but only comparable within a block.")
    return merged


def plot_merged(rows, out: Path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.error("matplotlib required for --plot")
        return
    fig, ax = plt.subplots(figsize=(10, 6))
    by_geom: dict[str, list] = {}
    for r in rows:
        by_geom.setdefault(f"{r['geometry']} blk{r['component']}", []).append(r)
    for label, rs in sorted(by_geom.items()):
        ax.plot([datetime.combine(r["epoch"], datetime.min.time()) for r in rs],
                [r["downslope_mm"] for r in rs], "o-", label=label, alpha=.85)
    ax.axvline(datetime(2026, 8, 26), color="crimson", ls="--", lw=1.2)
    ax.set_ylabel("Downslope displacement (mm)")
    ax.set_title("Multi-geometry merge: LOS projected onto slope direction\n"
                 "each block keeps its own zero", fontsize=11)
    ax.grid(alpha=.3); ax.legend(fontsize=9)
    fig.autofmt_xdate(); fig.tight_layout(); fig.savefig(out, dpi=140)
    logger.info("Wrote %s", out)


def main() -> int:
    ap = argparse.ArgumentParser(description="Merge ascending and descending LOS series")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--sensitivity", action="store_true",
                      help="which tracks can see motion here (needs no data)")
    mode.add_argument("--merge", action="store_true", help="merge a timeseries.py CSV")
    ap.add_argument("--lat", type=float, required=True)
    ap.add_argument("--lon", type=float, required=True)
    ap.add_argument("--ts", metavar="TS.csv", help="output of timeseries.py --csv")
    ap.add_argument("--min-sensitivity", type=float, default=0.3,
                    help="reject tracks below this |slope_hat . los_hat|")
    ap.add_argument("--csv", metavar="OUT.csv")
    ap.add_argument("--plot", metavar="OUT.png")
    args = ap.parse_args()

    if args.sensitivity:
        slope, aspect, elev = slope_aspect(args.lat, args.lon)
        rows = sensitivities(args.lat, args.lon, slope, aspect, args.min_sensitivity)
        report_sensitivity(args.lat, args.lon, slope, aspect, elev, rows, args.min_sensitivity)
        if args.csv:
            with open(args.csv, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]))
                w.writeheader(); w.writerows(rows)
            logger.info("Wrote %s", args.csv)
        return 0

    if not args.ts:
        ap.error("--merge needs --ts")
    rows = merge(Path(args.ts), args.lat, args.lon, args.min_sensitivity)
    if rows and args.csv:
        keys = ["epoch", "geometry", "component", "los_mm", "downslope_mm",
                "sensitivity", "weight", "error_mm"]
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow({**r, "epoch": str(r["epoch"])})
        logger.info("Wrote %s", args.csv)
    if rows and args.plot:
        plot_merged(rows, Path(args.plot))
    return 0


if __name__ == "__main__":
    sys.exit(main())
