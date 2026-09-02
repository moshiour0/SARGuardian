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


def los_unit(heading: float, incidence: float, left_looking: bool = True) -> np.ndarray:
    """
    Unit vector from the ground target toward the satellite, in ENU.

    The look side matters and is easy to get wrong. This formula was written
    for a right-looking sensor; NISAR looks LEFT, and its products say so in
    identification/lookDirection. Using the right-looking form on NISAR
    reverses both horizontal components while leaving the vertical untouched -
    so an AOI median, which is dominated by the vertical term, looks entirely
    reasonable while every east-west inference is backwards.

    Checked against the products' own losUnitVectorX/Y at 28.275 N:

        ASC 98   derived right-looking  E -0.6261  N -0.1053  U +0.7726
                 product                E +0.6161  N +0.1531  U +0.7727

    Prefer los_from_product() when a granule is at hand; this is the fallback
    for planning before anything is downloaded.
    """
    h, t = math.radians(heading), math.radians(incidence)
    side = -1.0 if left_looking else 1.0
    return np.array([side * -math.sin(t) * math.cos(h),
                     side * math.sin(t) * math.sin(h),
                     math.cos(t)])


def los_from_product(path, lat: float, lon: float) -> dict | None:
    """
    Read the look geometry the product itself carries, at one location.

    NISAR L2 ships losUnitVectorX/Y and incidenceAngle on a metadata cube, so
    the true geometry never has to be reconstructed from orbital elements at
    all. Z follows from the unit-length constraint, and the result is checked
    against cos(incidence) before being returned - if those disagree the
    convention is not what we think it is and guessing would be worse than
    failing.
    """
    try:
        import h5py
        from pyproj import Transformer
    except ImportError:
        logger.warning("h5py and pyproj are needed to read look geometry.")
        return None

    rg = "science/LSAR/GUNW/metadata/radarGrid"
    with h5py.File(path, "r") as f:
        if f"{rg}/losUnitVectorX" not in f:
            return None
        xs = np.asarray(f[f"{rg}/xCoordinates"][()])
        ys = np.asarray(f[f"{rg}/yCoordinates"][()])
        epsg = int(np.ravel(f[f"{rg}/projection"][()])[0]) if f"{rg}/projection" in f else 4326
        tx, ty = Transformer.from_crs(4326, epsg, always_xy=True).transform(lon, lat)
        j = int(np.argmin(np.abs(xs - tx)))
        i = int(np.argmin(np.abs(ys - ty)))

        def col(name):
            a = np.asarray(f[f"{rg}/{name}"][:, i, j], dtype=float)
            a = a[np.isfinite(a)]
            return float(a.mean()) if a.size else float("nan")

        lx, ly, inc = col("losUnitVectorX"), col("losUnitVectorY"), col("incidenceAngle")
        look = f["science/LSAR/identification/lookDirection"][()]
        look = look.decode() if isinstance(look, bytes) else str(look)

    if not all(np.isfinite([lx, ly, inc])):
        return None
    lz = math.sqrt(max(0.0, 1.0 - lx * lx - ly * ly))
    if abs(lz - math.cos(math.radians(inc))) > 0.02:
        logger.warning("losUnitVector Z (%.4f) disagrees with cos(incidence) "
                       "(%.4f) - convention unclear, not using it.",
                       lz, math.cos(math.radians(inc)))
        return None
    return {"los": np.array([lx, ly, lz]), "incidence": inc, "look": look}


def decompose(d_asc: float, d_desc: float, los_asc: np.ndarray,
              los_desc: np.ndarray) -> dict:
    """
    Two line-of-sight rates into east and vertical, assuming no north motion.

    Both geometries are nearly blind to north - the N components here are
    +0.15 and +0.18 against E of +0.62 and -0.59 - so solving for it would
    divide by almost nothing. Setting u_N = 0 is the standard reduction and
    the honest one.

    The reason to bother is that a path delay is not a vector. A delay changes
    both geometries by the same LOS amount, and because U is the largest
    component of both look vectors (0.77 and 0.79) it decomposes to almost
    pure vertical. So the vertical fraction is a diagnostic: a slope creeps
    parallel to its own surface, and a plunge steeper than the terrain is a
    delay wearing a vector's clothes.
    """
    M = np.array([[los_asc[0], los_asc[2]], [los_desc[0], los_desc[2]]])
    if abs(np.linalg.det(M)) < 1e-6:
        raise ValueError("The two geometries are too similar to separate.")
    u_e, u_u = np.linalg.solve(M, [d_asc, d_desc])
    mag = float(math.hypot(u_e, u_u))
    return {"east": float(u_e), "up": float(u_u), "magnitude": mag,
            "plunge_deg": float(math.degrees(math.atan2(-u_u, abs(u_e)))) if mag else 0.0,
            "vertical_fraction": float(abs(u_u) / mag) if mag else 0.0}


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
def elevation_grid(lats, lons, batch: int = 100):
    """
    One elevation raster for the whole AOI, in batched requests.

    slope_aspect() costs a request per point, which is fine for one location
    and hopeless for a map. OpenTopoData takes 100 locations per call, so a
    40x36 grid is fifteen requests instead of fourteen hundred.
    """
    pts = [(a, b) for a in lats for b in lons]
    out = []
    for k in range(0, len(pts), batch):
        chunk = pts[k:k + batch]
        loc = "|".join(f"{a:.6f},{b:.6f}" for a, b in chunk)
        url = "https://api.opentopodata.org/v1/srtm30m?locations=" + urllib.parse.quote(loc)
        for attempt in range(4):
            try:
                d = json.loads(urllib.request.urlopen(url, timeout=120).read().decode())
                break
            except Exception as exc:
                if attempt == 3:
                    raise
                logger.warning("retry (%s)", exc)
                time.sleep(4)
        out += [r["elevation"] if r["elevation"] is not None else np.nan
                for r in d["results"]]
        time.sleep(1.1)                      # the public endpoint allows 1/sec
        logger.info("elevation %d/%d", min(k + batch, len(pts)), len(pts))
    return np.array(out, dtype=float).reshape(len(lats), len(lons))


def slope_aspect_grid(z: np.ndarray, lats, lons):
    """Horn's method across a whole grid. Row 0 is north."""
    dy = abs(lats[1] - lats[0]) * 110570.0
    dx = abs(lons[1] - lons[0]) * 111320.0 * math.cos(math.radians(float(np.mean(lats))))
    zp = np.pad(z, 1, mode="edge")
    dzdx = ((zp[:-2, 2:] + 2 * zp[1:-1, 2:] + zp[2:, 2:])
            - (zp[:-2, :-2] + 2 * zp[1:-1, :-2] + zp[2:, :-2])) / (8 * dx)
    dzdy = ((zp[:-2, :-2] + 2 * zp[:-2, 1:-1] + zp[:-2, 2:])
            - (zp[2:, :-2] + 2 * zp[2:, 1:-1] + zp[2:, 2:])) / (8 * dy)
    slope = np.degrees(np.arctan(np.hypot(dzdx, dzdy)))
    aspect = np.degrees(np.arctan2(-dzdx, -dzdy)) % 360.0
    return slope, aspect


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


def coverage_map(aoi: str, step: float, min_sens: float, out_csv=None) -> int:
    """
    How much of an AOI could each geometry actually have seen?

    A non-detection is only as strong as the ground it covers. Quoting one
    sensitivity at one point says nothing about a valley whose aspect swings
    through every direction, and this is the number that qualifies the null:
    motion is invisible wherever the slope happens to lie across the look
    direction, however good the interferogram is.

    It also separates two different questions. One usable geometry gives a
    line-of-sight rate and nothing more; only where BOTH are usable can
    vertical be told from horizontal, and that is a much smaller area.
    """
    import gunw_reader as _g
    ring = _g.AOIS[aoi]
    lons = [p[0] for p in ring]
    lats = [p[1] for p in ring]
    LO = np.arange(min(lons), max(lons) + 1e-9, step)
    LA = np.arange(max(lats), min(lats) - 1e-9, -step)
    logger.info("%s: %d x %d grid at ~%.0f m", aoi, len(LA), len(LO), step * 111000)

    z = elevation_grid(LA, LO)
    slope, aspect = slope_aspect_grid(z, LA, LO)
    inside = np.array([[_g.point_in_ring(b, a, ring) for b in LO] for a in LA])
    ok = inside & np.isfinite(z)

    print(f"\n{'=' * 68}")
    print(f"{aoi.upper()}   {int(ok.sum())} cells inside the AOI")
    print("=" * 68)
    print(f"  elevation      {np.nanmin(z[ok]):.0f} - {np.nanmax(z[ok]):.0f} m")
    print(f"  slope          median {np.median(slope[ok]):.1f} deg, "
          f"{100 * np.mean(slope[ok] > 20):.0f}% steeper than 20 deg")

    sens = {}
    for label, ascending, inc in (("NISAR ASC 98", True, 39.4),
                                  ("NISAR DESC 48", False, 37.7)):
        h = heading_deg(98.4, float(np.mean(LA)), ascending)
        L = los_unit(h, inc)
        S = np.zeros_like(slope)
        for i in range(slope.shape[0]):
            for j in range(slope.shape[1]):
                if ok[i, j]:
                    S[i, j] = float(np.dot(downslope_unit(slope[i, j], aspect[i, j]), L))
        sens[label] = S
        steep = ok & (slope > 20)
        print(f"  {label:<15} usable over {100 * (np.abs(S[ok]) >= min_sens).mean():5.1f}% "
              f"of the AOI, {100 * (np.abs(S[steep]) >= min_sens).mean():5.1f}% of "
              f"slopes >20 deg")

    A, D = sens["NISAR ASC 98"], sens["NISAR DESC 48"]
    either = (np.abs(A) >= min_sens) | (np.abs(D) >= min_sens)
    both = (np.abs(A) >= min_sens) & (np.abs(D) >= min_sens)
    print(f"\n  EITHER geometry usable   {100 * either[ok].mean():5.1f}%  "
          f"<- the area a non-detection actually covers")
    print(f"  BOTH usable              {100 * both[ok].mean():5.1f}%  "
          f"<- the only area where vertical can be separated from horizontal")
    print(f"  NEITHER                  {100 * (~either[ok]).mean():5.1f}%  "
          f"<- blind, whatever the interferogram says")

    if out_csv:
        with open(out_csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["lat", "lon", "elev_m", "slope_deg", "aspect_deg",
                        "sens_asc", "sens_desc"])
            for i, a in enumerate(LA):
                for j, b in enumerate(LO):
                    if ok[i, j]:
                        w.writerow([f"{a:.5f}", f"{b:.5f}", f"{z[i,j]:.0f}",
                                    f"{slope[i,j]:.1f}", f"{aspect[i,j]:.0f}",
                                    f"{A[i,j]:.3f}", f"{D[i,j]:.3f}"])
        logger.info("Wrote %s", out_csv)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Merge ascending and descending LOS series")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--sensitivity", action="store_true",
                      help="which tracks can see motion here (needs no data)")
    mode.add_argument("--merge", action="store_true", help="merge a timeseries.py CSV")
    mode.add_argument("--map", metavar="AOI", choices=("langtang", "lhende"),
                      help="what fraction of an AOI each geometry can actually see")
    ap.add_argument("--lat", type=float)
    ap.add_argument("--lon", type=float)
    ap.add_argument("--step", type=float, default=0.003,
                    help="map grid spacing in degrees (0.003 ~ 333 m)")
    ap.add_argument("--ts", metavar="TS.csv", help="output of timeseries.py --csv")
    ap.add_argument("--min-sensitivity", type=float, default=0.3,
                    help="reject tracks below this |slope_hat . los_hat|")
    ap.add_argument("--csv", metavar="OUT.csv")
    ap.add_argument("--plot", metavar="OUT.png")
    args = ap.parse_args()

    if args.map:
        return coverage_map(args.map, args.step, args.min_sensitivity, args.csv)

    if args.sensitivity:
        if args.lat is None or args.lon is None:
            ap.error("--sensitivity needs --lat and --lon")
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
