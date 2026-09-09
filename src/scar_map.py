"""
scar_map.py
-----------
Map the detachment scar from Sentinel-2, and settle the aspect the SAR bound
depends on.

Why this exists
===============
Every sensitivity number in this project is a dot product between the radar
line of sight and the downslope vector, so the ASPECT of the failing surface
decides which tracks can see it and therefore what the bound means. Until this
module existed that aspect came from an SRTM pixel at a failure point which was
itself a guess, read off a report.

The guess was wrong. SRTM reads west-facing (273 deg) at 28.28771 N 85.52809 E,
while every published account puts the scar on the NORTH face of Langtang
Lirung. Both cannot be true, and the difference is a factor of three in the
bound:

    assumed surface        NISAR ASC 098      downslope bound
    west-facing 273 deg    -0.892 usable       45 mm/day
    due north     0 deg    -0.283 BLIND       143 mm/day

So it had to be measured. This does that from optical imagery, which sees the
scar directly rather than inferring it.

Method
======
1. Query the free Element84 Earth Search STAC for Sentinel-2 L2A over the AOI.
2. Score every scene by cloud over the AOI, not over the 110 km tile. They are
   very different numbers here: 2026-09-08 is 62.7% cloudy as a scene and 74.2%
   over the source zone, while 2026-08-12 is 18.7% and 2.6%.
3. Take the clearest pre-event scene, and composite the post-event scenes,
   most recent clear pixel winning. One post-event scene covers a quarter of
   the AOI; seven of them stacked reach half.
4. Difference the Normalised Difference Snow Index. The scar is ice and snow
   replaced by bare rock, so NDSI falls hard. Fresh September snowfall moves it
   the OTHER way, which is what makes this robust: weather cannot manufacture
   the signal, it can only hide it.
5. Region-grow from the strongest change, requiring each pixel to have been
   snow or ice beforehand.

What it found
=============
A feature of 0.68 km2 spanning 1,390 m north-south and 1,218 m east-west,
centred on 28.27802 N 85.52963 E:

    was snow or ice before      99.8%
    median NDSI change          -0.374
    scar-like pixel rate        90% inside, against 7.5% across the AOI
    elevation                   5,203 - 5,730 m
    slope                       27 - 52 deg, mean 39
    aspect                      298 - 354 deg, circular mean 342

Four independent checks agree with the published accounts, none of which went
into finding it: the detachment is described at about 5,200 m (measured 5,203
at its northern end), on the north face (measured 342 deg), on a rock face
under a hanging glacier (measured 99.8% ice-covered before), and roughly 1.4 km
across (measured 1.39 km).

The assumed failure point is 1.09 km away and is not on it. That location is
fully observed in both epochs - 121 usable pixels, no cloud to hide behind -
and its NDSI went UP by 0.379. Not one pixel fell. It was never glaciated. It
is bare ground that got a dusting of fresh snow, which is the opposite of a
scar in every respect.

What this changes
=================
At the mapped scar NISAR ASC 098 has a sensitivity of -0.47 at the centroid and
-0.57 on the mean aspect: usable, but two to four times noisier than the -0.892
the west-facing guess implied. NISAR descending is blind either way.

The line-of-sight floor also has to be re-measured, because the old one was
taken at the wrong place. Over the seven weeks before failure the three
covering intervals give 32.0, 26.0 and 7.8 mm/day AT THE SCAR, against 40.4,
19.2 and 34.3 at the abandoned point. A bound holding across a window is set by
its weakest interval, so the line-of-sight bound is 32.0 mm/day and the
downslope bound is

    32.0 / 0.472  =  68 mm/day

between the 45 the west-facing assumption licensed and the 143 the due-north
reading forced, and now resting on a surface that was measured rather than
guessed. The pixel counts are healthier too: 71, 60 and 60 of 169 at the scar
against 49, 65 and 26 at the old point.

The conclusion is untouched. The precursor Sentinel-1 recorded was about 0.33
mm/day, which is 200x below 68.

Limits, which are real
======================
Half the AOI is never seen cloud-free after the event, so this is the largest
feature that can be mapped with confidence, not provably the only one. The
comparison spans 12 August to 8 September, so a fortnight of ordinary seasonal
change sits inside it, and 7.5% of comparable ground looks scar-like on the
same test. The identification rests on the enrichment, the terrain and the
extent agreeing together - not on the change image alone. There is no published
scar polygon to check it against.

Usage
-----
    python src/scar_map.py --survey                 # what imagery exists
    python src/scar_map.py --map                    # find and measure the scar
    python src/scar_map.py --probe 28.28771 85.52809   # test one location
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import os
import sys
import urllib.request
from collections import deque
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("scar_map")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import resolve  # noqa: E402

STAC = "https://earth-search.aws.element84.com/v1/search"
COLLECTION = "sentinel-2-l2a"

# Source-zone bounding box, matching gunw_reader.SOURCE_RING.
AOI = (85.4645, 28.2453, 85.5562, 28.3529)

# Sentinel-2 scene classification. Cloud, cloud shadow, cirrus and no-data are
# all unusable; snow (11) emphatically is not, since snow is the thing that
# disappears.
UNUSABLE = {0, 1, 3, 8, 9, 10}
SNOW = 11

# Grid the analysis runs on: 20 m, the native posting of SWIR16.
GRID = (602, 457)

EVENT = "2026-08-26"


def _stac(body: dict) -> list[dict]:
    req = urllib.request.Request(STAC, data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as fh:
        return json.load(fh).get("features", [])


def search(start: str, end: str, limit: int = 100) -> list[dict]:
    """Every L2A scene intersecting the AOI in a date range, oldest first."""
    feats = _stac({"collections": [COLLECTION], "bbox": list(AOI),
                   "datetime": f"{start}T00:00:00Z/{end}T23:59:59Z",
                   "limit": limit})
    return sorted(feats, key=lambda f: f["properties"]["datetime"])


def by_id(ids: list[str]) -> dict[str, dict]:
    return {f["id"]: f for f in _stac({"collections": [COLLECTION],
                                       "ids": ids, "limit": len(ids) + 5})}


def _read(feature: dict, asset: str) -> np.ndarray:
    """One asset, clipped to the AOI and resampled onto the common grid."""
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.warp import transform_bounds
    from rasterio.windows import from_bounds
    href = feature["assets"][asset]["href"]
    with rasterio.open("/vsicurl/" + href) as s:
        b = transform_bounds("EPSG:4326", s.crs, *AOI)
        w = from_bounds(*b, transform=s.transform)
        return s.read(1, window=w, out_shape=GRID,
                      resampling=Resampling.nearest).astype(np.float32)


def grid_transform():
    """
    Affine for the common grid.

    Taken from the AOI bounds and the grid shape, NOT from a source window
    transform: the bands have different native postings (green 10 m, swir16
    20 m) and are resampled onto one grid here, so a transform lifted from the
    10 m band would be wrong by a factor of two and place the scar 600 m off.
    """
    from rasterio.transform import from_bounds as tf_from_bounds
    from rasterio.warp import transform_bounds
    b = transform_bounds("EPSG:4326", "EPSG:32645", *AOI)
    return tf_from_bounds(*b, width=GRID[1], height=GRID[0])


def aoi_cloud(feature: dict) -> dict:
    """Cloud over the AOI, which is not the scene-level figure."""
    scl = _read(feature, "scl").astype(int)
    bad = np.isin(scl, list(UNUSABLE))
    return {"id": feature["id"], "date": feature["properties"]["datetime"][:10],
            "scene_cloud": feature["properties"].get("eo:cloud_cover", float("nan")),
            "aoi_cloud_pct": 100.0 * bad.mean(),
            "aoi_clear_pct": 100.0 * (~bad).mean(),
            "snow_pct": 100.0 * (scl == SNOW).mean()}


def ndsi(feature: dict) -> tuple:
    """(NDSI, usable mask) on the common grid."""
    g = _read(feature, "green")
    sw = _read(feature, "swir16")
    scl = _read(feature, "scl").astype(int)
    den = g + sw
    n = np.where(den > 0, (g - sw) / np.maximum(den, 1.0), np.nan)
    return n, ~np.isin(scl, list(UNUSABLE))


def composite(features: list[dict]) -> tuple:
    """
    Most recent clear pixel wins. Returns (ndsi, filled).

    Pass features newest first. One post-event scene covers a quarter of this
    AOI, so a single date is not an option.
    """
    out = np.full(GRID, np.nan, np.float32)
    filled = np.zeros(GRID, bool)
    for f in features:
        n, ok = ndsi(f)
        take = ok & ~filled & np.isfinite(n)
        out[take] = n[take]
        filled |= take
        logger.info("%s added %7d px, composite now %.1f%%",
                    f["id"][:24], int(take.sum()), 100 * filled.mean())
    return out, filled


def grow(mask: np.ndarray, seed_rc: tuple) -> list[tuple]:
    """8-connected region grow from a seed. Plain BFS; scipy is not a dep."""
    H, W = mask.shape
    r0, c0 = seed_rc
    if not (0 <= r0 < H and 0 <= c0 < W) or not mask[r0, c0]:
        return []
    seen = np.zeros((H, W), bool)
    seen[r0, c0] = True
    q = deque([(r0, c0)])
    cells = []
    while q:
        r, c = q.popleft()
        cells.append((r, c))
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                rr, cc = r + dr, c + dc
                if 0 <= rr < H and 0 <= cc < W and mask[rr, cc] and not seen[rr, cc]:
                    seen[rr, cc] = True
                    q.append((rr, cc))
    return cells


def to_lonlat(cells: list[tuple], tf) -> tuple:
    from rasterio.transform import xy
    from rasterio.warp import transform as warp_transform
    rr = np.array([c[0] for c in cells])
    cc = np.array([c[1] for c in cells])
    xs, ys = xy(tf, rr, cc)
    lon, lat = warp_transform("EPSG:32645", "EPSG:4326",
                              np.atleast_1d(xs), np.atleast_1d(ys))
    return np.asarray(lon), np.asarray(lat)


def measure(cells, d, pre, tf) -> dict:
    """Extent, centroid and change statistics for a mapped feature."""
    lon, lat = to_lonlat(cells, tf)
    rr = np.array([c[0] for c in cells])
    cc = np.array([c[1] for c in cells])
    dv = d[rr, cc]
    w = np.where(np.isfinite(dv), -dv, 0.0)
    coslat = math.cos(math.radians(float(lat.mean())))
    return {
        "n_px": len(cells), "area_km2": len(cells) * 400 / 1e6,
        "centroid_lat": float(lat.mean()), "centroid_lon": float(lon.mean()),
        # Weighted by the strength of the change, so the centre follows the
        # deepest part of the scar rather than the outline of the mask.
        "weighted_lat": float((lat * w).sum() / w.sum()) if w.sum() else float("nan"),
        "weighted_lon": float((lon * w).sum() / w.sum()) if w.sum() else float("nan"),
        "extent_ns_m": float((lat.max() - lat.min()) * 111320),
        "extent_ew_m": float((lon.max() - lon.min()) * 111320 * coslat),
        "d_ndsi_median": float(np.nanmedian(dv)), "d_ndsi_min": float(np.nanmin(dv)),
        "was_snow_pct": float(100 * np.mean(pre[rr, cc] > 0.4)),
    }


def probe(lat, lon, d, both, pre, tf, radius=5) -> dict:
    """
    Change statistics in a small box about a point.

    The test that matters for a claimed failure location: is there a scar
    there, and if not, is that because it is hidden or because it is absent?
    Reporting the usable-pixel count separates those.
    """
    from rasterio.transform import rowcol
    from rasterio.warp import transform as warp_transform
    x, y = warp_transform("EPSG:4326", "EPSG:32645", [lon], [lat])
    r, c = rowcol(tf, x[0], y[0])
    H, W = d.shape
    r0, r1 = max(0, r - radius), min(H, r + radius + 1)
    c0, c1 = max(0, c - radius), min(W, c + radius + 1)
    sd, sb, sp = d[r0:r1, c0:c1], both[r0:r1, c0:c1], pre[r0:r1, c0:c1]
    n = int(sb.sum())
    if n == 0:
        return {"lat": lat, "lon": lon, "n_usable": 0, "observed": False}
    dv = sd[sb]
    return {"lat": lat, "lon": lon, "n_usable": n, "observed": True,
            "d_ndsi_median": float(np.nanmedian(dv)),
            "d_ndsi_min": float(np.nanmin(dv)),
            "was_snow_pct": float(100 * np.mean(sp[sb] > 0.4)),
            "scar_like_pct": float(100 * np.mean(dv < -0.35))}


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description="Map the detachment scar from Sentinel-2")
    m = ap.add_mutually_exclusive_group(required=True)
    m.add_argument("--survey", action="store_true",
                   help="list scenes and their cloud over the AOI")
    m.add_argument("--map", action="store_true", help="find and measure the scar")
    m.add_argument("--probe", nargs=2, type=float, metavar=("LAT", "LON"),
                   help="change statistics at one location")
    ap.add_argument("--pre", default="S2C_45RUM_20260812_0_L2A",
                    help="pre-event scene id (default: the clearest, 2.6%% AOI cloud)")
    ap.add_argument("--post", nargs="+",
                    default=["S2C_45RUM_20260908_0_L2A", "S2B_45RUM_20260906_0_L2A",
                             "S2B_45RUM_20260903_0_L2A", "S2C_45RUM_20260901_0_L2A",
                             "S2A_45RUM_20260831_0_L2A", "S2C_45RUM_20260829_0_L2A",
                             "S2B_45RUM_20260827_0_L2A"],
                    help="post-event scene ids, newest first")
    ap.add_argument("--seed", nargs=2, type=float, default=[28.27756, 85.52908],
                    metavar=("LAT", "LON"), help="region-grow seed")
    ap.add_argument("--drop", type=float, default=-0.25,
                    help="NDSI fall required to join the feature")
    ap.add_argument("--was-snow", type=float, default=0.35,
                    help="pre-event NDSI required to join the feature")
    ap.add_argument("--los-floor", type=float, default=32.0,
                    help="line-of-sight detection floor at the scar, mm/day, for "
                         "converting to a downslope bound. The default 32.0 is "
                         "the worst of the three intervals covering the seven "
                         "weeks before failure, measured AT the mapped scar by "
                         "local_floor.py - not the 40.4 measured at the "
                         "abandoned failure point 1.09 km away")
    ap.add_argument("--csv", metavar="OUT.csv")
    args = ap.parse_args()

    os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")

    if args.survey:
        feats = search("2026-05-01", "2026-09-30")
        print(f"\n{len(feats)} Sentinel-2 L2A scenes over the source zone")
        print(f"\n{'date':<12}{'scene%':>8}{'AOI cloud%':>12}{'AOI clear%':>12}  id")
        print("-" * 78)
        for f in feats:
            c = aoi_cloud(f)
            flag = "  <-- usable" if c["aoi_clear_pct"] > 60 else ""
            print(f'{c["date"]:<12}{c["scene_cloud"]:>8.1f}{c["aoi_cloud_pct"]:>12.1f}'
                  f'{c["aoi_clear_pct"]:>12.1f}  {c["id"]}{flag}')
        print("\nScene cloud is over the 110 km tile. Only the AOI columns matter.")
        return 0

    F = by_id([args.pre] + args.post)
    missing = [i for i in [args.pre] + args.post if i not in F]
    if missing:
        logger.error("scenes not found: %s", missing)
        return 1
    tf = grid_transform()
    logger.info("pre-event %s", args.pre)
    pre, pre_ok = ndsi(F[args.pre])
    logger.info("post-event composite from %d scenes", len(args.post))
    post, post_ok = composite([F[i] for i in args.post])

    both = pre_ok & post_ok & np.isfinite(pre) & np.isfinite(post)
    d = np.where(both, post - pre, np.nan)
    print(f"\ncomparable in both epochs: {int(both.sum()):,} px "
          f"({100*both.mean():.1f}% of the AOI)")
    base = float(100 * np.mean(d[both] < -0.35))
    print(f"scar-like across the AOI:  {base:.2f}%  "
          f"(the rate an isolated patch has to beat)")

    if args.probe:
        r = probe(args.probe[0], args.probe[1], d, both, pre, tf)
        print(f"\nPROBE {r['lat']:.5f} N {r['lon']:.5f} E")
        if not r["observed"]:
            print("  no usable pixels - cloud in one or both epochs, so this")
            print("  location is untested rather than clear.")
            return 0
        print(f"  usable pixels     {r['n_usable']}")
        print(f"  NDSI change       median {r['d_ndsi_median']:+.3f}, "
              f"min {r['d_ndsi_min']:+.3f}")
        print(f"  was snow or ice   {r['was_snow_pct']:.1f}%")
        print(f"  scar-like         {r['scar_like_pct']:.1f}%  "
              f"against {base:.2f}% across the AOI")
        if r["scar_like_pct"] < 2 * base and r["d_ndsi_median"] > 0:
            print("\n  NOT A SCAR. Observed, and the surface got brighter rather")
            print("  than darker - fresh snow on ground that was never glaciated.")
        return 0

    from rasterio.transform import rowcol
    from rasterio.warp import transform as warp_transform
    x, y = warp_transform("EPSG:4326", "EPSG:32645", [args.seed[1]], [args.seed[0]])
    seed = rowcol(tf, x[0], y[0])
    mask = both & (pre > args.was_snow) & (d < args.drop)
    cells = grow(mask, seed)
    if not cells:
        print("\nNo feature grown from that seed - it is not in the mask.")
        return 1
    s = measure(cells, d, pre, tf)
    print(f"\nMAPPED SCAR")
    print(f"  area              {s['area_km2']:.2f} km2 ({s['n_px']:,} px at 20 m)")
    print(f"  extent            {s['extent_ns_m']:.0f} m N-S x {s['extent_ew_m']:.0f} m E-W")
    print(f"  centroid          {s['centroid_lat']:.5f} N {s['centroid_lon']:.5f} E")
    print(f"  change-weighted   {s['weighted_lat']:.5f} N {s['weighted_lon']:.5f} E")
    print(f"  NDSI change       median {s['d_ndsi_median']:+.3f}, min {s['d_ndsi_min']:+.3f}")
    print(f"  was snow or ice   {s['was_snow_pct']:.1f}%")

    # Terrain sampled ACROSS the feature, not at its centroid. A 1.4 km scar
    # is not one aspect: this one runs 300 to 355 degrees, and the centroid
    # pixel alone gives -0.472 for NISAR ascending where the spread gives
    # -0.47 to -0.61. Quoting the centroid would repeat, at 1 km scale, the
    # single-pixel mistake that put the failure point on the wrong face.
    try:
        from geometry_merge import sensitivities, slope_aspect
        lon, lat = to_lonlat(cells, tf)
        la, lo = s["weighted_lat"], s["weighted_lon"]
        pts = [(la, lo),
               (float(np.percentile(lat, 85)), lo),
               (float(np.percentile(lat, 15)), lo),
               (la, float(np.percentile(lon, 85))),
               (la, float(np.percentile(lon, 15)))]
        print(f"\n  TERRAIN ACROSS THE FEATURE")
        print(f"  {'sample':<10}{'elev':>7}{'slope':>7}{'aspect':>8}")
        print("  " + "-" * 32)
        got = []
        for i, (a, b) in enumerate(pts):
            sl, asp, el = slope_aspect(a, b)
            got.append((sl, asp, el))
            print(f"  {'centroid' if i == 0 else f'edge {i}':<10}"
                  f"{el:>7.0f}{sl:>7.1f}{asp:>8.0f}")
        ang = np.radians([g[1] for g in got])
        masp = float(np.degrees(np.arctan2(np.sin(ang).mean(),
                                           np.cos(ang).mean())) % 360)
        msl = float(np.mean([g[0] for g in got]))
        print(f"\n  elevation {min(g[2] for g in got):.0f}-{max(g[2] for g in got):.0f} m, "
              f"slope {msl:.1f} deg, aspect {masp:.0f} deg (circular mean)")

        print(f"\n  {'track':<16}{'centroid':>10}{'mean aspect':>13}  verdict")
        print("  " + "-" * 52)
        c_rows = {r["track"]: r["sensitivity"]
                  for r in sensitivities(la, lo, got[0][0], got[0][1], 0.3)}
        worst = 0.0
        for r in sorted(sensitivities(la, lo, msl, masp, 0.3),
                        key=lambda r: -abs(r["sensitivity"])):
            t, v = r["track"], r["sensitivity"]
            both_ok = min(abs(v), abs(c_rows.get(t, v)))
            if "NISAR ASC" in t:
                worst = both_ok
            print(f"  {t:<16}{c_rows.get(t, float('nan')):>10.3f}{v:>13.3f}  "
                  f"{'usable' if both_ok >= 0.3 else 'BLIND'}")
        if worst:
            print(f"\n  A {args.los_floor:.1f} mm/day line-of-sight floor at the "
                  f"scar becomes a downslope")
            print(f"  bound of {args.los_floor/abs(worst):.0f} mm/day, taking "
                  f"the weaker of the two aspect samples.")
    except Exception as e:                       # network, mostly
        logger.warning("terrain lookup skipped: %s", e)

    if args.csv:
        import csv as _csv
        p = Path(resolve(args.csv))
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "w", newline="") as fh:
            w = _csv.DictWriter(fh, fieldnames=sorted(s))
            w.writeheader()
            w.writerow(s)
        logger.info("Wrote %s", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
