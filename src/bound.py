"""
bound.py
--------
The headline, computed in one place from the committed measurements.

Why this exists
===============
The headline used to be assembled by hand in four places - README, demo,
site, video script - and it drifted in all of them. It also compared unlike
things: a NISAR floor expressed as DOWNSLOPE motion against a Sentinel-1
precursor that the reports never say is downslope, and a floor for ONE PIXEL
against a rate that is an average over the slope. Both choices inflated the
gap, and together they turned "about 30-100x" into "180-360x".

So this module states each comparison like for like and computes it:

    line of sight against line of sight   the precursor as reported, against
                                          the NISAR ascending LOS floor
    per pixel                             the pessimistic end: one pixel's
                                          3-sigma floor at the candidate
    window median                         the optimistic end: 3 MADs of the
                                          1 km block medians across the AOI,
                                          which includes the correlated error
                                          a single window cannot see

and it does so at EVERY detachment-like candidate, because a floor measured
at one candidate is not the floor at the others. The candidate that matches
the published detachment elevation is flagged, not silently chosen.

Inputs are external facts stated as constants with their source; everything
else is read from outputs/.

    python src/bound.py            # table, and outputs/bound_source.csv
"""

from __future__ import annotations

import csv
import math
import sys
from datetime import date, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from geometry_merge import sensitivities  # noqa: E402
from local_floor import point_in_aoi  # noqa: E402

OUT = ROOT / "outputs"
EVENT = date(2026, 8, 26)

# External inputs, with where they come from. Not measurements of this project.
#
# Shirzaei (Virginia Tech), Sentinel-1 InSAR, 8 Jan - 18 Aug 2026, as reported
# in the press in September 2026: "roughly 10 millimetres per month", with
# acceleration over the final weeks. The reports do not say whether it is line
# of sight or downslope, or which track; it is taken here as line of sight,
# the form an InSAR rate is normally reported in, and compared as such.
PRECURSOR_MM_DAY = 10.0 / 30.0
# Published accounts: the detachment was on the north face at about
# 5,200-5,400 m (Wikipedia, 2026 Nepal-Tibet floods, citing USGS and GFZ).
PUBLISHED_ELEV_M = (5200.0, 5400.0)
MIN_ABS_SENSITIVITY = 0.25  # below this the LOS-to-downslope conversion is ill-conditioned

COVER_DAYS = 60            # the seven weeks before failure, plus the 24-day pair
ASC_REFS = {"20251128", "20251210", "20251222", "20260103",
            "20260702", "20260714", "20260726", "20260819"}


def _d(s: str) -> date:
    return datetime.strptime(s, "%Y%m%d").date()


def pair_of(name: str) -> tuple[date, date]:
    """(reference, secondary) from 'GOFF_YYYYMMDD_YYYYMMDD_...' or 'GUNW_...'."""
    parts = Path(name).name.split("_")
    return _d(parts[1]), _d(parts[2])


def covers_pre_event(name: str) -> bool:
    """Ends before the collapse, and starts inside the pre-failure window."""
    ref, sec = pair_of(name)
    return sec < EVENT and (EVENT - ref).days <= COVER_DAYS


def load(name: str) -> list[dict]:
    with open(OUT / name, newline="") as fh:
        return list(csv.DictReader(fh))


def worst(rows: list[dict], key: str, *, positive_only: bool = False) -> float:
    """Largest eligible finite value; NaN/Inf must never become a reported bound.

    A detection floor is a strictly positive noise estimate. Zero is not a
    scientifically valid "perfect instrument" floor; negative floors are
    invalid too. Set positive_only for any quantity used as a noise floor.
    """
    vals = []
    for row in rows:
        try:
            value = float(row.get(key, "nan"))
        except (TypeError, ValueError, OverflowError):
            continue
        if math.isfinite(value) and (not positive_only or value > 0):
            vals.append(value)
    return max(vals) if vals else float("nan")


def eligible_floor_row(row: dict, floor_key: str) -> bool:
    """Require sample support and a finite positive floor; ignore stale usable flags."""
    if str(row.get("usable", "")).strip().lower() not in {"true", "1", "yes"}:
        return False
    try:
        n = int(float(row.get("window_px", "nan")))
        floor = float(row.get(floor_key, "nan"))
    except (TypeError, ValueError, OverflowError):
        return False
    return n >= 6 and math.isfinite(floor) and floor > 0


def eligible_source_candidates(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """Only candidate centres inside the exact source polygon may enter this bound."""
    eligible, excluded = [], []
    for row in rows:
        try:
            lat, lon = float(row["lat"]), float(row["lon"])
        except (KeyError, TypeError, ValueError, OverflowError):
            excluded.append(row)
            continue
        if point_in_aoi(lon, lat, "source"):
            eligible.append(row)
        else:
            excluded.append(row)
    return eligible, excluded


def assess(cand: dict, goff: list[dict], phase: list[dict], cref: list[dict]) -> dict:
    """Every number the headline needs, for one candidate."""
    cid = str(cand["id"])
    lat, lon = float(cand["lat"]), float(cand["lon"])
    if not point_in_aoi(lon, lat, "source"):
        raise ValueError(f"Candidate {cid} is outside the exact source AOI; refusing to compute its source-zone bound.")
    slope, aspect, elev = (float(cand[k]) for k in ("slope_deg", "aspect_deg", "elev_m"))
    sens_rows = [r for r in sensitivities(lat, lon, slope, aspect, 0.3)
                 if r["track"] == "NISAR ASC 98"]
    if len(sens_rows) != 1:
        raise ValueError(f"Expected exactly one NISAR ASC 98 sensitivity for candidate {cid}; found {len(sens_rows)}.")
    sens = float(sens_rows[0]["sensitivity"])
    if not math.isfinite(sens):
        sens = float("nan")

    g = [
        r for r in goff
        if str(r.get("target_id")) == cid
        and eligible_floor_row(r, "local_floor_mm_day")
        and str(r.get("file", ""))[5:13] in ASC_REFS
        and covers_pre_event(str(r.get("file", "")))
    ]
    p = [
        r for r in phase
        if str(r.get("target_id")) == cid
        and eligible_floor_row(r, "local_floor_mm_day")
    ]
    c = []
    for row in cref:
        if str(row.get("target_id")) != cid:
            continue
        try:
            window_px = int(row.get("window_px", 0))
            common_velocity = float(row.get("target_mm_day_common_ref", "nan"))
        except (TypeError, ValueError, OverflowError):
            continue
        if window_px >= 6 and math.isfinite(common_velocity):
            c.append(row)

    px = worst(g, "local_floor_mm_day", positive_only=True)
    blk = worst(g, "block_median_floor_mm_day", positive_only=True)
    common_velocities = []
    for row in c:
        try:
            value = abs(float(row["target_mm_day_common_ref"]))
        except (KeyError, TypeError, ValueError, OverflowError):
            continue
        if math.isfinite(value):
            common_velocities.append(value)
    v_common = max(common_velocities, default=float("nan"))
    downslope_resolved = math.isfinite(sens) and abs(sens) >= MIN_ABS_SENSITIVITY
    return {
        "id": cid, "lat": lat, "lon": lon, "elev_m": elev, "aspect_deg": round(aspect, 1),
        "matches_published_elevation": PUBLISHED_ELEV_M[0] <= elev <= PUBLISHED_ELEV_M[1],
        "nisar_asc_sensitivity": round(sens, 3),
        "goff_intervals": len(g),
        "los_floor_pixel_mm_day": round(px, 2),
        "los_floor_block_mm_day": round(blk, 2),
        "downslope_geometry_resolved": downslope_resolved,
        "downslope_bound_mm_day": (round(px / abs(sens), 1)
                                    if g and downslope_resolved and math.isfinite(px)
                                    else float("nan")),
        "ratio_los_pixel": round(px / PRECURSOR_MM_DAY, 1) if g else float("nan"),
        "ratio_los_block": round(blk / PRECURSOR_MM_DAY, 1) if g else float("nan"),
        "max_speed_common_datum_mm_day": round(v_common, 2),
        "phase_pairs": len(p),
        "phase_floor_pixel_mm_day": round(worst(p, "local_floor_mm_day", positive_only=True), 2),
        "phase_floor_block_mm_day": round(worst(p, "block_median_floor_mm_day", positive_only=True), 2),
        "precursor_rate_assumption": "10 mm/month interpreted as LOS; source reports do not specify direction",
    }


def main() -> int:
    all_candidates = [r for r in load("scar_candidates_source.csv")
                       if r.get("reading") == "DETACHMENT-like"]
    cands, excluded = eligible_source_candidates(all_candidates)
    if excluded:
        labels = ", ".join(str(row.get("id", "(unknown id)")) for row in excluded)
        print(f"\nEXCLUDED {len(excluded)} candidate(s) outside/malformed for the exact source AOI: {labels}")
    if not cands:
        print("ERROR: no detachment-like candidate has a valid centre inside the exact source AOI; no bound will be computed.")
        return 1
    goff = load("local_floor_candidates.csv")
    phase = load("local_floor_phase_candidates.csv")
    cref = load("common_ref_summer.csv")
    rows = [assess(c, goff, phase, cref) for c in cands]

    print(f"\nPRECURSOR  {PRECURSOR_MM_DAY:.2f} mm/day (Sentinel-1, as reported; "
          f"taken as line of sight)\n")
    print(f"  {'#':>2}{'elev':>6}{'pub?':>6}{'sens':>7}{'LOS px':>8}{'LOS blk':>9}"
          f"{'x px':>7}{'x blk':>7}{'downsl':>8}{'|v| cd':>8}{'phase':>8}")
    print("  " + "-" * 78)
    for r in rows:
        ph = (f"{r['phase_floor_pixel_mm_day']:.1f}" if r["phase_pairs"] else "none")
        if r["goff_intervals"]:
            print(f"  {r['id']:>2}{r['elev_m']:>6.0f}{'yes' if r['matches_published_elevation'] else 'no':>6}"
                  f"{r['nisar_asc_sensitivity']:>7.3f}{r['los_floor_pixel_mm_day']:>8.1f}"
                  f"{r['los_floor_block_mm_day']:>9.1f}{r['ratio_los_pixel']:>7.0f}"
                  f"{r['ratio_los_block']:>7.0f}{r['downslope_bound_mm_day']:>8.0f}"
                  f"{r['max_speed_common_datum_mm_day']:>8.1f}{ph:>8}")
        else:
            print(f"  {r['id']:>2}{r['elev_m']:>6.0f}{'yes' if r['matches_published_elevation'] else 'no':>6}"
                  f"{r['nisar_asc_sensitivity']:>7.3f}   no valid offsets - UNOBSERVED{ph:>17}")
    print("\n  LOS px   worst per-pixel 3-sigma floor over the intervals covering the")
    print("           seven weeks before failure, at the candidate (pessimistic)")
    print("  LOS blk  worst floor on a 1 km window median, from AOI block scatter")
    print("           (optimistic)")
    print("  x        each divided by the precursor - LOS against LOS")
    print("  |v| cd   fastest covering interval at the candidate on a common datum")
    print("  phase    worst winter per-pixel phase floor; 'none' = no valid phase")
    print("           pixel in any GUNW pair at that candidate")

    with open(OUT / "bound_source.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\n  wrote {OUT / 'bound_source.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
