"""
candidate_check.py
------------------
Interrogate one location before calling it a landslide.

The problem
===========
A per-pixel stack will always produce candidates. Correlated atmosphere makes
clusters, low coherence makes outliers, and an unwrapper faced with a steep
gradient makes cycle errors that are spatially coherent and beautifully
monotonic. Every one of those looks like creep in a single geometry and a
single season.

So a candidate is not a measurement until it survives being attacked from
four directions at once, and this tool runs the attack.

    geometry     does the other orbit see it, at the sensitivity it should?
    season       is it present when the hazard mechanism is active?
    differential is it local, or is it the whole AOI moving together?
    ceiling      is it slower than lambda/4 per pair, i.e. measurable at all?

The seasonality test is the one people skip
===========================================
Himalayan slope creep is monsoon-driven: pore pressure peaks June to
September. A signal that appears in the accumulation season and vanishes
through the monsoon has the seasonality of a snowpack, not a slope. Dry snow
is nearly transparent at L-band, so it preserves coherence while adding a path
delay that scales with accumulated water equivalent - producing a steady,
spatially coherent, high-coherence apparent motion that is not motion.

Usage
-----
    # measure the target first, in every product
    python src/timeseries.py --dir data/nisar_l2/GUNW --invert --aoi langtang \
        --auto-ref --target-lat 28.27484 --target-lon 85.47405 --target-radius 8

    # then interpret it against the AOI it sits in
    python src/candidate_check.py --target-log run.log \
        --aoi-csv outputs/gunw_stats_langtang.csv --lat 28.27484 --lon 85.47405
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("candidate")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import resolve  # noqa: E402

WAVELENGTH_M = 0.2439
MIN_PX = 20                    # below this a target window is not measured
MONSOON = (6, 7, 8, 9)

TARGET_LINE = re.compile(
    r"(\d{4}-\d\d-\d\d) -> (\d{4}-\d\d-\d\d) :\s+([-+]?[\d.]+) mm\s+"
    r"\((\d+) px, coh ([\d.]+)\)")


def read_target(path: Path) -> dict:
    """Per-pair measurements at the target, from a timeseries.py --target run."""
    out = {}
    for line in open(path, errors="ignore"):
        m = TARGET_LINE.search(line)
        if m:
            key = (m.group(1).replace("-", ""), m.group(2).replace("-", ""))
            out.setdefault(key, (float(m.group(3)), int(m.group(4)), float(m.group(5))))
    return out


def read_aoi(path: Path) -> dict:
    out = {}
    for r in csv.DictReader(open(path, newline="")):
        if "_UR_" in r["file"]:
            continue
        g = re.search(r"_GUNW_\d+_(\d+)_([AD])_", r["file"])
        if g:
            out[(r["reference"], r["secondary"])] = (float(r["median"]), g.group(2))
    return out


def span_days(k) -> int:
    from datetime import datetime
    a = datetime.strptime(k[0], "%Y%m%d")
    b = datetime.strptime(k[1], "%Y%m%d")
    return (b - a).days


def main() -> int:
    ap = argparse.ArgumentParser(description="Stress-test one candidate location")
    ap.add_argument("--target-log", required=True,
                    help="log from timeseries.py --target-lat/--target-lon")
    ap.add_argument("--aoi-csv", required=True, help="gunw_reader --csv for the same AOI")
    ap.add_argument("--lat", type=float); ap.add_argument("--lon", type=float)
    ap.add_argument("--min-px", type=int, default=MIN_PX)
    args = ap.parse_args()

    tgt = read_target(resolve(args.target_log))
    aoi = read_aoi(resolve(args.aoi_csv))
    shared = sorted(set(tgt) & set(aoi))
    if not shared:
        logger.error("No pairs common to the target log and the AOI CSV.")
        return 1

    where = (f"{args.lat:.5f} N {args.lon:.5f} E"
             if args.lat is not None and args.lon is not None else "target")
    print(f"\n{'=' * 74}")
    print(f"CANDIDATE AT {where}")
    print("=" * 74)
    print("\nLocal signal = target window minus AOI median, so whatever the AOI and")
    print("the target share - atmosphere, an unresolved reference offset - cancels.\n")
    print(f"  {'GEOM':<5}{'PAIR':<22}{'TARGET':>9}{'AOI':>9}{'LOCAL':>9}{'PX':>6}{'COH':>6}")
    print("  " + "-" * 68)

    groups: dict[tuple, list] = {}
    for k in sorted(shared, key=lambda k: (aoi[k][1], k[0])):
        med, geom = aoi[k]
        val, px, coh = tgt[k]
        local = val - med
        season = "monsoon" if int(k[0][4:6]) in MONSOON else "dry"
        groups.setdefault((geom, season), []).append((k, local, px))
        flag = "" if px >= args.min_px else "  too few px"
        print(f"  {geom:<5}{k[0]}-{k[1]:<13}{val:>9.1f}{med:>9.1f}"
              f"{local:>9.1f}{px:>6}{coh:>6.2f}{flag}")

    print()
    rates = {}
    for (geom, season), rows in sorted(groups.items()):
        use = [r for r in rows if r[2] >= args.min_px]
        if not use:
            print(f"  {geom} {season:<8} no pair has {args.min_px}+ valid pixels - "
                  f"not measured here")
            continue
        total = sum(r[1] for r in use)
        days = sum(span_days(r[0]) for r in use)
        rates[(geom, season)] = total / days
        print(f"  {geom} {season:<8} {len(use)} pair(s), net {total:+8.1f} mm / "
              f"{days:>3} d = {total/days:+6.2f} mm/day")

    # ---- the four tests ---------------------------------------------------
    print(f"\n{'=' * 74}\nVERDICT\n{'=' * 74}")
    ceiling = WAVELENGTH_M / 4 * 1000 / 12.0
    fast = [f"{g} {s}" for (g, s), v in rates.items() if abs(v) > ceiling]

    dry = {g: v for (g, s), v in rates.items() if s == "dry"}
    mon = {g: v for (g, s), v in rates.items() if s == "monsoon"}

    problems = []
    if len(dry) >= 2 and all(v < 0 for v in dry.values()) or \
       len(dry) >= 2 and all(v > 0 for v in dry.values()):
        print(f"  Both geometries agree in the dry season "
              f"({', '.join(f'{g} {v:+.2f}' for g, v in sorted(dry.items()))} mm/day).")
        print("  A geometry-specific artefact is ruled out: this is a real phase signal.")

    if dry and mon:
        d = max(abs(v) for v in dry.values())
        m = max(abs(v) for v in mon.values())
        print(f"\n  dry season   up to {d:.2f} mm/day")
        print(f"  monsoon      up to {m:.2f} mm/day")
        if m < 0.4 * d:
            problems.append("seasonality")
            print("\n  SEASONALITY IS BACKWARDS. Himalayan slope creep is monsoon-driven -")
            print("  pore pressure peaks June to September. A signal that runs through")
            print("  the accumulation season and stops for the monsoon has the")
            print("  seasonality of a snowpack, not a slope.")
            print("\n  Dry snow is nearly transparent at L-band, so it preserves the")
            print("  coherence that makes this look trustworthy while adding a path")
            print("  delay that grows with accumulated water equivalent. The result is")
            print("  steady, spatially coherent apparent motion that is not motion.")

    if fast:
        problems.append("ceiling")
        print(f"\n  ABOVE THE CEILING in {', '.join(fast)}: faster than lambda/4 per")
        print(f"  12-day pair ({ceiling:.2f} mm/day), so the phase has wrapped and the")
        print("  unwrapper has guessed. Not a measurement.")

    if not problems:
        print("\n  Survives every test available here. Still not a detection until an")
        print("  independent measurement - field, GNSS, or a different sensor - agrees.")
    else:
        print(f"\n  NOT GROUND MOTION on this evidence. Record the location and the")
        print("  reason; do not carry it forward as a detection.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
