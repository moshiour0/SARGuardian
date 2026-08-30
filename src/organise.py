"""
organise.py
-----------
Sort a flat pile of NISAR L2 downloads into a structure a team can share, and
write a manifest so anyone can see what is present without opening 20 GB.

NISAR granule names encode everything needed:

    NISAR_L2_PR_GUNW_006_048_D_074_007_2000_SH_20251125T125813_..._20251207T...
              |   |    |   |  |            |   reference date   secondary date
              |   |    |   |  |            polarisation
              |   |    |   |  frame
              |   |    |   orbit direction
              |   |    track/path
              |   product type
              processing type: PR = routine, UR = urgent response

Layout produced:

    data/nisar_l2/
        GUNW/2025-11_winter/...
        GOFF/2026-07_summer/...
        GOFF/co_event/...          anything whose secondary date is after the event
        MANIFEST.csv

Moves by default (these files are 1-2 GB; copying wastes disk). Dry-run unless
--apply is given.

Usage
-----
    python src/organise.py --src src/workspace_langtang/nisar_l2
    python src/organise.py --src src/workspace_langtang/nisar_l2 --apply
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import shutil
import sys
from datetime import date, datetime
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("organise")

EVENT_DATE = date(2026, 8, 26)
DEST_DEFAULT = Path("data/nisar_l2")

NAME = re.compile(
    r"NISAR_L2_(?P<proc>[A-Z]{2})_(?P<product>[A-Z]{4})_"
    r"(?P<cycle>\d+)_(?P<path>\d+)_(?P<direction>[AD])_(?P<frame>\d+)_"
    r"(?P<sub>\d+)_(?P<bw>\d+)_(?P<pol>[A-Z]{2})_"
)


def parse(path: Path) -> dict | None:
    m = NAME.match(path.name)
    if not m:
        return None
    stamps = re.findall(r"_(\d{8})T\d{6}", path.name)
    if len(stamps) < 4:
        return None
    ref = datetime.strptime(stamps[0], "%Y%m%d").date()
    sec = datetime.strptime(stamps[2], "%Y%m%d").date()
    d = m.groupdict()
    return {
        "file": path.name,
        "product": d["product"],
        "processing": "urgent" if d["proc"] == "UR" else "routine",
        "path": int(d["path"]),
        "direction": d["direction"],
        "polarisation": d["pol"],
        "reference": ref,
        "secondary": sec,
        "span_days": (sec - ref).days,
        "days_before_event": (EVENT_DATE - sec).days,
        "size_mb": round(path.stat().st_size / 1024 / 1024),
        "source": path,
    }


def bucket(rec: dict) -> str:
    """Which subfolder this pair belongs in."""
    if rec["secondary"] > EVENT_DATE:
        return "co_event"
    month = rec["reference"].month
    season = "winter" if month in (11, 12, 1, 2, 3) else "summer"
    return f"{rec['reference']:%Y-%m}_{season}"


def main() -> int:
    ap = argparse.ArgumentParser(description="Organise NISAR L2 downloads")
    ap.add_argument("--src", required=True, help="folder holding the flat .h5 pile")
    ap.add_argument("--dest", default=str(DEST_DEFAULT))
    ap.add_argument("--apply", action="store_true", help="actually move (default: dry run)")
    ap.add_argument("--copy", action="store_true", help="copy instead of move")
    args = ap.parse_args()

    src = Path(args.src)
    dest = Path(args.dest)
    if not src.is_dir():
        logger.error("No such folder: %s", src)
        return 1

    files = sorted(src.glob("*.h5"))
    if not files:
        logger.error("No .h5 files in %s", src)
        return 1

    recs, skipped = [], []
    for fp in files:
        r = parse(fp)
        (recs if r else skipped).append(r or fp.name)

    for s in skipped:
        logger.warning("Unparsed, left in place: %s", s)

    recs.sort(key=lambda r: (r["product"], r["reference"], r["path"]))
    total_gb = sum(r["size_mb"] for r in recs) / 1024

    print(f"\n{len(recs)} products, {total_gb:.1f} GB\n")
    print(f"{'PRODUCT':<8}{'PROC':<9}{'DIR':>4}{'PATH':>6}{'POL':>5}  "
          f"{'REFERENCE':<12}{'SECONDARY':<12}{'SPAN':>6}{'PRE-EV':>8}{'MB':>7}  -> BUCKET")
    print("-" * 104)
    for r in recs:
        b = bucket(r)
        flag = "  <== SPANS THE COLLAPSE" if b == "co_event" else ""
        print(f"{r['product']:<8}{r['processing']:<9}{r['direction']:>4}{r['path']:>6}"
              f"{r['polarisation']:>5}  {str(r['reference']):<12}{str(r['secondary']):<12}"
              f"{r['span_days']:>4} d{r['days_before_event']:>6} d{r['size_mb']:>7}"
              f"  -> {r['product']}/{b}{flag}")

    # summary by product and bucket
    from collections import Counter
    print()
    for prod in sorted({r["product"] for r in recs}):
        sub = [r for r in recs if r["product"] == prod]
        gb = sum(r["size_mb"] for r in sub) / 1024
        buckets = Counter(bucket(r) for r in sub)
        print(f"  {prod}: {len(sub)} files, {gb:.1f} GB  "
              + ", ".join(f"{k}={v}" for k, v in sorted(buckets.items())))

    if not args.apply:
        print(f"\nDRY RUN. Nothing moved. Re-run with --apply to {'copy' if args.copy else 'move'}.")
        return 0

    dest.mkdir(parents=True, exist_ok=True)
    moved = 0
    for r in recs:
        target_dir = dest / r["product"] / bucket(r)
        target_dir.mkdir(parents=True, exist_ok=True)
        target = target_dir / r["file"]
        if target.exists():
            logger.info("Already there: %s", r["file"][:60])
            continue
        if args.copy:
            shutil.copy2(r["source"], target)
        else:
            shutil.move(str(r["source"]), str(target))
        moved += 1
        logger.info("%s -> %s", r["file"][:52], target_dir)

    manifest = dest / "MANIFEST.csv"
    keys = ["product", "processing", "direction", "path", "polarisation",
            "reference", "secondary", "span_days", "days_before_event",
            "size_mb", "bucket", "file"]
    with open(manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in recs:
            w.writerow({**r, "bucket": bucket(r),
                        "reference": str(r["reference"]), "secondary": str(r["secondary"])})
    logger.info("Moved %d files. Manifest: %s", moved, manifest)

    leftovers = list(src.glob("*"))
    if not leftovers:
        logger.info("Source folder is empty - safe to delete: %s", src)
    return 0


if __name__ == "__main__":
    sys.exit(main())
