"""
timeseries.py
-------------
Assemble NISAR GUNW pairs into a cumulative LOS displacement time series.

Each interferogram measures a *difference*:

    phi(i,j)  =  d(t_j) - d(t_i)

Given a network of pairs, invert for d(t) at every acquisition epoch. This is
the small-baseline (SBAS) idea, reduced to what NISAR L2 actually needs -
the products are already unwrapped and geocoded, so all that is left is the
temporal inversion.

Three things this tool refuses to fudge
=======================================
1. It never mixes orbit geometries. Ascending path 98 and descending path 48
   measure different components of the same 3-D motion. They are inverted
   separately and reported separately.

2. It never bridges a disconnected network. If no chain of interferograms links
   two blocks of epochs, their displacements are not comparable and no amount
   of least squares will make them so. Each connected component gets its own
   series with its own zero.

3. It reports redundancy honestly. A component with as many pairs as unknowns
   has an exact solution and therefore no residual and no error estimate. That
   is not a good fit - it is no fit at all.

Usage
-----
    # Network structure straight from the catalogue, before downloading
    python src/timeseries.py --from-catalogue

    # Network structure of what you already have
    python src/timeseries.py --dir data/nisar_l2 --network

    # Region-mean inversion (fast, robust - start here)
    python src/timeseries.py --dir data/nisar_l2 --invert \
        --ref-lat 28.21 --ref-lon 85.47 --csv outputs/ts.csv --plot outputs/ts.png

Requires: numpy.  --dir modes also need h5py (via gunw_reader).
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import sys
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("timeseries")

EVENT_DATE = date(2026, 8, 26)


# ---------------------------------------------------------------------------
# Pair inventory
# ---------------------------------------------------------------------------
class Pair:
    """One interferogram: a measured difference between two epochs."""

    __slots__ = ("path", "direction", "ref", "sec", "value", "n_px", "coherence", "source")

    def __init__(self, path, direction, ref, sec, value=None, n_px=None, coherence=None, source=""):
        self.path = path
        self.direction = direction
        self.ref = ref
        self.sec = sec
        self.value = value          # mm, region mean, filled by --invert
        self.n_px = n_px
        self.coherence = coherence
        self.source = source

    @property
    def span_days(self) -> int:
        return (self.sec - self.ref).days

    def __repr__(self):
        return f"<{self.direction}{self.path} {self.ref}->{self.sec}>"


def pairs_from_granule_names(names: list[str], paths: list, directions: list) -> list[Pair]:
    out = []
    for name, path, direction in zip(names, paths, directions):
        stamps = re.findall(r"_(\d{8})T\d{6}", name)
        if len(stamps) < 4:
            continue
        out.append(
            Pair(
                path=path,
                direction=(direction or "?")[0],
                ref=datetime.strptime(stamps[0], "%Y%m%d").date(),
                sec=datetime.strptime(stamps[2], "%Y%m%d").date(),
                source=name,
            )
        )
    return out


def pairs_from_catalogue() -> list[Pair]:
    """Query ASF for GUNW pairs over the AOI without downloading anything."""
    import nisar_acquisition as acq

    rows = acq.query_by_year("2025-01-01", date.today().isoformat(),
                             platform="NISAR", processingLevel="GUNW")
    return pairs_from_granule_names(
        [r["granuleName"] for r in rows],
        [r.get("path") for r in rows],
        [r.get("flightDirection") for r in rows],
    )


def pairs_from_dir(directory: Path) -> list[Pair]:
    files = sorted(p for p in directory.glob("*.h5") if "GUNW" in p.name.upper())
    if not files:
        logger.error("No *GUNW*.h5 in %s", directory)
        return []
    out = []
    for fp in files:
        stamps = re.findall(r"_(\d{8})T\d{6}", fp.name)
        m = re.search(r"_GUNW_\d+_(\d+)_([AD])_", fp.name)
        if len(stamps) < 4 or not m:
            logger.warning("Cannot parse %s", fp.name)
            continue
        out.append(
            Pair(
                path=int(m.group(1)),
                direction=m.group(2),
                ref=datetime.strptime(stamps[0], "%Y%m%d").date(),
                sec=datetime.strptime(stamps[2], "%Y%m%d").date(),
                source=str(fp),
            )
        )
    return out


# ---------------------------------------------------------------------------
# Network structure
# ---------------------------------------------------------------------------
def connected_components(pairs: list[Pair]) -> list[list[date]]:
    """Group epochs into blocks linked by interferograms (union-find)."""
    parent: dict[date, date] = {}

    def find(a):
        parent.setdefault(a, a)
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for p in pairs:
        union(p.ref, p.sec)

    groups: dict[date, list[date]] = defaultdict(list)
    for epoch in parent:
        groups[find(epoch)].append(epoch)
    return sorted((sorted(v) for v in groups.values()), key=lambda g: g[0])


def describe_network(pairs: list[Pair]) -> dict:
    by_geom: dict[tuple, list[Pair]] = defaultdict(list)
    for p in pairs:
        by_geom[(p.direction, p.path)].append(p)

    print(f"\n{len(pairs)} interferograms across {len(by_geom)} geometries\n")
    summary = {}

    for (direction, path), group in sorted(by_geom.items()):
        comps = connected_components(group)
        epochs = sorted({e for p in group for e in (p.ref, p.sec)})
        label = f"{'ASCENDING' if direction == 'A' else 'DESCENDING'} path {path}"
        print("=" * 68)
        print(f"{label}   {len(group)} pairs, {len(epochs)} epochs, "
              f"{len(comps)} connected component(s)")
        print("=" * 68)

        if len(comps) > 1:
            print("  !! NETWORK IS DISCONNECTED - these blocks cannot be joined.")
            print("     Displacements in one block are NOT comparable to another.")

        comp_info = []
        for k, comp in enumerate(comps, 1):
            inside = [p for p in group if p.ref in comp and p.sec in comp]
            unknowns = len(comp) - 1
            redundancy = len(inside) - unknowns
            print(f"\n  Component {k}: {comp[0]} .. {comp[-1]}  "
                  f"({(comp[-1]-comp[0]).days} d, {len(comp)} epochs, {len(inside)} pairs)")
            print(f"    unknowns {unknowns}, redundancy {redundancy}", end="")
            if redundancy == 0:
                print("  <- exactly determined: no residual, no error estimate")
            elif redundancy < 0:
                print("  <- UNDER-DETERMINED")
            else:
                print()
            for p in sorted(inside, key=lambda x: x.sec):
                flag = "  <- last before event" if p.sec <= EVENT_DATE and \
                    p.sec == max(q.sec for q in inside if q.sec <= EVENT_DATE) else ""
                print(f"      {p.ref} -> {p.sec}  ({p.span_days:>2} d){flag}")
            comp_info.append({"epochs": [str(e) for e in comp], "n_pairs": len(inside),
                              "redundancy": redundancy})
        summary[label] = comp_info
        print()
    return summary


# ---------------------------------------------------------------------------
# Inversion
# ---------------------------------------------------------------------------
def invert_component(pairs: list[Pair], epochs: list[date]) -> dict:
    """
    Least-squares solve for cumulative displacement at each epoch.

    d(epochs[0]) is fixed at 0 - every series is relative to its own first
    acquisition, which is why disconnected components can never be compared.
    """
    epochs = sorted(epochs)
    index = {e: i for i, e in enumerate(epochs)}
    n_unknown = len(epochs) - 1

    usable = [p for p in pairs if p.value is not None]
    if not usable:
        return {"error": "no pairs carry a measured value"}
    if n_unknown == 0:
        return {"error": "single epoch"}

    G = np.zeros((len(usable), n_unknown))
    obs = np.zeros(len(usable))
    for r, p in enumerate(usable):
        i, j = index[p.ref], index[p.sec]
        if i > 0:
            G[r, i - 1] = -1.0
        if j > 0:
            G[r, j - 1] = +1.0
        obs[r] = p.value

    solution, *_ = np.linalg.lstsq(G, obs, rcond=None)
    cumulative = np.concatenate([[0.0], solution])
    residuals = G @ solution - obs

    dof = len(usable) - n_unknown
    if dof > 0:
        sigma = float(np.sqrt(np.sum(residuals ** 2) / dof))
        try:
            cov = sigma ** 2 * np.linalg.inv(G.T @ G)
            errors = np.concatenate([[0.0], np.sqrt(np.diag(cov))])
        except np.linalg.LinAlgError:
            errors = np.full(len(epochs), np.nan)
    else:
        sigma = float("nan")
        errors = np.full(len(epochs), np.nan)

    # linear velocity over the component
    days = np.array([(e - epochs[0]).days for e in epochs], dtype=float)
    velocity = velocity_err = float("nan")
    if len(epochs) >= 3:
        A = np.vstack([days, np.ones_like(days)]).T
        coef, *_ = np.linalg.lstsq(A, cumulative, rcond=None)
        velocity = float(coef[0])
        fit_res = A @ coef - cumulative
        if len(epochs) > 2:
            s = np.sqrt(np.sum(fit_res ** 2) / (len(epochs) - 2))
            try:
                velocity_err = float(s * np.sqrt(np.linalg.inv(A.T @ A)[0, 0]))
            except np.linalg.LinAlgError:
                pass

    return {
        "epochs": epochs,
        "cumulative_mm": cumulative,
        "error_mm": errors,
        "residual_rms_mm": float(np.sqrt(np.mean(residuals ** 2))) if len(residuals) else float("nan"),
        "sigma_mm": sigma,
        "dof": dof,
        "velocity_mm_per_day": velocity,
        "velocity_err_mm_per_day": velocity_err,
        "n_pairs": len(usable),
    }


def measure_pairs(pairs: list[Pair], coh_threshold: float,
                  ref_lat: float | None, ref_lon: float | None,
                  clip_aoi: bool, flip_sign: bool,
                  target_lat: float | None = None, target_lon: float | None = None,
                  target_radius_px: int = 5) -> None:
    """
    Fill Pair.value with a displacement measurement from each GUNW.

    Without a target, the median over the whole valid AOI is used. That is only
    meaningful if the AOI is dominated by the signal you care about - a small
    landslide inside a large stable AOI will be averaged into nothing.

    With --target-lat/--target-lon, the median over a small window at the target
    is used instead. Combined with a reference point on stable ground, that is
    the displacement OF the target RELATIVE TO the reference, which is the only
    thing InSAR can actually measure.
    """
    from gunw_reader import read_gunw

    for p in pairs:
        if not p.source or not Path(p.source).exists():
            continue
        try:
            r = read_gunw(Path(p.source), coherence_threshold=coh_threshold,
                          ref_lat=ref_lat, ref_lon=ref_lon,
                          clip_aoi=clip_aoi, flip_sign=flip_sign)
        except Exception as exc:
            logger.error("%s: %s", Path(p.source).name, exc)
            continue
        valid = r["valid"]
        disp = r["displacement_mm"]

        if target_lat is not None and target_lon is not None and r["xs"] is not None:
            xs, ys = r["xs"], r["ys"]
            tx, ty = target_lon, target_lat
            if r["epsg"] and r["epsg"] != 4326:
                try:
                    from pyproj import Transformer
                    tx, ty = Transformer.from_crs(4326, r["epsg"], always_xy=True).transform(
                        target_lon, target_lat)
                except ImportError:
                    logger.warning("pyproj missing - cannot place target on a projected grid.")
            j = int(np.argmin(np.abs(xs - tx)))
            i = int(np.argmin(np.abs(ys - ty)))
            t = target_radius_px
            sel = np.zeros_like(valid)
            sel[max(0, i - t): i + t + 1, max(0, j - t): j + t + 1] = True
            valid = valid & sel

        n = int(valid.sum())
        if n == 0:
            logger.warning("%s: no valid pixels in the measurement window, excluded",
                           Path(p.source).name)
            continue
        p.value = float(np.median(disp[valid]))
        p.n_px = n
        p.coherence = (float(r["coherence"][valid].mean())
                       if r["coherence"] is not None else None)
        logger.info("%s -> %s : %+8.2f mm  (%d px, coh %.2f)",
                    p.ref, p.sec, p.value, n,
                    p.coherence if p.coherence is not None else float("nan"))


def report_series(label: str, comp_index: int, result: dict) -> list[dict]:
    if "error" in result:
        print(f"\n{label} component {comp_index}: {result['error']}")
        return []

    print(f"\n{label}  -  component {comp_index}")
    print(f"  {result['n_pairs']} pairs, {len(result['epochs'])} epochs, dof {result['dof']}")
    if result["dof"] > 0:
        print(f"  residual RMS {result['residual_rms_mm']:.2f} mm  (sigma {result['sigma_mm']:.2f} mm)")
    else:
        print("  exactly determined - residual and error are undefined")

    print(f"\n  {'EPOCH':<12}{'CUMULATIVE mm':>15}{'+/- mm':>10}{'DAYS':>7}")
    print("  " + "-" * 44)
    t0 = result["epochs"][0]
    rows = []
    for e, d, err in zip(result["epochs"], result["cumulative_mm"], result["error_mm"]):
        errs = "   n/a" if not np.isfinite(err) else f"{err:>10.2f}"
        print(f"  {str(e):<12}{d:>15.2f}{errs}{(e-t0).days:>7}")
        rows.append({"geometry": label, "component": comp_index, "epoch": str(e),
                     "days_from_start": (e - t0).days,
                     "cumulative_mm": round(float(d), 3),
                     "error_mm": None if not np.isfinite(err) else round(float(err), 3)})

    v, ve = result["velocity_mm_per_day"], result["velocity_err_mm_per_day"]
    if np.isfinite(v):
        span = (result["epochs"][-1] - t0).days
        pm = f" +/- {ve:.3f}" if np.isfinite(ve) else ""
        print(f"\n  linear velocity {v:+.3f}{pm} mm/day   "
              f"({v*365:+.0f} mm/yr, over {span} d)")
        if np.isfinite(ve) and ve > 0 and abs(v) < 2 * ve:
            print("  NOT significant at 2 sigma - consistent with no motion.")
    return rows


# ---------------------------------------------------------------------------
def plot_series(all_rows: list[dict], out: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.error("matplotlib required for --plot:  pip install matplotlib")
        return

    groups = defaultdict(list)
    for r in all_rows:
        groups[(r["geometry"], r["component"])].append(r)

    fig, ax = plt.subplots(figsize=(10, 6))
    for (geom, comp), rows in sorted(groups.items()):
        xs = [datetime.strptime(r["epoch"], "%Y-%m-%d") for r in rows]
        ys = [r["cumulative_mm"] for r in rows]
        es = [r["error_mm"] or 0 for r in rows]
        ax.errorbar(xs, ys, yerr=es, marker="o", capsize=3,
                    label=f"{geom} (block {comp})")
    ax.axvline(datetime(2026, 8, 26), color="crimson", ls="--", lw=1.2)
    ax.text(datetime(2026, 8, 26), ax.get_ylim()[1], " 26 Aug event",
            color="crimson", va="top", fontsize=9)
    ax.set_ylabel("Cumulative LOS displacement (mm)\n+ve away from satellite")
    ax.set_title("NISAR L-band cumulative displacement\n"
                 "each block has its own zero - blocks are not comparable",
                 fontsize=11)
    ax.grid(alpha=.3)
    ax.legend(fontsize=9)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    logger.info("Wrote %s", out)


def main() -> int:
    ap = argparse.ArgumentParser(description="GUNW pairs -> cumulative LOS time series")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--dir", metavar="DIR", help="directory of GUNW .h5 files")
    src.add_argument("--from-catalogue", action="store_true",
                     help="analyse network from ASF without downloading")
    ap.add_argument("--network", action="store_true", help="network structure only")
    ap.add_argument("--invert", action="store_true", help="run the inversion (needs --dir)")
    ap.add_argument("--coh-threshold", type=float, default=0.3)
    ap.add_argument("--ref-lat", type=float)
    ap.add_argument("--ref-lon", type=float)
    ap.add_argument("--target-lat", type=float, help="measure here instead of the AOI median")
    ap.add_argument("--target-lon", type=float)
    ap.add_argument("--target-radius", type=int, default=5, help="target window half-width in px")
    ap.add_argument("--no-clip", action="store_true")
    ap.add_argument("--flip-sign", action="store_true")
    ap.add_argument("--csv", metavar="OUT.csv")
    ap.add_argument("--plot", metavar="OUT.png")
    ap.add_argument("--json", metavar="OUT.json", help="dump the network summary")
    args = ap.parse_args()

    pairs = pairs_from_catalogue() if args.from_catalogue else pairs_from_dir(Path(args.dir))
    if not pairs:
        logger.error("No pairs found.")
        return 1

    summary = describe_network(pairs)
    if args.json:
        Path(args.json).write_text(json.dumps(summary, indent=1))
        logger.info("Wrote %s", args.json)

    if args.network or not args.invert:
        if args.from_catalogue:
            print("Nothing downloaded yet - run nisar_acquisition.py --download GUNW "
                  "then re-run with --dir to invert.")
        return 0

    if args.from_catalogue:
        logger.error("--invert needs local files. Use --dir.")
        return 1

    if args.ref_lat is None or args.ref_lon is None:
        logger.warning(
            "No --ref-lat/--ref-lon. Each interferogram keeps its own arbitrary "
            "offset, so the assembled series will be meaningless. Strongly advised."
        )

    logger.info("Measuring displacement for %d pairs...", len(pairs))
    if args.target_lat is None:
        logger.warning(
            "No --target-lat/--target-lon: using the median over the whole AOI. "
            "A localised signal will be averaged away. Give a target for hazard work."
        )
    measure_pairs(pairs, args.coh_threshold, args.ref_lat, args.ref_lon,
                  not args.no_clip, args.flip_sign,
                  args.target_lat, args.target_lon, args.target_radius)

    by_geom = defaultdict(list)
    for p in pairs:
        by_geom[(p.direction, p.path)].append(p)

    all_rows: list[dict] = []
    for (direction, path), group in sorted(by_geom.items()):
        label = f"{'ASC' if direction == 'A' else 'DESC'} path {path}"
        for k, comp in enumerate(connected_components(group), 1):
            inside = [p for p in group if p.ref in comp and p.sec in comp]
            all_rows += report_series(label, k, invert_component(inside, comp))

    if args.csv and all_rows:
        keys = ["geometry", "component", "epoch", "days_from_start", "cumulative_mm", "error_mm"]
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(all_rows)
        logger.info("Wrote %s (%d rows)", args.csv, len(all_rows))

    if args.plot and all_rows:
        plot_series(all_rows, Path(args.plot))

    return 0


if __name__ == "__main__":
    sys.exit(main())
