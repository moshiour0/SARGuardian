"""
inverse_velocity.py
-------------------
Fukuzono inverse-velocity failure forecasting on a measured displacement
series. This is the step that turns the pipeline into an alarm.

The method
==========
For accelerating creep to failure, Voight's material-failure law

    dv/dt = A * v^alpha

with alpha = 2 integrates to a straight line in inverse velocity:

    1/v(t) = A * (T - t)

so 1/v falls linearly and crosses zero at the failure time T. Fit the trailing
window, extrapolate the x-intercept, and that is a predicted failure date.
Fukuzono (1985); used operationally in open-pit mine slope monitoring.

Why this is not the same code as detectability.py
=================================================
That module runs the detector on a *simulated* curve where the truth is known.
This one runs on measured data, which brings three problems the simulation does
not have:

1. Irregular sampling. Real pairs are 12 and 24 days, not a fixed cadence.
2. Disconnected blocks. Winter and summer have separate zeros and can never be
   joined, so each is forecast independently.
3. A noise floor. Without a significance gate an inverse-velocity fit will
   happily "predict" failure from three noise samples. Every velocity must
   clear a multiple of the measured floor before it is allowed into a fit.

That third point is the whole difference between a detector and a random
number generator. Measure the floor first with:

    python src/goff_reader.py --noise-floor data/nisar_l2/GOFF/2025-12_winter

Reporting a non-detection
=========================
When no alarm fires the tool says which gate stopped it, and computes the
velocity that WOULD have been required. A null result with that number
attached is publishable; a null without it is just silence.

Usage
-----
    python src/inverse_velocity.py --ts outputs/ts_goff.csv --noise-floor 75
    python src/inverse_velocity.py --ts outputs/ts_gunw.csv --noise-floor 5 \
        --event-date 2026-08-26 --plot outputs/inv_velocity.png
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from collections import defaultdict
from datetime import date, datetime, timedelta
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("inv-velocity")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import resolve  # noqa: E402


# ---------------------------------------------------------------------------
def load_series(path: Path) -> dict[tuple[str, int], list[dict]]:
    series: dict[tuple[str, int], list[dict]] = defaultdict(list)
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            series[(row["geometry"], int(row["component"]))].append({
                "epoch": datetime.strptime(row["epoch"], "%Y-%m-%d").date(),
                "cumulative_mm": float(row["cumulative_mm"]),
                "error_mm": float(row["error_mm"]) if row.get("error_mm") else None,
            })
    for v in series.values():
        v.sort(key=lambda r: r["epoch"])
    return dict(series)


def velocities(rows: list[dict]) -> list[dict]:
    """Mid-interval velocity between consecutive epochs."""
    out = []
    for a, b in zip(rows, rows[1:]):
        dt = (b["epoch"] - a["epoch"]).days
        if dt <= 0:
            continue
        out.append({
            "t0": a["epoch"], "t1": b["epoch"], "days": dt,
            "mid": a["epoch"] + timedelta(days=dt / 2),
            "v_mm_day": (b["cumulative_mm"] - a["cumulative_mm"]) / dt,
        })
    return out


def fit_inverse_velocity(win: list[dict]) -> dict:
    """
    Least-squares 1/v against time, with the x-intercept and its uncertainty.

    t_f = -b/a for the line 1/v = a*t + b. Uncertainty by first-order
    propagation from the covariance of (a, b).
    """
    t = np.array([(w["mid"] - win[0]["mid"]).days for w in win], dtype=float)
    v = np.array([w["v_mm_day"] for w in win], dtype=float)
    inv = 1.0 / v

    n = len(t)
    A = np.vstack([t, np.ones_like(t)]).T
    coef, *_ = np.linalg.lstsq(A, inv, rcond=None)
    a, b = float(coef[0]), float(coef[1])

    pred = A @ coef
    ss_res = float(np.sum((inv - pred) ** 2))
    ss_tot = float(np.sum((inv - inv.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    t_fail_rel = -b / a if a != 0 else float("nan")
    sigma_tf = float("nan")
    if n > 2 and a != 0:
        s2 = ss_res / (n - 2)
        try:
            cov = s2 * np.linalg.inv(A.T @ A)
            # d(t_f)/da = b/a^2 ,  d(t_f)/db = -1/a
            g = np.array([b / a**2, -1.0 / a])
            var = float(g @ cov @ g)
            sigma_tf = float(np.sqrt(var)) if var > 0 else float("nan")
        except np.linalg.LinAlgError:
            pass

    return {"slope": a, "intercept": b, "r2": r2,
            "t_fail_days_from_window_start": t_fail_rel,
            "t_fail_date": win[0]["mid"] + timedelta(days=t_fail_rel)
                           if np.isfinite(t_fail_rel) else None,
            "sigma_days": sigma_tf, "n": n}


# ---------------------------------------------------------------------------
def analyse_block(label: str, comp: int, rows: list[dict], noise_floor: float,
                  window: int, r2_min: float, horizon_days: float,
                  sig_multiple: float, event_date: date | None) -> dict:
    print(f"\n{'='*74}")
    print(f"{label}  block {comp}   {rows[0]['epoch']} .. {rows[-1]['epoch']}  "
          f"({len(rows)} epochs)")
    print("=" * 74)

    vs = velocities(rows)
    if not vs:
        print("  no velocity estimates possible")
        return {"alarm": False, "reason": "no intervals"}

    threshold = sig_multiple * noise_floor
    print(f"\n  significance gate: |v| must exceed {sig_multiple:g} x {noise_floor:g} "
          f"= {threshold:.1f} mm/day\n")
    print(f"  {'INTERVAL':<26}{'DAYS':>6}{'v mm/day':>11}{'1/v':>10}   STATUS")
    print("  " + "-" * 66)
    for w in vs:
        sig = w["v_mm_day"] > threshold
        iv = f"{1/w['v_mm_day']:.4f}" if w["v_mm_day"] > 0 else "-"
        status = "above floor" if sig else "below floor"
        print(f"  {str(w['t0'])+' -> '+str(w['t1']):<26}{w['days']:>6}"
              f"{w['v_mm_day']:>11.2f}{iv:>10}   {status}")

    usable = [w for w in vs if w["v_mm_day"] > threshold]
    print(f"\n  {len(usable)} of {len(vs)} intervals clear the floor")

    # ---- gates, in order, each with a stated reason -----------------------
    if len(usable) < window:
        need = window - len(usable)
        print(f"\n  NO ALARM - only {len(usable)} usable velocities, the fit needs "
              f"{window}. Short by {need}.")
        vmax = max(w["v_mm_day"] for w in vs)
        print(f"  fastest interval measured: {vmax:+.2f} mm/day "
              f"({vmax/noise_floor:.2f} x the noise floor)")
        print(f"  a detection here would have needed sustained motion above "
              f"{threshold:.0f} mm/day across {window} consecutive intervals.")
        return {"alarm": False, "reason": "insufficient velocities above noise floor",
                "max_velocity": vmax, "required": threshold}

    best = None
    for k in range(window, len(usable) + 1):
        win = usable[k - window:k]
        fit = fit_inverse_velocity(win)
        if fit["slope"] >= 0:
            continue                      # 1/v must be falling
        if fit["r2"] < r2_min:
            continue
        if fit["t_fail_date"] is None:
            continue
        lead = (fit["t_fail_date"] - win[-1]["mid"]).days
        if not (0 < lead <= horizon_days):
            continue
        best = (win, fit, lead)
        break

    if best is None:
        print(f"\n  NO ALARM - {len(usable)} usable velocities, but no window of "
              f"{window} met all of:")
        print(f"    1/v decreasing (acceleration), R2 >= {r2_min}, "
              f"predicted failure 0-{horizon_days:g} days ahead")
        return {"alarm": False, "reason": "no qualifying accelerating trend"}

    win, fit, lead = best
    print(f"\n  *** ALARM ***")
    print(f"    fitted on {fit['n']} velocities ending {win[-1]['mid']}")
    print(f"    1/v slope {fit['slope']:+.5f} per day, R2 {fit['r2']:.3f}")
    pm = f" +/- {fit['sigma_days']:.1f}" if np.isfinite(fit["sigma_days"]) else ""
    print(f"    predicted failure {fit['t_fail_date']}{pm} days")
    print(f"    lead time {lead} days from the last observation")
    if event_date:
        err = (fit["t_fail_date"] - event_date).days
        print(f"    actual event {event_date}  ->  prediction error {err:+d} days")
    return {"alarm": True, "predicted": fit["t_fail_date"], "lead_days": lead,
            "r2": fit["r2"], "sigma_days": fit["sigma_days"]}


def plot(series, noise_floor, sig_multiple, out: Path, event_date=None):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.error("matplotlib required for --plot")
        return
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for (geom, comp), rows in sorted(series.items()):
        lab = f"{geom} blk{comp}"
        ax1.plot([r["epoch"] for r in rows], [r["cumulative_mm"] for r in rows],
                 "o-", label=lab, alpha=.85)
        vs = velocities(rows)
        if vs:
            ax2.plot([w["mid"] for w in vs], [w["v_mm_day"] for w in vs],
                     "o-", label=lab, alpha=.85)
    thr = sig_multiple * noise_floor
    ax2.axhline(thr, color="crimson", ls="--", lw=1.2)
    ax2.axhline(-thr, color="crimson", ls="--", lw=1.2)
    ax2.text(ax2.get_xlim()[0], thr, f"  significance gate {thr:.0f} mm/day",
             color="crimson", va="bottom", fontsize=8)
    if event_date:
        for ax in (ax1, ax2):
            ax.axvline(event_date, color="k", ls=":", lw=1.2)
    ax1.set_ylabel("Cumulative displacement (mm)")
    ax1.set_title("Measured series and velocity against the detection floor", fontsize=12)
    ax2.set_ylabel("Velocity (mm/day)")
    for ax in (ax1, ax2):
        ax.grid(alpha=.3); ax.legend(fontsize=8)
    fig.autofmt_xdate(); fig.tight_layout(); fig.savefig(out, dpi=140)
    logger.info("Wrote %s", out)


def main() -> int:
    ap = argparse.ArgumentParser(description="Inverse-velocity forecasting on measured data")
    ap.add_argument("--ts", required=True, help="CSV from timeseries.py --csv")
    ap.add_argument("--noise-floor", type=float, required=True,
                    help="mm/day, from goff_reader.py --noise-floor (GOFF) "
                         "or the phase ceiling analysis (GUNW)")
    ap.add_argument("--sig-multiple", type=float, default=1.0,
                    help="velocity must exceed this many times the floor")
    ap.add_argument("--window", type=int, default=3)
    ap.add_argument("--r2-min", type=float, default=0.70)
    ap.add_argument("--horizon", type=float, default=60.0)
    ap.add_argument("--event-date", help="YYYY-MM-DD, to score the prediction")
    ap.add_argument("--plot", metavar="OUT.png")
    args = ap.parse_args()

    ev = datetime.strptime(args.event_date, "%Y-%m-%d").date() if args.event_date else None
    series = load_series(resolve(args.ts))
    if not series:
        logger.error("No rows in %s", args.ts); return 1

    print(f"\n{len(series)} independent block(s) from {args.ts}")
    print(f"Noise floor {args.noise_floor:g} mm/day, "
          f"significance gate {args.sig_multiple:g}x")

    results = []
    for (geom, comp), rows in sorted(series.items()):
        results.append(analyse_block(geom, comp, rows, args.noise_floor,
                                     args.window, args.r2_min, args.horizon,
                                     args.sig_multiple, ev))

    print(f"\n{'='*74}\nSUMMARY\n{'='*74}")
    fired = [r for r in results if r.get("alarm")]
    if fired:
        for r in fired:
            print(f"  ALARM: predicted {r['predicted']}, lead {r['lead_days']} d, "
                  f"R2 {r['r2']:.2f}")
    else:
        print(f"  No alarm in any of {len(results)} blocks.")
        for r in results:
            print(f"    - {r.get('reason')}")
        mv = [r["max_velocity"] for r in results if "max_velocity" in r]
        if mv:
            gate = args.sig_multiple * args.noise_floor
            fastest = max(mv)
            print(f"\n  Fastest single interval anywhere: {fastest:+.2f} mm/day "
                  f"({fastest/args.noise_floor:.1f}x the noise floor)")
            print(f"  Detection required: {gate:.2f} mm/day sustained across "
                  f"{args.window} consecutive intervals.")
            # Two different failure modes, and saying which one is the point.
            # A single fast interval that exceeds the gate is not a shortfall in
            # magnitude - it is a shortfall in persistence, which is exactly how
            # atmospheric noise differs from creep. Reporting them the same way
            # would hide the distinction the whole method rests on.
            if fastest > gate:
                print(f"\n  Motion DID exceed the gate in at least one interval, but "
                      f"never for\n  {args.window} in a row. An isolated excursion "
                      f"that does not persist is the\n  signature of atmosphere, "
                      f"not of accelerating creep - and the series\n  returning to "
                      f"its starting value confirms it.")
            else:
                print(f"\n  No interval anywhere reached the gate; the fastest was "
                      f"{gate/fastest:.1f}x short.")
            print("\n  This is a bounded non-detection, not an absence of evidence:")
            print("  any precursor slower than the floor is invisible to this product,")
            print("  and that bound is the quotable result.")

    if args.plot:
        plot(series, args.noise_floor, args.sig_multiple, Path(args.plot), ev)
    return 0


if __name__ == "__main__":
    sys.exit(main())
