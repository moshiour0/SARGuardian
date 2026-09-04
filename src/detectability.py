"""
detectability.py
----------------
How much warning does a slope-failure early-warning system actually get,
as a function of satellite revisit interval?

This is a simulation. It downloads nothing and needs no SAR data.

Method
======
1. Model accelerating creep to failure with the Voight relation

       dv/dt = A * v^alpha

   For alpha = 2 - the usual value for brittle slope failure - this
   integrates to a straight line in inverse velocity:

       1/v(t) = A * (T - t)          v(t) = 1 / (A (T - t))
       x(t)   = (1/A) * ln( T / (T - t) )

   so 1/v falls linearly to zero at the failure time T. That is the
   Fukuzono construction, used operationally in open-pit mine slope
   monitoring, and it is what turns a displacement curve into a
   *predicted failure date*.

2. Sample that curve at a fixed revisit interval, adding realistic InSAR
   line-of-sight noise.

3. Run an inverse-velocity detector on the samples: fit 1/v against time
   over a trailing window, extrapolate to 1/v = 0, and raise an alarm when
   the fit is good, the trend is accelerating, and the predicted failure is
   in the near future.

4. Repeat over many noise realisations and many revisit intervals, and
   report how much warning time survives.

The competing effects, which is why the answer is not obvious
============================================================
Velocity is estimated as dx/dt between consecutive samples, so its noise
scales as sigma*sqrt(2)/dt. A LONGER revisit therefore gives a QUIETER
velocity estimate. But it also gives FEWER samples inside the precursory
window, and the detector needs a minimum number of them to fit a trend.
Those pull in opposite directions and the optimum depends on how long the
precursor lasts.

Calibration warning
===================
The absolute warning times below are only as good as the assumed precursor
duration, total creep, and measurement noise. Defaults are plausible for a
Blatten-class alpine glacier/rock failure but are NOT fitted to that event.
Calibrate against a documented record before quoting any number.

Usage
-----
    python src/detectability.py --sweep --plot outputs/detectability.png
    python src/detectability.py --sweep --precursor 10 20 40 --csv outputs/det.csv
    python src/detectability.py --demo --revisit 4
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("detectability")


# ---------------------------------------------------------------------------
# Forward model
# ---------------------------------------------------------------------------
@dataclass
class Precursor:
    """
    Accelerating creep to failure.

    duration_days   time from the onset of detectable acceleration to failure
    creep_to_1d_mm  cumulative displacement from onset until one day before
                    failure. Sets the amplitude scale.
    """
    duration_days: float = 10.0
    creep_to_1d_mm: float = 300.0

    @property
    def A(self) -> float:
        # x(T-1) = ln(T)/A  ->  A = ln(T)/creep
        return float(np.log(self.duration_days) / self.creep_to_1d_mm)

    def displacement(self, t: np.ndarray) -> np.ndarray:
        """Cumulative displacement (mm) at times t (days since onset)."""
        T = self.duration_days
        t = np.clip(t, 0.0, T - 1e-6)
        return np.log(T / (T - t)) / self.A

    def velocity(self, t: np.ndarray) -> np.ndarray:
        T = self.duration_days
        t = np.clip(t, 0.0, T - 1e-6)
        return 1.0 / (self.A * (T - t))


# ---------------------------------------------------------------------------
# Detector
# ---------------------------------------------------------------------------
@dataclass
class Detector:
    window: int = 3              # velocity estimates in the trailing fit
    r2_min: float = 0.70         # goodness of the 1/v trend
    horizon_days: float = 30.0   # ignore predictions further out than this
    min_velocity_mm_day: float = 1.0   # ignore velocities buried in noise

    def run(self, times: np.ndarray, disp: np.ndarray, failure_day: float) -> dict:
        """
        Returns the first alarm and the warning time it bought.

        times/disp are the sampled observations. failure_day is the truth,
        used only to score the result.
        """
        if len(times) < self.window + 1:
            return {"alarm": False, "reason": "too few samples"}

        dt = np.diff(times)
        vel = np.diff(disp) / dt
        tau = times[:-1] + dt / 2.0        # velocity is a mid-interval quantity

        for k in range(self.window, len(vel) + 1):
            wt = tau[k - self.window:k]
            wv = vel[k - self.window:k]

            if np.any(wv < self.min_velocity_mm_day):
                continue                    # not yet moving measurably

            inv = 1.0 / wv
            slope, intercept = np.polyfit(wt, inv, 1)
            if slope >= 0:
                continue                    # 1/v must be falling

            pred = -intercept / slope
            now = wt[-1]
            if not (now < pred <= now + self.horizon_days):
                continue

            fit = slope * wt + intercept
            ss_res = float(np.sum((inv - fit) ** 2))
            ss_tot = float(np.sum((inv - inv.mean()) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
            if r2 < self.r2_min:
                continue

            return {
                "alarm": True,
                "alarm_day": float(now),
                "warning_days": float(failure_day - now),
                "predicted_failure_day": float(pred),
                "prediction_error_days": float(pred - failure_day),
                "r2": r2,
                "n_samples_used": int(self.window),
            }

        return {"alarm": False, "reason": "no qualifying trend before failure"}


# ---------------------------------------------------------------------------
# Experiment
# ---------------------------------------------------------------------------
def simulate(precursor: Precursor, detector: Detector, revisit_days: float,
             noise_mm: float, n_trials: int, rng: np.random.Generator,
             lead_in_days: float = 30.0, wavelength_m: float | None = None,
             n_null: int = 1000) -> dict:
    """
    One (precursor, revisit) cell of the sweep.

    lead_in_days of stable ground before onset, so the detector has to cope
    with a quiet period first rather than starting conveniently at onset.

    wavelength_m, if given, imposes the PHASE UNWRAPPING CEILING. Interferometric
    phase is ambiguous once displacement between passes exceeds lambda/4, so the
    measurement is lost exactly when the slope accelerates hardest. Blatten
    reached 0.65 m/day six days before failure, which is 11x above even the
    1-day L-band ceiling - so phase saturates long before failure and only
    offset tracking survives. Pass None to model offset tracking (no ceiling,
    but use a much larger noise_mm).
    """
    T = precursor.duration_days
    ceiling_mm = (wavelength_m / 4.0) * 1000.0 if wavelength_m else None
    warnings, errors, alarms, saturated, premature = [], [], 0, 0, 0

    for _ in range(n_trials):
        # acquisitions on a fixed cadence with a random phase relative to onset
        offset = rng.uniform(0, revisit_days)
        t = np.arange(-lead_in_days + offset, T, revisit_days)
        t = t[t < T]                       # nothing observed after failure
        if len(t) < detector.window + 1:
            continue

        truth = np.where(t < 0, 0.0, precursor.displacement(np.maximum(t, 0.0)))

        if ceiling_mm is not None:
            # First interval whose true displacement exceeds lambda/4 aliases.
            # Motion only accelerates, so everything after it is lost too.
            step = np.abs(np.diff(truth))
            bad = np.nonzero(step > ceiling_mm)[0]
            if len(bad):
                cut = int(bad[0]) + 1
                t, truth = t[:cut], truth[:cut]
                saturated += 1
                if len(t) < detector.window + 1:
                    continue

        obs = truth + rng.normal(0.0, noise_mm, size=len(t))

        result = detector.run(t, obs, failure_day=T)
        if result["alarm"]:
            # An alarm inside the stable lead-in is a FALSE alarm, whatever the
            # trial contains later. The detector scans from the earliest window
            # and returns the first qualifying trend, and the earliest windows
            # sit before onset where the truth is identically zero. Counting
            # those as detections inflates the rate AND the warning time, since
            # warning is scored as T - now and `now` is negative. The tell is a
            # sweep where warning time RISES as detection collapses.
            if result["alarm_day"] < 0.0:
                premature += 1
            else:
                alarms += 1
                warnings.append(result["warning_days"])
                errors.append(result["prediction_error_days"])

    # A detection rate with no false-alarm rate is not a performance figure: a
    # detector that always fires scores 100%. The null costs almost nothing to
    # measure, because every trial already carries a stable lead-in - so run
    # the same detector, same cadence, same noise, on a slope that never moves.
    false_alarms = 0
    for _ in range(n_null):
        offset = rng.uniform(0, revisit_days)
        t = np.arange(-lead_in_days + offset, T, revisit_days)
        t = t[t < T]
        if len(t) < detector.window + 1:
            continue
        obs = rng.normal(0.0, noise_mm, size=len(t))
        if detector.run(t, obs, failure_day=T)["alarm"]:
            false_alarms += 1

    sat_rate = saturated / n_trials
    far = false_alarms / n_null if n_null else float("nan")
    base = {"revisit_days": revisit_days, "saturation_rate": sat_rate,
            "n_trials": n_trials, "false_alarm_rate": far, "n_null": n_null,
            "premature_rate": premature / n_trials}
    if alarms == 0:
        return {**base, "detection_rate": 0.0, "warning_median": float("nan"),
                "warning_p25": float("nan"), "warning_p75": float("nan"),
                "abs_pred_error_median": float("nan")}

    w = np.array(warnings)
    return {
        **base,
        "detection_rate": alarms / n_trials,
        "warning_median": float(np.median(w)),
        "warning_p25": float(np.percentile(w, 25)),
        "warning_p75": float(np.percentile(w, 75)),
        "abs_pred_error_median": float(np.median(np.abs(errors))),
    }


REVISITS = [1, 2, 3, 4, 6, 8, 12, 16, 24]


def sweep(precursor_days: list[float], revisits: list[float], noise_mm: float,
          creep_mm: float, detector: Detector, n_trials: int, seed: int,
          wavelength_m: float | None = None, n_null: int = 1000) -> list[dict]:
    rng = np.random.default_rng(seed)
    rows = []
    for T in precursor_days:
        p = Precursor(duration_days=T, creep_to_1d_mm=creep_mm)
        print(f"\n{'='*74}")
        mode = (f"PHASE, lambda={wavelength_m*100:.1f}cm" if wavelength_m
                else "OFFSET TRACKING, no ceiling")
        print(f"PRECURSOR {T:g} days   ({creep_mm:g} mm creep to 1 day before failure, "
              f"noise {noise_mm:g} mm)")
        print(f"MEASUREMENT: {mode}")
        print(f"{'='*74}")
        print(f"{'REVISIT':>8}{'SAMPLES':>8}{'SATUR':>7}{'DETECT':>8}{'FALSE':>7}"
              f"{'EARLY':>7}{'WARNING (days)':>22}{'|ERR|':>8}")
        print(f"{'days':>8}{'in T':>8}{'rate':>7}{'rate':>8}{'ALARM':>7}"
              f"{'rej.':>7}{'median [p25-p75]':>22}{'days':>8}")
        print("-" * 82)
        for dt in revisits:
            r = simulate(p, detector, dt, noise_mm, n_trials, rng,
                         wavelength_m=wavelength_m, n_null=n_null)
            r["precursor_days"] = T
            r["noise_mm"] = noise_mm
            r["creep_mm"] = creep_mm
            rows.append(r)
            n_in = int(T // dt)
            rate = r["detection_rate"]
            sat = f"{r.get('saturation_rate', 0):.0%}"
            far = f"{r['false_alarm_rate']:.1%}"
            early = f"{r['premature_rate']:.0%}"
            if rate == 0:
                print(f"{dt:>8g}{n_in:>8}{sat:>7}{'never':>8}{far:>7}{early:>7}"
                      f"{'-':>22}{'-':>8}")
            else:
                # a handful of lucky alarms is not detection - do not round it to 0%
                shown = "<1%" if rate < 0.005 else f"{rate:.0%}"
                band = f"{r['warning_median']:.1f} [{r['warning_p25']:.1f}-{r['warning_p75']:.1f}]"
                print(f"{dt:>8g}{n_in:>8}{sat:>7}{shown:>8}{far:>7}{early:>7}"
                      f"{band:>22}{r['abs_pred_error_median']:>8.1f}")
        # where does it stop working
        # Detection is only real when it beats its own false-alarm rate. A cell
        # detecting 10% against an 11.9% null has found nothing, however
        # respectable 10% looks in isolation.
        cells = [r for r in rows if r["precursor_days"] == T]
        useless = [r for r in cells
                   if r["detection_rate"] <= r["false_alarm_rate"]]
        if useless:
            print("\n  Detection at or below the false-alarm rate - no signal at all:")
            for r in useless:
                print(f"    revisit {r['revisit_days']:g} d: "
                      f"{r['detection_rate']:.0%} detection vs "
                      f"{r['false_alarm_rate']:.1%} false alarm")
        dead = [r["revisit_days"] for r in cells if r["detection_rate"] < 0.5]
        if dead:
            limit = min(dead)
            print(f"\n  Detection collapses at revisit >= {limit:g} days "
                  f"(= precursor / {T/limit:.1f}).")
            print(f"  Usable revisit for a {T:g}-day precursor: <= {T/3:.1f} days.")
        else:
            print("\n  Detected at better than 50% for every revisit tested.")
    return rows


# ---------------------------------------------------------------------------
def demo(precursor: Precursor, detector: Detector, revisit: float,
         noise_mm: float, seed: int) -> None:
    """Single realisation, printed sample by sample - useful for a talk."""
    rng = np.random.default_rng(seed)
    T = precursor.duration_days
    t = np.arange(-20.0, T, revisit)
    t = t[t < T]
    truth = np.where(t < 0, 0.0, precursor.displacement(np.maximum(t, 0.0)))
    obs = truth + rng.normal(0, noise_mm, len(t))

    print(f"\nPrecursor {T:g} d, revisit {revisit:g} d, noise {noise_mm:g} mm")
    print(f"{'DAY':>8}{'TRUE mm':>11}{'OBS mm':>11}{'v mm/d':>10}{'1/v':>10}")
    print("-" * 50)
    dt = np.diff(t); vel = np.diff(obs) / dt; tau = t[:-1] + dt / 2
    print(f"{t[0]:>8.1f}{truth[0]:>11.1f}{obs[0]:>11.1f}{'-':>10}{'-':>10}")
    for i in range(len(vel)):
        iv = f"{1/vel[i]:.3f}" if vel[i] > 0 else "-"
        print(f"{t[i+1]:>8.1f}{truth[i+1]:>11.1f}{obs[i+1]:>11.1f}{vel[i]:>10.2f}{iv:>10}")

    r = detector.run(t, obs, failure_day=T)
    print()
    if r["alarm"]:
        print(f"ALARM on day {r['alarm_day']:+.1f}  ->  {r['warning_days']:.1f} days of warning")
        print(f"  predicted failure day {r['predicted_failure_day']:+.2f} "
              f"(true {T:+.2f}, error {r['prediction_error_days']:+.2f} d), R2 {r['r2']:.3f}")
    else:
        print(f"NO ALARM before failure - {r['reason']}")


def plot(rows: list[dict], out: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.error("matplotlib required for --plot:  pip install matplotlib")
        return

    groups: dict[float, list[dict]] = {}
    for r in rows:
        groups.setdefault(r["precursor_days"], []).append(r)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 9), sharex=True)
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(groups)))

    for c, (T, rs) in zip(colors, sorted(groups.items())):
        rs = sorted(rs, key=lambda r: r["revisit_days"])
        x = [r["revisit_days"] for r in rs]
        ax1.plot(x, [r["detection_rate"] * 100 for r in rs], "o-", color=c,
                 label=f"{T:g}-day precursor")
        # Only draw warning time where detection is actually reliable. Below
        # 50% the median is computed from a handful of lucky alarms and would
        # imply a capability that does not exist.
        ok = [r for r in rs if r["detection_rate"] >= 0.5]
        if ok:
            xo = [r["revisit_days"] for r in ok]
            ax2.plot(xo, [r["warning_median"] for r in ok], "o-", color=c,
                     label=f"{T:g}-day precursor")
            ax2.fill_between(xo, [r["warning_p25"] for r in ok],
                             [r["warning_p75"] for r in ok], color=c, alpha=.15)
            ax2.plot(xo[-1], ok[-1]["warning_median"], "x", color=c,
                     ms=11, mew=2.5)

    ax1.axhline(50, color="grey", ls=":", lw=1)
    ax1.set_ylabel("Detection rate (%)")
    ax1.set_title("Slope-failure early warning vs satellite revisit interval\n"
                  "inverse-velocity detector on simulated accelerating creep",
                  fontsize=12)
    ax1.grid(alpha=.3); ax1.legend(fontsize=9); ax1.set_ylim(-5, 105)

    for dt, lab in [(1, "commercial\ntasking"), (4, "NISAR+S1\ncombined"), (12, "single\ntrack")]:
        for ax in (ax1, ax2):
            ax.axvline(dt, color="crimson", ls="--", lw=1, alpha=.5)
        ax1.text(dt, 103, lab, color="crimson", fontsize=8, ha="center", va="top")

    ax2.set_ylabel("Warning time before failure (days)\nmedian, shaded = IQR")
    ax2.set_xlabel("Revisit interval (days)")
    ax2.grid(alpha=.3); ax2.legend(fontsize=9)
    ax2.set_xscale("log"); ax2.set_xticks(REVISITS)
    ax2.set_xticklabels([str(r) for r in REVISITS])

    fig.tight_layout(); fig.savefig(out, dpi=140)
    logger.info("Wrote %s", out)


def main() -> int:
    ap = argparse.ArgumentParser(description="Warning time vs satellite revisit interval")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--sweep", action="store_true")
    mode.add_argument("--demo", action="store_true")
    ap.add_argument("--precursor", type=float, nargs="+", default=[5, 10, 20, 40],
                    help="precursor durations in days")
    ap.add_argument("--revisit", type=float, nargs="+", default=REVISITS)
    ap.add_argument("--noise", type=float, default=5.0, help="LOS noise, mm 1-sigma")
    ap.add_argument("--creep", type=float, default=300.0,
                    help="mm of creep from onset to 1 day before failure")
    ap.add_argument("--window", type=int, default=3, help="velocity samples in the fit")
    ap.add_argument("--r2-min", type=float, default=0.70)
    ap.add_argument("--trials", type=int, default=400)
    ap.add_argument("--null-trials", type=int, default=1000,
                    help="zero-signal trials per cell, for the false-alarm rate")
    ap.add_argument("--seed", type=int, default=20260829)
    ap.add_argument("--wavelength", type=float, default=None,
                    help="radar wavelength in m; imposes the lambda/4 phase ceiling. "
                         "0.2384 = NISAR L-band, 0.0555 = Sentinel-1 C-band. "
                         "Omit to model offset tracking (no ceiling).")
    ap.add_argument("--blatten", action="store_true",
                    help="calibrated Blatten preset: 7-day rapid phase, 10 m/day at failure")
    ap.add_argument("--csv", metavar="OUT.csv")
    ap.add_argument("--plot", metavar="OUT.png")
    args = ap.parse_args()

    det = Detector(window=args.window, r2_min=args.r2_min)

    if args.demo:
        demo(Precursor(args.precursor[0], args.creep), det,
             args.revisit[0], args.noise, args.seed)
        return 0

    rows = sweep(args.precursor, args.revisit, args.noise, args.creep,
                 det, args.trials, args.seed, args.wavelength,
                 n_null=args.null_trials)

    print("\n" + "=" * 74)
    print("These are MODEL results. Absolute warning times depend on the assumed")
    print("precursor duration, creep amplitude and measurement noise. Calibrate")
    print("against a documented failure before quoting any number.")
    print("=" * 74)

    if args.csv:
        keys = ["precursor_days", "revisit_days", "detection_rate",
                "false_alarm_rate", "premature_rate", "warning_median",
                "warning_p25", "warning_p75", "abs_pred_error_median",
                "noise_mm", "creep_mm", "n_trials"]
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
            w.writeheader(); w.writerows(rows)
        logger.info("Wrote %s (%d rows)", args.csv, len(rows))

    if args.plot:
        plot(rows, Path(args.plot))
    return 0


if __name__ == "__main__":
    sys.exit(main())
