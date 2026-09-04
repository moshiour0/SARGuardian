"""
Does the detector detect?

The suite had 23 tests and not one of them ran the alarm. Every real dataset
returned "no alarm", which is the correct answer for a null result and also
the output of a detector that cannot alarm at all - the two are
indistinguishable without a positive control.

A referee found the difference. The gate read `v > threshold` on a SIGNED
line-of-sight velocity, so a slope moving toward the satellite was discarded
entirely. Downslope motion projects negative in ascending on a west-facing
slope, which is the dominant configuration at this site (sensitivity -0.908),
so the detector was blind to precisely the case it exists for.

These tests build a textbook Fukuzono failure - 1/v falling linearly to zero -
and require the same answer in both directions.
"""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

from inverse_velocity import analyse_block, fit_inverse_velocity, velocities

ROOT = Path(__file__).resolve().parent.parent
FAIL_DAY = 60.0
A = 0.004                       # 1/v = A * (T - t)


def failure_series(sign: int, days=(0, 6, 12, 18, 24, 30, 36, 42, 48)):
    """Voight/Fukuzono acceleration, projected into LOS with a given sign."""
    from datetime import date, timedelta
    start = date(2026, 6, 1)
    rows, cum, prev = [], 0.0, None
    for d in days:
        v = 1.0 / (A * (FAIL_DAY - d))
        if prev is not None:
            cum += sign * v * (d - prev)
        prev = d
        rows.append({"epoch": start + timedelta(days=d), "cumulative_mm": cum,
                     "error_mm": None})
    return rows


def run(rows, noise_floor=0.5):
    return analyse_block("TEST", 1, rows, noise_floor=noise_floor, window=3,
                         r2_min=0.8, horizon_days=120, sig_multiple=1.0,
                         event_date=None)


# ---------------------------------------------------------------------------
def test_detector_alarms_on_a_textbook_failure():
    """The positive control that never existed. If this fails, nothing else matters."""
    r = run(failure_series(+1))
    assert r["alarm"], "no alarm on a perfect Fukuzono acceleration"
    assert r["r2"] > 0.99


@pytest.mark.parametrize("sign", [+1, -1])
def test_alarm_does_not_depend_on_look_direction(sign):
    """
    Away from the satellite or toward it, the ground is doing the same thing.
    Which sign that produces is a property of the orbit, not of the slope.
    """
    r = run(failure_series(sign))
    assert r["alarm"], (
        f"sign {sign:+d} produced no alarm - the gate is reading signed velocity")


def test_both_directions_predict_the_same_failure_date():
    """A sign convention must not move the answer."""
    a = run(failure_series(+1))
    b = run(failure_series(-1))
    assert a["predicted"] == b["predicted"]
    assert a["lead_days"] == b["lead_days"]


def test_reported_peak_is_the_largest_magnitude():
    """
    `max()` on signed velocities reports the largest motion AWAY from the
    satellite and silently discards everything moving toward it. On the
    committed GOFF summer series that understated the peak by a factor of ten,
    -15.86 mm/day reported as +1.55.
    """
    rows = [{"epoch": d, "cumulative_mm": c, "error_mm": None} for d, c in
            zip(*[[__import__("datetime").date(2026, 7, x) for x in (1, 13, 25)],
                  [0.0, -190.0, -170.0]])]
    r = run(rows, noise_floor=100.0)          # floor high enough to block the fit
    assert not r["alarm"]
    assert r["max_velocity"] == pytest.approx(-15.833, abs=0.01), (
        f"reported {r['max_velocity']:+.2f}; the largest magnitude is -15.83")


def test_inverse_velocity_fit_uses_speed_not_signed_velocity():
    """
    1/v built from a negative velocity rises toward zero instead of falling to
    it, so the `slope < 0` test rejects a real failure. The fit has to work in
    speed and carry direction separately.
    """
    for sign in (+1, -1):
        rows = failure_series(sign)
        vs = velocities(rows)
        fit = fit_inverse_velocity(vs[-3:])
        assert fit["slope"] < 0, f"sign {sign:+d}: 1/v is not falling"
        assert fit["r2"] > 0.99


def test_direction_reversal_is_not_treated_as_acceleration():
    """
    Taking |v| everywhere would let noise either side of zero look like a
    consistent run. A slope approaching failure does not reverse.
    """
    from datetime import date, timedelta
    start = date(2026, 7, 1)
    cum, rows = 0.0, []
    for k, v in enumerate([0.0, +30.0, -30.0, +30.0, -30.0]):
        cum += v
        rows.append({"epoch": start + timedelta(days=12 * k),
                     "cumulative_mm": cum, "error_mm": None})
    r = run(rows, noise_floor=1.0)
    assert not r["alarm"], "alternating noise produced an alarm"


def test_cli_runs_end_to_end(tmp_path):
    """The module is used from the command line; that path must work too."""
    p = tmp_path / "ts.csv"
    rows = failure_series(-1)
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["geometry", "component", "epoch", "days_from_start",
                    "cumulative_mm", "error_mm"])
        for i, r in enumerate(rows):
            w.writerow(["DESC path 48", 1, r["epoch"], i,
                        round(r["cumulative_mm"], 3), ""])
    out = subprocess.run([sys.executable, str(ROOT / "src" / "inverse_velocity.py"),
                          "--ts", str(p), "--noise-floor", "0.5"],
                         capture_output=True, text=True, cwd=ROOT)
    assert out.returncode == 0, out.stderr[-500:]
    assert "ALARM" in out.stdout


# ---------------------------------------------------------------------------
# Per-pair gating.
#
# One scalar floor across a stack whose per-pair floors vary twelvefold is not a
# threshold, it is an average of thresholds. Measured on the source zone, the
# descending interval 2026-06-29 -> 2026-07-11 reads +19.60 mm/day and clears a
# global 18.6 mm/day gate at 1.05x, while that pair's own 3-sigma floor is
# 114.3 mm/day - against which the same number is 0.17x, plainly inside the
# noise. A gate that manufactures excursions on the noisy geometry is worse
# than no gate, because it looks like a measurement.
# ---------------------------------------------------------------------------
import csv as _csv
from datetime import date as _date

from inverse_velocity import load_floors, velocities


def _series(pairs):
    """[(epoch, cumulative_mm), ...] -> rows in load_series shape."""
    return [{"epoch": e, "cumulative_mm": v, "error_mm": None} for e, v in pairs]


def test_each_interval_carries_the_floor_of_the_pair_that_produced_it():
    rows = _series([(_date(2026, 6, 29), 0.0),
                    (_date(2026, 7, 11), 235.1),
                    (_date(2026, 7, 23), 138.8)])
    floors = {(_date(2026, 6, 29), _date(2026, 7, 11)): 114.3,
              (_date(2026, 7, 11), _date(2026, 7, 23)): 86.2}
    vs = velocities(rows, floors)
    assert [w["floor"] for w in vs] == [114.3, 86.2]
    # and the first interval, which clears a 18.6 global gate, does not clear
    # its own
    assert abs(vs[0]["v_mm_day"]) > 18.6
    assert abs(vs[0]["v_mm_day"]) < vs[0]["floor"]


def test_an_interval_with_no_measured_floor_reports_none():
    rows = _series([(_date(2026, 6, 29), 0.0), (_date(2026, 7, 11), 10.0)])
    assert velocities(rows, {})[0]["floor"] is None
    assert velocities(rows)[0]["floor"] is None


def test_duplicate_processings_keep_the_larger_floor(tmp_path):
    """
    Routine and urgent processing of the same acquisitions both appear in the
    stats CSV. A bound must not be improved by reprocessing the same data, so
    the pessimistic floor wins.
    """
    p = tmp_path / "goff.csv"
    with open(p, "w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=["file", "layer", "reference",
                                            "secondary", "detect_floor_mm_day"])
        w.writeheader()
        w.writerow({"file": "..._PR_...", "layer": "HH/layer2",
                    "reference": "20260816", "secondary": "20260828",
                    "detect_floor_mm_day": "85.6"})
        w.writerow({"file": "..._UR_...", "layer": "HH/layer2",
                    "reference": "20260816", "secondary": "20260828",
                    "detect_floor_mm_day": "44.6"})
    got = load_floors(p)
    assert got[(_date(2026, 8, 16), _date(2026, 8, 28))] == 85.6
