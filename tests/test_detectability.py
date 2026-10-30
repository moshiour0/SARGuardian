"""
The significance gate is the result, not a detail of the detector.

Why this exists
===============
This module had no tests, and it carries the project's procurement conclusion -
the advice an agency would act on when deciding what revisit to buy. Two faults
lived in it undisturbed.

`--blatten` was declared, documented as a calibrated preset, and never read.
Passing it changed nothing and said nothing.

The larger one: the detector's threshold was a hardcoded 1.0 mm/day that was
not exposed on the command line and did not scale with noise. Velocity noise is
sigma*sqrt(2)/dt, so at 5 mm displacement noise it is 7.1 mm/day at daily
revisit and 0.6 mm/day at 12-day - the flat 1 mm/day gate sat an order of
magnitude BELOW the noise at short revisit and above it at long revisit. The
detector was therefore admitting pure noise as signal, and doing it worst
exactly where the published conclusion said daily sampling was dangerous.

Measured both ways on the same sweep:

    10-day precursor, 1-day revisit    fixed gate      noise gate
      detection                            90%             98%
      false alarm                          12.0%            0.0%
      premature alarm                       9%              0%
      prediction error                      4.2 d           0.7 d

So "daily revisit is worse, not better" was an artefact of a detector with no
significance gate - the very thing inverse_velocity.py exists to apply to
measured data. What actually survives is a trade, not a cliff: a noise-aware
gate is HIGHER at short dt, so a fast-sampling detector must wait for a faster
slope, buying accuracy and losing lead time.

These tests hold the gate's behaviour, the preset's, and the reproducibility of
a single cell.
"""

from __future__ import annotations

from contextlib import redirect_stdout
from datetime import date, timedelta
from io import StringIO
from types import SimpleNamespace
import warnings

import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from detectability import Detector, Precursor, simulate
from pixel_stack import _correlation_length_km, chain, parse_export, report

# ---------------------------------------------------------------------------
# GUNW pixel-stack network and spatial-noise gates
# ---------------------------------------------------------------------------
def _pixel_stack_geo(pixel_size_m: float = 200.0) -> dict:
    return {"transform": SimpleNamespace(a=pixel_size_m)}


def _sparse_pixel_positions():
    return [(5, 5), (15, 17), (25, 29), (35, 41), (45, 53), (55, 65),
            (65, 8), (8, 70), (70, 70), (40, 10), (10, 40), (30, 60)]


def test_gunw_chain_does_not_join_different_subswaths_or_crids():
    """Same track/look and consecutive dates do not make different networks compatible."""
    first = date(2026, 6, 8)
    second = date(2026, 6, 20)
    third = date(2026, 7, 2)
    rasters = [
        {
            "geom": "A098",
            "network_id": "A098-F016-SW023-BW4000-SH-X05026",
            "ref": first,
            "sec": second,
        },
        {
            "geom": "A098",
            "network_id": "A098-F016-SW024-BW4000-SH-P05023",
            "ref": second,
            "sec": third,
        },
    ]

    runs = chain(rasters)
    assert len(runs) == 2
    assert all(len(run) == 1 for run in runs)
    assert {run[0]["network_id"] for run in runs} == {
        "A098-F016-SW023-BW4000-SH-X05026",
        "A098-F016-SW024-BW4000-SH-P05023",
    }


def test_parse_export_preserves_complete_network_identity(tmp_path):
    """The GeoTIFF source tag must retain frame/swath/bandwidth/pol/CRID."""
    import rasterio
    from rasterio.transform import from_origin

    path = tmp_path / "GUNW_20260608_20260620_PR.tif"
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=2,
        width=2,
        count=2,
        dtype="float32",
        crs="EPSG:32645",
        transform=from_origin(500000, 3000000, 80, 80),
    ) as ds:
        ds.write(np.ones((2, 2), dtype="float32"), 1)
        ds.write(np.ones((2, 2), dtype="float32"), 2)
        ds.update_tags(
            source=(
                "NISAR_L2_PR_GUNW_022_098_A_016_023_4000_SH_"
                "20260608T233922_20260608T233956_20260620T233921_"
                "20260620T233956_P05023_N_F_J_001.h5"
            )
        )

    parsed = parse_export(path)
    assert parsed is not None
    assert parsed["geom"] == "A098"
    assert parsed["network_id"] == "A098-F016-SW023-BW4000-SH-P05023"


def test_sparse_complete_pixels_withhold_spatial_noise_correlation():
    """A handful of isolated valid pixels cannot estimate a spatial null."""
    velocity = np.full((80, 80), np.nan, dtype="float32")
    valid = np.zeros(velocity.shape, dtype=bool)
    for row, col in _sparse_pixel_positions():
        velocity[row, col] = 1.0
        valid[row, col] = True

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        assert _correlation_length_km(velocity, valid, _pixel_stack_geo()) is None


def test_dense_valid_pixels_can_estimate_spatial_correlation():
    rng = np.random.default_rng(23)
    velocity = rng.normal(0.0, 1.0, size=(80, 80)).astype("float32")
    valid = np.ones(velocity.shape, dtype=bool)

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        estimate = _correlation_length_km(velocity, valid, _pixel_stack_geo())

    assert estimate is not None
    assert estimate >= 0.2


def test_sparse_stack_verdict_is_not_testable_not_no_detection():
    """Withheld noise calibration must not be presented as evidence of no motion."""
    velocity = np.full((80, 80), np.nan, dtype="float32")
    sigma = np.full_like(velocity, np.nan)
    valid = np.zeros(velocity.shape, dtype=bool)
    observed = np.zeros(velocity.shape, dtype="uint16")
    for row, col in _sparse_pixel_positions():
        velocity[row, col] = 0.5
        sigma[row, col] = 0.01
        valid[row, col] = True
        observed[row, col] = 3

    first = date(2026, 6, 1)
    run = [
        {"ref": first, "sec": first + timedelta(days=12)},
        {"ref": first + timedelta(days=12), "sec": first + timedelta(days=24)},
        {"ref": first + timedelta(days=24), "sec": first + timedelta(days=36)},
    ]
    result = {
        "velocity": velocity,
        "sigma": sigma,
        "valid": valid,
        "run": run,
        "days": np.array([0.0, 12.0, 24.0, 36.0]),
        "observed_pair_count": observed,
        "common_mode_policy": "PAIRWISE_MEDIAN_REMOVED",
        "uncertainty_status": "RESIDUAL_BASED_DIAGNOSTIC_NOT_INDEPENDENTLY_VALIDATED",
        "geo": _pixel_stack_geo(),
    }

    output = StringIO()
    with redirect_stdout(output), warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        verdict = report(result, "A098 chain 1")

    rendered = output.getvalue()
    assert verdict is None
    assert result["spatial_null_test_status"] == "WITHHELD_INSUFFICIENT_SPATIAL_NEIGHBOR_SUPPORT"
    assert "VERDICT: NOT TESTABLE" in rendered
    assert "NO DETECTION" not in rendered
    assert "noise-only cluster test        : WITHHELD" in rendered


ROOT = Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------------------
# The gate
# ---------------------------------------------------------------------------
def test_noise_gate_scales_with_the_velocity_noise_it_guards():
    """
    A velocity is a difference of two noisy displacements over dt, so its
    1-sigma is sigma*sqrt(2)/dt. The gate must track that, which means it rises
    as revisit shortens - the whole reason short revisit is not free.
    """
    d = Detector(gate="noise", noise_mm=5.0, sig_multiple=1.0)
    for dt in (1.0, 2.0, 4.0, 12.0):
        assert d.threshold(dt) == pytest.approx(3 * 5.0 * math.sqrt(2) / dt)
    assert d.threshold(1.0) > d.threshold(12.0)
    assert d.threshold(1.0) / d.threshold(12.0) == pytest.approx(12.0)


def test_fixed_gate_ignores_noise_which_is_the_point_of_naming_it():
    d = Detector(gate="fixed", min_velocity_mm_day=1.0, noise_mm=5.0)
    assert d.threshold(1.0) == d.threshold(12.0) == 1.0


def test_the_flat_gate_sat_below_the_noise_at_short_revisit():
    """
    The specific arithmetic behind the artefact. At 5 mm noise the flat
    1 mm/day threshold is a tenth of the velocity noise at daily revisit, so
    it admitted essentially every sample.
    """
    noise_mm, dt = 5.0, 1.0
    velocity_sigma = noise_mm * math.sqrt(2) / dt
    assert Detector(gate="fixed").threshold(dt) < velocity_sigma / 5
    assert Detector(gate="noise", noise_mm=noise_mm).threshold(dt) > velocity_sigma


def test_a_slope_that_never_moves_does_not_alarm_under_the_noise_gate():
    """Pure noise, many trials, and the false-alarm rate must be negligible."""
    rng = np.random.default_rng(7)
    det = Detector(gate="noise", noise_mm=5.0)
    r = simulate(Precursor(10.0, 300.0), det, revisit_days=1.0, noise_mm=5.0,
                 n_trials=1, rng=rng, n_null=400)
    assert r["false_alarm_rate"] < 0.02, (
        f"false alarm rate {r['false_alarm_rate']:.1%} on a slope that never "
        "moved - the gate is not gating")


def test_the_flat_gate_reproduces_the_false_alarms_it_used_to_produce():
    """
    The comparison has to be runnable or the correction is just an assertion.
    Same noise, same cadence, gate the only difference.
    """
    rng = np.random.default_rng(7)
    loose = simulate(Precursor(10.0, 300.0), Detector(gate="fixed"),
                     revisit_days=1.0, noise_mm=5.0, n_trials=1,
                     rng=rng, n_null=400)
    rng = np.random.default_rng(7)
    tight = simulate(Precursor(10.0, 300.0), Detector(gate="noise", noise_mm=5.0),
                     revisit_days=1.0, noise_mm=5.0, n_trials=1,
                     rng=rng, n_null=400)
    assert loose["false_alarm_rate"] > tight["false_alarm_rate"]
    assert loose["false_alarm_rate"] > 0.02


def test_a_real_acceleration_still_detected_with_the_gate_on():
    """The positive control. A gate that never fires is not an improvement."""
    rng = np.random.default_rng(3)
    r = simulate(Precursor(10.0, 300.0), Detector(gate="noise", noise_mm=5.0),
                 revisit_days=2.0, noise_mm=5.0, n_trials=200, rng=rng, n_null=50)
    assert r["detection_rate"] > 0.5
    assert r["warning_median"] > 0


# ---------------------------------------------------------------------------
# The preset that did nothing
# ---------------------------------------------------------------------------
def _run(*args) -> str:
    return subprocess.run(
        [sys.executable, str(ROOT / "src" / "detectability.py"), *args],
        capture_output=True, text=True, check=True).stdout


def test_blatten_preset_actually_changes_the_run():
    """
    It was declared, documented, and never read: `args.blatten` appeared in no
    expression. Passing it produced a run on the defaults, silently.
    """
    out = _run("--sweep", "--blatten", "--noise", "300",
               "--trials", "5", "--null-trials", "5")
    assert "PRECURSOR 7 days" in out
    assert "27000 mm creep" in out
    default = _run("--sweep", "--precursor", "40", "--revisit", "12",
                   "--trials", "5", "--null-trials", "5")
    assert "PRECURSOR 7 days" not in default


def test_an_explicit_flag_beats_the_preset():
    """A preset that overrides what the user typed is worse than no preset."""
    out = _run("--sweep", "--blatten", "--precursor", "20", "--revisit", "4",
               "--noise", "300", "--trials", "5", "--null-trials", "5")
    assert "PRECURSOR 20 days" in out
    assert "PRECURSOR 7 days" not in out


# ---------------------------------------------------------------------------
# Reproducibility of one cell
# ---------------------------------------------------------------------------
def test_a_cell_reproduces_regardless_of_what_it_was_swept_alongside():
    """
    The generator used to be created once and consumed in sequence across the
    whole sweep, so a cell's numbers depended on which --revisit list it
    happened to travel with and no single cell could be reproduced on its own.
    """
    alone = _run("--sweep", "--precursor", "10", "--revisit", "4",
                 "--trials", "40", "--null-trials", "40")
    together = _run("--sweep", "--precursor", "10", "--revisit", "1", "2", "4",
                    "--trials", "40", "--null-trials", "40")

    def row(text, dt="4"):
        for line in text.splitlines():
            if line.strip().startswith(dt + " "):
                return line.split()
        raise AssertionError(f"no revisit {dt} row in:\n{text}")

    assert row(alone)[:5] == row(together)[:5]


# ---------------------------------------------------------------------------
# Forward model and the phase ceiling
# ---------------------------------------------------------------------------
def test_precursor_velocity_is_the_derivative_of_its_own_displacement():
    p = Precursor(10.0, 300.0)
    t = np.linspace(0.0, 9.0, 200)
    numeric = np.gradient(p.displacement(t), t)
    assert np.allclose(numeric[5:-5], p.velocity(t)[5:-5], rtol=0.02)


def test_creep_amplitude_is_what_was_asked_for():
    """x(T-1) = creep_to_1d_mm, by construction. A wrong A silently rescales."""
    for T, creep in ((10.0, 300.0), (7.0, 27000.0), (40.0, 300.0)):
        p = Precursor(T, creep)
        assert float(p.displacement(np.array([T - 1.0]))[0]) == pytest.approx(creep, rel=1e-6)


def test_gradient_fraction_relaxes_the_phase_ceiling():
    """
    lambda/4 limits the phase difference between ADJACENT PIXELS, not the
    displacement of one pixel between passes. Applying it to the temporal step
    assumes the whole step falls across one cell boundary, which is the most
    pessimistic reading. A gentler gradient must saturate less often.
    """
    def sat(frac):
        rng = np.random.default_rng(11)
        return simulate(Precursor(10.0, 300.0), Detector(gate="noise", noise_mm=5.0),
                        revisit_days=2.0, noise_mm=5.0, n_trials=120, rng=rng,
                        wavelength_m=0.2384, n_null=10,
                        gradient_fraction=frac)["saturation_rate"]

    # 300 mm over 10 days gives inter-pass steps of 26 to 260 mm at 2-day
    # revisit, against a 59.6 mm L-band quarter-wavelength - so the pessimistic
    # reading saturates and a tenth of the gradient largely does not.
    assert sat(1.0) > 0.9, "the pessimistic case must saturate"
    assert sat(0.1) < 0.5 * sat(1.0), (
        "a gradual displacement gradient must saturate phase far less often "
        "than a scarp that abuts stable ground within one pixel")

    # Blatten is so fast that no plausible gradient saves phase, which is the
    # half of the argument that does survive.
    rng = np.random.default_rng(11)
    blatten = simulate(Precursor(7.0, 27000.0), Detector(gate="noise", noise_mm=5.0),
                       revisit_days=1.0, noise_mm=5.0, n_trials=40, rng=rng,
                       wavelength_m=0.2384, n_null=5,
                       gradient_fraction=0.02)["saturation_rate"]
    assert blatten > 0.9, (
        "at 27 m of creep in 7 days phase saturates even at 2% of the gradient")


# ---------------------------------------------------------------------------
# Rates must divide by the trials they describe.
# ---------------------------------------------------------------------------
def test_skipped_trials_leave_the_denominator():
    """
    A cadence that puts too few samples in the record never reaches the
    detector, and those trials used to stay in the denominator anyway. At a
    7-day precursor, 12-day revisit and a 30-day lead-in that is 11 trials in
    12, so a saturation rate of ~100% was published as ~9% - and it moved with
    --lead-in, a parameter that changes no physics.
    """
    from detectability import Detector, Precursor, simulate
    rng = np.random.default_rng(1)
    r = simulate(Precursor(7.0, 27000.0), Detector(noise_mm=5.0), 12.0,
                 5.0, 200, rng, lead_in_days=30.0, wavelength_m=0.2384,
                 n_null=100)
    assert r["n_ran"] <= r["n_eligible"] <= r["n_trials"]
    assert r["n_eligible"] < r["n_trials"]          # trials really are dropped
    # the rate describes the trials that ran, so it is not diluted to ~1/12
    assert r["saturation_rate"] > 0.5


def test_the_window_count_scales_inversely_with_revisit():
    """
    Every trailing window is a separate opportunity to alarm. The count goes
    as 1/revisit, so a per-window criterion applied 37 times at daily sampling
    is not the same test as one applied once at 12-day. The multiplicity is
    reported rather than left implicit.
    """
    from detectability import Detector, Precursor, simulate
    rng = np.random.default_rng(2)
    fast = simulate(Precursor(10.0), Detector(noise_mm=5.0), 1.0, 5.0, 50, rng,
                    n_null=50)
    slow = simulate(Precursor(10.0), Detector(noise_mm=5.0), 4.0, 5.0, 50, rng,
                    n_null=50)
    assert fast["windows_per_trial"] > 3 * slow["windows_per_trial"]


def test_the_detector_reports_how_many_windows_it_tested():
    """Without this the per-window false-alarm rate cannot be formed at all."""
    from detectability import Detector
    d = Detector(noise_mm=5.0)
    t = np.arange(0.0, 20.0, 1.0)
    out = d.run(t, np.zeros_like(t), failure_day=20.0)
    assert out["tested"] == len(t) - 1 - d.window + 1


def test_too_few_samples_reports_zero_windows_tested():
    from detectability import Detector
    d = Detector(noise_mm=5.0)
    t = np.array([0.0, 1.0])
    assert d.run(t, np.zeros_like(t), failure_day=5.0)["tested"] == 0
