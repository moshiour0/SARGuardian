"""
mutate.py
---------
Put each historical bug back and check the suite notices.

Why this exists
===============
A test that passes proves nothing on its own - it may assert something that
was never in danger. The only evidence a regression test works is that it
fails when the regression returns.

So this reintroduces each real fault by editing the source, runs the suite,
and checks that the test written for that fault is among the failures. Every
mutation is reverted in a finally block, and the run refuses to start if the
working tree is dirty, because a crash mid-mutation would otherwise leave a
sabotaged reader on disk. That nearly happened once; hence the guard.

    python tests/mutate.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"

# (label, file, find, replace, test that must fail)
MUTATIONS = [
    ("reference may cross connected components",
     "gunw_reader.py",
     "comp_ok = (~cut_v | (cut_c == target_comp)).all(axis=(1, 3))",
     "comp_ok = np.ones_like(full)",
     "test_reference_never_crosses_a_connected_component"),

    ("mask polarity inverted",
     "gunw_reader.py",
     "usable = (water == 0) & (ref_sub > 0) & (sec_sub > 0) & (losmask != 255)",
     "usable = (losmask == 0)",
     "test_mask_polarity_rejects_water_and_missing_subswaths"),

    ("reference requires a fully-valid block",
     "gunw_reader.py",
     "for level in (1.0, 0.8, 0.6, 0.4):",
     "for level in (1.0,):",
     "test_reference_survives_a_partly_invalid_scene"),

    ("displacement sign flipped",
     "gunw_reader.py",
     "scale = -(wavelength / (4.0 * math.pi)) * 1000.0",
     "scale = (wavelength / (4.0 * math.pi)) * 1000.0",
     "test_round_trip_holds_across_sign_and_magnitude"),

    ("ionosphere screen not removed",
     "gunw_reader.py",
     "phase = phase - np.nan_to_num(iono, nan=0.0)",
     "phase = phase",
     "test_ionosphere_screen_is_subtracted"),

    ("consistency check reads the wrong keys",
     "gunw_reader.py",
     'groups.setdefault((r["reference"], r["secondary"]), []).append(r)',
     'groups.setdefault((r["reference_date"], r["secondary_date"]), []).append(r)',
     "test_consistency_check_runs_on_reader_output_keys"),

    ("broken products pass verification",
     "gunw_reader.py",
     "            good.append(f)",
     "            good.append(f)\n        except Exception:\n            good.append(f)",
     "test_truncated_product_is_reported_not_skipped"),

    ("velocity gate reads sign, not magnitude",
     "inverse_velocity.py",
     'usable = [w for w in vs if abs(w["v_mm_day"]) > w["gate"]]',
     'usable = [w for w in vs if w["v_mm_day"] > w["gate"]]',
     "test_alarm_does_not_depend_on_look_direction"),

    ("inverse-velocity fit uses signed velocity",
     "inverse_velocity.py",
     'v = np.abs(np.array([w["v_mm_day"] for w in win], dtype=float))',
     'v = np.array([w["v_mm_day"] for w in win], dtype=float)',
     "test_inverse_velocity_fit_uses_speed_not_signed_velocity"),

    ("reported peak is the signed maximum",
     "inverse_velocity.py",
     'wmax = max(vs, key=lambda w: abs(w["v_mm_day"]))',
     'wmax = max(vs, key=lambda w: w["v_mm_day"])',
     "test_reported_peak_is_the_largest_magnitude"),

    ("duplicate products averaged instead of resolved",
     "timeseries.py",
     "out.append(best[1])",
     "best[1].value = sum(p.value for _, p in cands) / len(cands)\n        out.append(best[1])",
     "test_duplicate_products_are_resolved_not_averaged"),

    ("export grid clipped to valid pixels, not to the AOI",
     "gunw_reader.py",
     "    if col_off > 0.01 or row_off > 0.01:",
     "    if False:",
     "test_a_product_off_the_lattice_is_refused_not_resampled"),

    ("AOI grid not anchored to the absolute lattice",
     "gunw_reader.py",
     "    x0 = math.floor(min(rx) / res_x) * res_x - pad_px * res_x",
     "    x0 = min(rx) - pad_px * res_x",
     "test_two_products_covering_different_ground_get_the_same_grid"),

    ("per-pair floor dropped, gate falls back to the scalar",
     "inverse_velocity.py",
     '"floor": (floors or {}).get((a["epoch"], b["epoch"])),',
     '"floor": None,',
     "test_each_interval_carries_the_floor_of_the_pair_that_produced_it"),

    ("duplicate processing improves the bound",
     "inverse_velocity.py",
     "floors[(a, b)] = max(f, floors.get((a, b), 0.0))",
     "floors[(a, b)] = f",
     "test_duplicate_processings_keep_the_larger_floor"),

    ("geometry test decides the wrong way round",
     "candidate_check.py",
     "    if d_motion < d_delay:",
     "    if d_delay < d_motion:",
     "test_same_signed_rates_are_a_path_delay_where_sensitivities_oppose"),

    ("blind track still gets a verdict",
     "candidate_check.py",
     "    if abs(sens_asc) < MIN_SENS or abs(sens_desc) < MIN_SENS:",
     "    if False:",
     "test_a_blind_track_yields_no_verdict"),

    ("gentle terrain still gets a verdict",
     "candidate_check.py",
     "    if slope_deg is not None and slope_deg < MIN_SLOPE_DEG:",
     "    if False:",
     "test_gentle_terrain_yields_no_verdict"),

    ("indistinguishable geometry still gets a verdict",
     "candidate_check.py",
     "    if abs(ratio_motion - 1.0) < RATIO_MARGIN:",
     "    if False:",
     "test_no_verdict_when_the_two_hypotheses_predict_the_same_ratio"),

    ("elevation fit chases the outlier tails",
     "troposphere.py",
     '        if not robust:',
     '        if True:',
     "test_outliers_do_not_steer_the_slope"),

    ("trend removed in mm per metre, not mm per km",
     "troposphere.py",
     '    return disp_mm - fit["slope_mm_per_km"] * (elev_m / 1000.0)',
     '    return disp_mm - fit["slope_mm_per_km"] * elev_m',
     "test_removal_does_not_raise_the_core_scatter"),

    ("intercept subtracted along with the gradient",
     "troposphere.py",
     '    return disp_mm - fit["slope_mm_per_km"] * (elev_m / 1000.0)',
     '    return disp_mm - fit["slope_mm_per_km"] * (elev_m / 1000.0) - fit["intercept_mm"]',
     "test_intercept_is_not_subtracted"),

    ("extrapolated fit not flagged",
     "troposphere.py",
     '               leverage=(rng / iqr) if iqr > 0 else float("inf"),',
     '               leverage=1.0,',
     "test_leverage_is_flagged_when_the_fit_is_extrapolated"),

    ("flat terrain given a slope anyway",
     "troposphere.py",
     '    if np.ptp(e_all) <= 0:',
     '    if False:',
     "test_flat_terrain_is_refused"),
]


def failing_tests() -> set[str]:
    r = subprocess.run(
        [sys.executable, "-m", "pytest", str(ROOT / "tests"), "-q", "--no-header",
         "-p", "no:cacheprovider"],
        capture_output=True, text=True, cwd=ROOT)
    out = r.stdout + r.stderr
    names = set()
    for line in out.splitlines():
        if line.startswith("FAILED ") or line.startswith("ERROR "):
            part = line.split(" ", 1)[1].split(" ")[0]
            if "::" in part:
                names.add(part.split("::")[-1].split("[")[0])
    return names


def main() -> int:
    dirty = subprocess.run(["git", "status", "--porcelain", "src"],
                           capture_output=True, text=True, cwd=ROOT).stdout.strip()
    if dirty:
        print("src/ has uncommitted changes. Commit or stash first - this script\n"
              "edits source in place and must be able to restore it exactly.\n")
        print(dirty)
        return 2

    baseline = failing_tests()
    if baseline:
        print(f"Suite is not green to begin with: {sorted(baseline)}")
        return 2
    print("baseline: suite green\n")
    print(f"  {'MUTATION':<46}{'CAUGHT BY':<12}RESULT")
    print("  " + "-" * 74)

    escaped = []
    for label, fname, find, repl, expect in MUTATIONS:
        path = SRC / fname
        original = path.read_text(encoding="utf-8")
        try:
            if find not in original:
                print(f"  {label:<46}{'-':<12}TARGET GONE - update mutate.py")
                escaped.append(label)
                continue
            path.write_text(original.replace(find, repl, 1), encoding="utf-8")
            failed = failing_tests()
            if expect in failed:
                print(f"  {label:<46}{'yes':<12}caught ({len(failed)} failed)")
            else:
                print(f"  {label:<46}{'NO':<12}ESCAPED")
                escaped.append(label)
        finally:
            path.write_text(original, encoding="utf-8")
            for junk in SRC.rglob("__pycache__"):
                for f in junk.glob("*"):
                    f.unlink()

    print()
    if escaped:
        print(f"{len(escaped)} mutation(s) escaped - those tests do not protect "
              f"what they claim:")
        for e in escaped:
            print(f"  - {e}")
        return 1
    print(f"All {len(MUTATIONS)} mutations caught. The suite bites.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
