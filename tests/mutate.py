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

    ("duplicate products averaged instead of resolved",
     "timeseries.py",
     "out.append(best[1])",
     "best[1].value = sum(p.value for _, p in cands) / len(cands)\n        out.append(best[1])",
     "test_duplicate_products_are_resolved_not_averaged"),
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
