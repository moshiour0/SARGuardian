from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from paths import EXCLUDE_DIRS, _is_excluded, find_products, resolve


GUNW = "NISAR_L2_PR_GUNW_009_098_A_016_010_4000_SH_20260103T233921_20260103T233956_20260115T233922_20260115T233956_X05010_N_F_J_001.h5"
GOFF = "NISAR_L2_PR_GOFF_009_098_A_016_010_4000_SH_20260103T233921_20260103T233956_20260115T233922_20260115T233956_X05010_N_F_J_001.h5"


def _touch(root: Path, relative: str) -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"inventory fixture")
    return path


def test_excluded_directory_names_are_never_searched(tmp_path):
    excluded = _touch(tmp_path, f".venv/site-packages/{GUNW}")
    assert _is_excluded(excluded)
    assert _is_excluded(tmp_path / "node_modules" / GUNW)
    assert _is_excluded(tmp_path / ".git" / "objects" / GUNW)
    assert set(("site-packages", "node_modules", ".git")) <= EXCLUDE_DIRS


def test_find_products_filters_by_product_and_deduplicates_overlapping_roots(tmp_path):
    g = _touch(tmp_path, f"data/nisar_l2/{GUNW}")
    o = _touch(tmp_path, f"data/nisar_l2/_incoming/nested/{GOFF}")
    _touch(tmp_path, "data/nisar_l2/_incoming/not-a-product.h5")
    _touch(tmp_path, f"data/nisar_l2/.venv/{GUNW}")

    # The same file appears through both roots; its canonical path must occur once.
    found = find_products(roots=[tmp_path, tmp_path / "data" / "nisar_l2"])
    assert found.count(g.resolve()) == 1
    assert found.count(o.resolve()) == 1
    assert len(found) == 2

    only_goff = find_products("GOFF", roots=[tmp_path])
    assert only_goff == [o.resolve()]
    only_gunw = find_products("GUNW", roots=[tmp_path])
    assert only_gunw == [g.resolve()]


def test_strict_filename_filter_is_explicit(tmp_path):
    nisar = _touch(tmp_path, GUNW)
    not_nisar = _touch(tmp_path, "random_fixture.h5")
    assert find_products(roots=[tmp_path], strict=True) == [nisar.resolve()]
    assert set(find_products(roots=[tmp_path], strict=False)) == {
        nisar.resolve(), not_nisar.resolve()
    }


def test_resolve_honours_absolute_paths_and_prefers_existing_cwd_relative_paths(
    tmp_path, monkeypatch
):
    absolute = _touch(tmp_path, "absolute.txt")
    assert resolve(absolute) == absolute

    monkeypatch.chdir(tmp_path)
    local = _touch(tmp_path, "local/relative.txt")
    assert resolve(Path("local/relative.txt")) == local.resolve()
