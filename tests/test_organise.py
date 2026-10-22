from __future__ import annotations

import csv
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from organise import (bucket, is_already_filed, main, parse)


BETA = (
    "NISAR_L2_PR_GUNW_009_098_A_016_010_4000_SH_"
    "20260103T233921_20260103T233956_20260115T233922_20260115T233956_"
    "X05010_N_F_J_001.h5"
)
PROVISIONAL = (
    "NISAR_L2_PR_GUNW_024_098_A_016_025_4000_SH_"
    "20260702T233920_20260702T233955_20260714T233920_20260714T233954_"
    "P05023_F_F_J_001.h5"
)


def _put(root: Path, name: str, body: bytes = b"synthetic HDF5 placeholder") -> Path:
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(body)
    return path


def test_parse_records_crid_maturity_geometry_and_span(tmp_path):
    beta_path = _put(tmp_path, BETA)
    row = parse(beta_path)
    assert row is not None
    assert row["release_crid"] == "X05010"
    assert row["maturity"] == "BETA"
    assert row["product"] == "GUNW"
    assert row["direction"] == "A"
    assert row["path"] == 98
    assert row["frame"] == 16
    assert row["subswath"] == 10
    assert row["bandwidth"] == 4000
    assert row["span_days"] == 12

    provisional_path = _put(tmp_path, PROVISIONAL)
    prow = parse(provisional_path)
    assert prow is not None
    assert prow["release_crid"] == "P05023"
    assert prow["maturity"] == "PROVISIONAL"


def test_parse_rejects_invalid_dates_and_nonpositive_temporal_baselines(tmp_path):
    malformed = BETA.replace("20260103", "20261303", 1)
    same_date = BETA.replace("20260115", "20260103", 1)
    assert parse(_put(tmp_path, malformed)) is None
    assert parse(_put(tmp_path, same_date)) is None


def test_bucket_and_filed_path_classification_are_explicit(tmp_path):
    row = parse(_put(tmp_path, PROVISIONAL))
    assert row is not None
    assert bucket(row) == "2026-07_summer"
    dest = tmp_path / "data" / "nisar_l2"
    assert is_already_filed(dest / "GUNW" / "2026-07_summer" / PROVISIONAL, dest)
    assert is_already_filed(dest / "GOFF" / "co_event" / "file.h5", dest)
    assert not is_already_filed(dest / "_incoming" / PROVISIONAL, dest)
    assert not is_already_filed(dest / "misc" / "GUNW" / PROVISIONAL, dest)
    assert not is_already_filed(dest / "GUNW" / "misc" / PROVISIONAL, dest)


def test_incremental_organising_preserves_every_product_in_manifest(tmp_path, monkeypatch):
    source = tmp_path / "incoming"
    dest = tmp_path / "data" / "nisar_l2"
    _put(source, BETA, b"first synthetic granule")
    monkeypatch.setattr(sys, "argv", ["organise.py", "--src", str(source),
                                      "--dest", str(dest), "--apply"])
    assert main() == 0
    manifest = dest / "MANIFEST.csv"
    first = list(csv.DictReader(manifest.open(newline="", encoding="utf-8")))
    assert len(first) == 1
    assert first[0]["release_crid"] == "X05010"
    assert first[0]["maturity"] == "BETA"

    later = PROVISIONAL.replace("20260702", "20260726").replace(
        "20260714", "20260819"
    )
    # Keep a valid four-stamp name, but make this a separate pair/product.
    _put(source, later, b"second synthetic granule")
    monkeypatch.setattr(sys, "argv", ["organise.py", "--src", str(source),
                                      "--dest", str(dest), "--apply"])
    assert main() == 0
    final_rows = list(csv.DictReader(manifest.open(newline="", encoding="utf-8")))
    assert len(final_rows) == 2
    assert {row["release_crid"] for row in final_rows} == {"X05010", "P05023"}
    assert {row["storage_bucket"] for row in final_rows} == {"2026-01_winter", "2026-07_summer"}


def test_different_content_with_same_basename_fails_before_moving_anything(
    tmp_path, monkeypatch
):
    dest = tmp_path / "data" / "nisar_l2"
    target = _put(dest / "GUNW" / "2026-01_winter", BETA, b"committed bytes")
    source = tmp_path / "incoming"
    src = _put(source, BETA, b"modified bytes!")
    original_target = target.read_bytes()
    original_source = src.read_bytes()

    monkeypatch.setattr(sys, "argv", ["organise.py", "--src", str(source),
                                      "--dest", str(dest), "--apply"])
    assert main() == 1
    assert target.read_bytes() == original_target
    assert src.read_bytes() == original_source


def test_rerun_on_already_organised_folder_refreshes_full_manifest(tmp_path, monkeypatch):
    dest = tmp_path / "data" / "nisar_l2"
    _put(dest / "GUNW" / "2026-01_winter", BETA)
    monkeypatch.setattr(sys, "argv", ["organise.py", "--src", str(dest),
                                      "--dest", str(dest), "--apply"])
    assert main() == 0
    rows = list(csv.DictReader((dest / "MANIFEST.csv").open(newline="", encoding="utf-8")))
    assert len(rows) == 1
    assert rows[0]["relative_path"] == f"GUNW/2026-01_winter/{BETA}"
