from __future__ import annotations

import sys

ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.download_target49_products import main


def test_single_label_dry_run_selects_exact_requested_granule(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(sys, "argv", [
        "download_target49_products.py",
        "--label", "GUNW-009-PROVISIONAL",
        "--outdir", str(tmp_path),
        "--dry-run",
    ])
    assert main() == 0
    out = capsys.readouterr().out
    assert "products=1" in out
    assert "GUNW-009-PROVISIONAL" in out
    assert "GUNW-007-PROVISIONAL" not in out
    assert "NISAR_L2_PR_GUNW_024_098_A_016_025" in out
    assert "DRY_RUN: no network calls" in out


def test_multiple_labels_dry_run_does_not_need_credentials(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(sys, "argv", [
        "download_target49_products.py",
        "--label", "GUNW-009-PROVISIONAL",
        "--label", "GOFF-027-PROVISIONAL",
        "--outdir", str(tmp_path),
        "--dry-run",
    ])
    assert main() == 0
    out = capsys.readouterr().out
    assert "products=2" in out
    assert "GUNW-009-PROVISIONAL" in out
    assert "GOFF-027-PROVISIONAL" in out
    assert "GUNW-007-PROVISIONAL" not in out
