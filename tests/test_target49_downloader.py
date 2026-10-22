from __future__ import annotations

import sys

ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.download_target49_products import direct_product_url, download_exact, main, product_info


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



def test_direct_product_url_selects_the_exact_maturity_collection():
    beta = "NISAR_L2_PR_GUNW_009_098_A_016_010_4000_SH_20260103T233921_20260103T233956_20260115T233922_20260115T233956_X05010_N_F_J_001.h5"
    provisional = "NISAR_L2_PR_GUNW_024_098_A_016_025_4000_SH_20260702T233920_20260702T233955_20260714T233920_20260714T233954_P05023_F_F_J_001.h5"
    assert product_info("GUNW-006-BETA-control") == ("GUNW", "BETA")
    assert product_info("GUNW-009-PROVISIONAL") == ("GUNW", "PROVISIONAL")
    assert direct_product_url(beta, "GUNW", "BETA") == (
        "https://nisar.asf.earthdatacloud.nasa.gov/NISAR/"
        "NISAR_L2_GUNW_BETA_V1/" + beta[:-3] + "/" + beta
    )
    assert direct_product_url(provisional, "GUNW", "PROVISIONAL").split("/")[5] == "NISAR_L2_GUNW_PROVISIONAL_V1"


def test_direct_download_streams_and_hashes_without_storing_url(tmp_path):
    import hashlib

    class FakeResponse:
        headers = {"Content-Length": "5"}
        def raise_for_status(self):
            return None
        def iter_content(self, chunk_size):
            assert chunk_size == 8 * 1024 * 1024
            yield b"he"
            yield b"llo"
        def close(self):
            self.closed = True

    class FakeSession:
        def __init__(self):
            self.response = FakeResponse()
            self.called_url = None
        def get(self, url, **kwargs):
            self.called_url = url
            assert kwargs["stream"] is True
            assert kwargs["allow_redirects"] is True
            assert kwargs["timeout"] == (30, 120)
            return self.response

    session = FakeSession()
    dest = tmp_path / "small-test.h5"
    nbytes, digest = download_exact(
        session,
        "https://nisar.asf.earthdatacloud.nasa.gov/NISAR/NISAR_L2_GUNW_BETA_V1/x/x.h5",
        dest,
        retries=1,
    )
    assert nbytes == 5
    assert dest.read_bytes() == b"hello"
    assert digest == hashlib.sha256(b"hello").hexdigest()
    # The direct signed URL is passed to the session, never written to the file or manifest.
    assert dest.read_bytes() != session.called_url.encode()


def test_direct_download_rejects_noncanonical_initial_host(tmp_path):
    import pytest
    from tools.download_target49_products import download_exact

    class UnusedSession:
        def get(self, *args, **kwargs):
            raise AssertionError("must reject URL before making a request")

    with pytest.raises(ValueError, match="canonical"):
        download_exact(UnusedSession(), "https://example.com/a.h5", tmp_path / "a.h5", retries=1)
