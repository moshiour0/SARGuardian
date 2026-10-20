from __future__ import annotations

from urllib.parse import parse_qs, urlparse

import nisar_acquisition as acquisition


def test_asf_search_uses_a_documented_json_output_format(monkeypatch):
    captured = {}

    def fake_get(url):
        parsed = urlparse(url)
        captured.update(parse_qs(parsed.query))
        return {"results": [{"granuleName": "NISAR_L2_PR_GUNW_example"}]}

    monkeypatch.setattr(acquisition, "_get", fake_get)
    rows = acquisition.asf_query(
        warn_on_cap=False,
        platform="NISAR",
        processingLevel="GUNW",
    )

    assert captured["output"] == ["json"]
    assert captured["platform"] == ["NISAR"]
    assert captured["processingLevel"] == ["GUNW"]
    assert rows == [{"granuleName": "NISAR_L2_PR_GUNW_example"}]
    assert "jsonlite" not in captured["output"]
