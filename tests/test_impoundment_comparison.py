from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.compare_impoundment_geojson import compare_geojson, main as compare_main


def _write(path: Path, features):
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}), encoding="utf-8")
    return path


def _feature(lon, lat, *, threshold=10, sill=3500, volume=2.5):
    by_height = {}
    for height in (10, 25, 50, 100, 150):
        if height < threshold:
            continue
        by_height[str(float(height))] = {
            "volume_Mm3": volume * height / 10,
            "area_km2": 0.1 * height / 10,
            "backwater_km": 0.25 * height / 10,
            "truncated": False,
        }
    return {
        "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [lon, lat]},
        "properties": {
            "name": "site",
            "threshold_height_m": threshold,
            "sill_m": sill,
            "max_volume_Mm3": volume,
            "by_height": by_height,
        },
    }


def test_exactly_equal_geojson_is_reported_as_semantic_match(tmp_path):
    baseline = _write(tmp_path / "baseline.geojson", [_feature(85.5, 28.3), _feature(85.51, 28.31, threshold=25)])
    candidate = _write(tmp_path / "candidate.geojson", [_feature(85.5, 28.3), _feature(85.51, 28.31, threshold=25)])
    summary, changes = compare_geojson(baseline, candidate)
    assert summary["status"] == "EXACT_SEMANTIC_MATCH"
    assert summary["exact_coordinate_matches"] == 2
    assert summary["response_counts_by_threshold_height_m"]["baseline"]["10"] == 1
    assert summary["response_counts_by_threshold_height_m"]["candidate"]["10"] == 1
    assert changes == []


def test_shifted_site_is_not_falsely_paired_to_its_nearest_neighbour(tmp_path):
    baseline = _write(tmp_path / "baseline.geojson", [_feature(85.5, 28.3)])
    candidate = _write(tmp_path / "candidate.geojson", [_feature(85.5001, 28.3001)])
    summary, changes = compare_geojson(baseline, candidate)
    assert summary["status"] == "REVIEW_REQUIRED"
    assert summary["exact_coordinate_matches"] == 0
    assert summary["baseline_sites_without_exact_match"] == 1
    assert summary["candidate_sites_without_exact_match"] == 1
    assert {row["status"] for row in changes} == {"BASELINE_SITE_MISSING", "NEW_CANDIDATE_SITE"}


def test_property_mismatch_at_same_exact_site_is_reported(tmp_path):
    baseline = _write(tmp_path / "baseline.geojson", [_feature(85.5, 28.3, sill=3500)])
    candidate = _write(tmp_path / "candidate.geojson", [_feature(85.5, 28.3, sill=3525)])
    summary, changes = compare_geojson(baseline, candidate)
    assert summary["status"] == "REVIEW_REQUIRED"
    assert summary["exact_coordinate_matches"] == 1
    assert summary["property_change_count_at_exactly_matched_sites"] == 1
    assert any(row["field"] == "sill_m" for row in changes)


def test_small_float_rounding_is_tolerated_but_not_coordinate_changes(tmp_path):
    baseline = _write(tmp_path / "baseline.geojson", [_feature(85.5, 28.3, volume=2.50000000001)])
    candidate = _write(tmp_path / "candidate.geojson", [_feature(85.5, 28.3, volume=2.5)])
    summary, changes = compare_geojson(baseline, candidate, tolerance=1e-9)
    assert summary["status"] == "EXACT_SEMANTIC_MATCH"
    assert changes == []


def test_duplicate_exact_coordinates_are_rejected(tmp_path):
    duplicated = _write(tmp_path / "duplicate.geojson", [
        _feature(85.5, 28.3), _feature(85.5, 28.3)
    ])
    with pytest.raises(ValueError, match="duplicate exact site coordinate"):
        compare_geojson(duplicated, duplicated)


def test_response_counts_use_by_height_results_not_new_top_level_threshold():
    """Historical GeoJSONs lack threshold_height_m; compare using by_height."""
    from tools.compare_impoundment_geojson import response_counts

    legacy_site = {
        "geometry": {"type": "Point", "coordinates": [85.5, 28.3]},
        "properties": {
            "by_height": {
                "50.0": {"volume_Mm3": 5.0, "area_km2": 0.2, "backwater_km": 0.4},
                "100": {"volume_Mm3": 10.0, "area_km2": 0.3, "backwater_km": 0.6},
            },
        },
    }
    new_schema_site = {
        "geometry": {"type": "Point", "coordinates": [85.6, 28.4]},
        "properties": {
            "threshold_height_m": 10.0,
            "by_height": {
                "10.0": {"volume_Mm3": 1.0, "area_km2": 0.1, "backwater_km": 0.2},
                "25.0": {"volume_Mm3": 2.0, "area_km2": 0.2, "backwater_km": 0.3},
            },
        },
    }
    counts = response_counts([legacy_site, new_schema_site], [10, 25, 50, 100])
    assert counts == {"10": 1, "25": 1, "50": 1, "100": 1}


def test_response_count_skips_partial_or_nonfinite_height_results():
    from tools.compare_impoundment_geojson import response_counts

    rows = [
        {"properties": {"threshold_height_m": 10.0,
                        "by_height": {"10": {"volume_Mm3": 1.0}}}},
        {"properties": {"threshold_height_m": 10.0,
                        "by_height": {"10": {
                            "volume_Mm3": float("nan"), "area_km2": 0.1, "backwater_km": 0.2
                        }}}},
    ]
    assert response_counts(rows, [10]) == {"10": 0}


def test_nearest_neighbour_is_review_only_and_never_an_exact_site_match(tmp_path):
    baseline = _write(tmp_path / "baseline.geojson", [_feature(85.5, 28.3, sill=3500)])
    candidate = _write(tmp_path / "candidate.geojson", [_feature(85.5001, 28.3001, sill=3525)])
    summary, changes = compare_geojson(baseline, candidate)

    assert summary["exact_coordinate_matches"] == 0
    assert summary["status"] == "REVIEW_REQUIRED"
    review = summary["nearest_neighbour_review_only"]
    assert review["baseline_unmatched_to_candidate"]["count"] == 1
    assert review["baseline_unmatched_to_candidate"]["max_m"] < 100
    assert review["baseline_unmatched_to_candidate"]["within_100m"] == 1
    missing = next(row for row in changes if row["status"] == "BASELINE_SITE_MISSING")
    assert missing["nearest_neighbour_is_identity_match"] is False
    assert missing["nearest_candidate_distance_m"] < 100
    assert missing["nearest_candidate_sill_m"] == 3525
    assert missing["nearest_candidate_sill_m_minus_source"] == 25


def test_truncated_height_results_are_reported_separately_as_lower_bounds():
    from tools.compare_impoundment_geojson import response_counts, truncated_response_counts

    rows = [
        {"properties": {"by_height": {"10.0": {
            "volume_Mm3": 1.0, "area_km2": 0.1,
            "backwater_km": 0.2, "truncated": True,
        }}}},
        {"properties": {"by_height": {"10.0": {
            "volume_Mm3": 2.0, "area_km2": 0.2,
            "backwater_km": 0.4, "truncated": False,
        }}}},
        {"properties": {"by_height": {"10.0": {
            "volume_Mm3": 3.0, "area_km2": 0.3,
            "backwater_km": 0.6, "truncated": "true",
        }}}},
    ]
    assert response_counts(rows, [10]) == {"10": 3}
    assert truncated_response_counts(rows, [10]) == {"10": 1}


def test_coordinate_rounding_is_separate_from_exact_identity(tmp_path):
    baseline = _write(tmp_path / "baseline.geojson", [_feature(85.500005, 28.3)])
    candidate = _write(tmp_path / "candidate.geojson", [_feature(85.5, 28.3)])
    summary, changes = compare_geojson(baseline, candidate)
    assert summary["exact_coordinate_matches"] == 0
    assert summary["coordinate_rounding_matches"] == 1
    assert summary["baseline_sites_unmatched_after_coordinate_tolerance"] == 0
    assert summary["candidate_sites_unmatched_after_coordinate_tolerance"] == 0
    assert summary["status"] == "SEMANTIC_MATCH_WITH_COORDINATE_ROUNDING"
    assert summary["review_row_count"] == 0
    assert changes[0]["status"] == "COORDINATE_ROUNDING_MATCH"
    assert changes[0]["nearest_neighbour_is_identity_match"] is False


def test_property_change_at_rounded_coordinate_pair_requires_review(tmp_path):
    baseline = _write(tmp_path / "baseline.geojson", [_feature(85.500005, 28.3, sill=3500)])
    candidate = _write(tmp_path / "candidate.geojson", [_feature(85.5, 28.3, sill=3525)])
    summary, changes = compare_geojson(baseline, candidate)
    assert summary["coordinate_rounding_matches"] == 1
    assert summary["property_change_count_at_rounded_coordinate_matches"] == 1
    assert summary["status"] == "REVIEW_REQUIRED"
    assert any(row["status"] == "PROPERTY_CHANGED_AT_ROUNDED_COORDINATE_MATCH"
               and row["field"] == "sill_m" for row in changes)


def test_cli_csv_keeps_nearest_neighbour_review_columns(tmp_path, monkeypatch):
    import csv
    import sys

    baseline = _write(tmp_path / "baseline.geojson", [_feature(85.5, 28.3)])
    candidate = _write(tmp_path / "candidate.geojson", [_feature(85.5001, 28.3001)])
    out_json, out_csv = tmp_path / "compare.json", tmp_path / "compare.csv"
    monkeypatch.setattr(sys, "argv", [
        "compare_impoundment_geojson.py",
        "--baseline", str(baseline),
        "--candidate", str(candidate),
        "--json", str(out_json),
        "--csv", str(out_csv),
    ])
    assert compare_main() == 0
    with out_csv.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
        assert "nearest_candidate_distance_m" in rows[0]
        assert "nearest_baseline_distance_m" in rows[0]
        assert "nearest_candidate_sill_m" in rows[0]
        assert "nearest_candidate_sill_m_minus_source" in rows[0]
        assert "nearest_neighbour_is_identity_match" in rows[0]
    summary = json.loads(out_json.read_text(encoding="utf-8"))
    assert summary["status"] == "REVIEW_REQUIRED"
    assert summary["exact_coordinate_matches"] == 0
