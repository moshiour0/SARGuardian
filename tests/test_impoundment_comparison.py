from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.compare_impoundment_geojson import compare_geojson


def _write(path: Path, features):
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}), encoding="utf-8")
    return path


def _feature(lon, lat, *, threshold=10, sill=3500, volume=2.5):
    return {
        "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [lon, lat]},
        "properties": {
            "name": "site",
            "threshold_height_m": threshold,
            "sill_m": sill,
            "max_volume_Mm3": volume,
            "by_height": {
                "10": {"volume_Mm3": volume if threshold <= 10 else 0.0},
                "25": {"volume_Mm3": volume * 2},
            },
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
