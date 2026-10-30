"""Regression tests for spatial-grid and per-pixel support invariants."""

from __future__ import annotations

from datetime import date
from pathlib import Path

import numpy as np
import pytest
import rasterio
from affine import Affine

from pixel_stack import align, stack, write_geotiff


def _write_raster(path: Path, values: np.ndarray, transform: Affine) -> None:
    with rasterio.open(
        path, "w", driver="GTiff", width=values.shape[1], height=values.shape[0],
        count=2, dtype="float32", crs="EPSG:32645", transform=transform,
        nodata=np.nan,
    ) as dst:
        dst.write(values.astype("float32"), 1)
        dst.write(np.full(values.shape, 0.8, dtype="float32"), 2)


def test_align_refuses_half_pixel_origin_instead_of_rounding_it(tmp_path):
    values = np.ones((4, 4), dtype="float32")
    a = tmp_path / "a.tif"
    b = tmp_path / "b.tif"
    _write_raster(a, values, Affine(80, 0, 500000, 0, -80, 3200000))
    _write_raster(b, values, Affine(80, 0, 500040, 0, -80, 3200000))
    rasters = [{"path": a}, {"path": b}]
    with pytest.raises(SystemExit, match="not aligned to the common pixel lattice"):
        align(rasters)


def test_stack_refuses_one_pair_because_velocity_uncertainty_has_no_dof():
    with pytest.raises(ValueError, match="at least two consecutive pair intervals"):
        stack([{"path": Path("unused"), "ref": date(2026, 1, 1), "sec": date(2026, 1, 13)}])


def test_missing_pair_pixels_are_reported_as_partial_not_filled_with_zero(tmp_path):
    base_transform = Affine(80, 0, 500000, 0, -80, 3200000)
    first = np.full((3, 3), 10.0, dtype="float32")
    second = np.full((3, 3), 10.0, dtype="float32")
    first[0, 0] = np.nan
    second[1, 1] = np.nan
    p1 = tmp_path / "20260101_20260113.tif"
    p2 = tmp_path / "20260113_20260125.tif"
    _write_raster(p1, first, base_transform)
    _write_raster(p2, second, base_transform)
    run = [
        {"path": p1, "ref": date(2026, 1, 1), "sec": date(2026, 1, 13), "geom": "A098"},
        {"path": p2, "ref": date(2026, 1, 13), "sec": date(2026, 1, 25), "geom": "A098"},
    ]
    result = stack(run, remove_common_mode=False)
    assert result["observed_pair_count"][0, 0] == 1
    assert result["observed_pair_count"][1, 1] == 1
    assert not result["valid"][0, 0]
    assert not result["valid"][1, 1]
    assert np.isnan(result["cumulative"][1, 0, 0])
    assert np.isnan(result["cumulative"][2, 1, 1])
    assert result["complete_series_pixels"] == 7
    assert result["total_grid_pixels"] == 9


def test_partial_support_counts_are_exported_separately_from_velocity(tmp_path):
    base_transform = Affine(80, 0, 500000, 0, -80, 3200000)
    first = np.full((3, 3), 10.0, dtype="float32")
    second = np.full((3, 3), 10.0, dtype="float32")
    first[0, 0] = np.nan
    second[1, 1] = np.nan
    p1 = tmp_path / "20260101_20260113.tif"
    p2 = tmp_path / "20260113_20260125.tif"
    _write_raster(p1, first, base_transform)
    _write_raster(p2, second, base_transform)
    run = [
        {"path": p1, "ref": date(2026, 1, 1), "sec": date(2026, 1, 13), "geom": "A098"},
        {"path": p2, "ref": date(2026, 1, 13), "sec": date(2026, 1, 25), "geom": "A098"},
    ]
    result = stack(run, remove_common_mode=False)
    assert result["observed_pair_count"][0, 0] == 1
    assert result["longest_contiguous_pair_count"][0, 0] == 1
    assert result["partial_support"][0, 0]
    assert not result["valid"][0, 0]

    out = tmp_path / "stack.tif"
    write_geotiff(result, out)
    with rasterio.open(out) as src:
        assert src.count == 4
        assert "complete connected chain only" in src.descriptions[0]
        assert "Observed valid pair count" in src.descriptions[2]
        assert "Longest contiguous valid pair count" in src.descriptions[3]
        assert src.read(3)[0, 0] == 1.0
        assert src.read(4)[0, 0] == 1.0
        assert src.read(1, masked=True).mask[0, 0]
        assert int(src.tags()["partial_support_pixels"]) == 2
        assert src.tags()["support_policy"] == "velocity bands require complete connected pair chain"
        assert src.tags()["common_mode_policy"] == "PAIRWISE_MEDIAN_RETAINED"
        assert src.tags()["common_mode_removed"] == "false"
        assert int(src.tags()["source_pair_count"]) == 2
        import json
        pairs = json.loads(src.tags()["pair_records_json"])
        assert [p["reference_date"] for p in pairs] == ["2026-01-01", "2026-01-13"]


def test_common_mode_removal_is_explicit_and_records_each_removed_pair_median(tmp_path):
    base_transform = Affine(80, 0, 500000, 0, -80, 3200000)
    first = np.full((3, 3), 10.0, dtype="float32")
    second = np.full((3, 3), 12.0, dtype="float32")
    p1 = tmp_path / "20260101_20260113.tif"
    p2 = tmp_path / "20260113_20260125.tif"
    _write_raster(p1, first, base_transform)
    _write_raster(p2, second, base_transform)
    run = [
        {"path": p1, "ref": date(2026, 1, 1), "sec": date(2026, 1, 13), "geom": "A098"},
        {"path": p2, "ref": date(2026, 1, 13), "sec": date(2026, 1, 25), "geom": "A098"},
    ]
    removed = stack(run, remove_common_mode=True)
    retained = stack(run, remove_common_mode=False)
    assert removed["common_mode_policy"] == "PAIRWISE_MEDIAN_REMOVED"
    assert removed["common_mode_pair_medians_mm"] == [10.0, 12.0]
    assert retained["common_mode_policy"] == "PAIRWISE_MEDIAN_RETAINED"
    assert retained["common_mode_pair_medians_mm"] == [10.0, 12.0]
    # These two choices cannot be silently conflated; broad common motion is
    # removed in one branch and retained in the other diagnostic branch.
    assert np.allclose(removed["cumulative"][:, 1, 1], [0.0, 0.0, 0.0])
    assert np.allclose(retained["cumulative"][:, 1, 1], [0.0, 10.0, 22.0])
