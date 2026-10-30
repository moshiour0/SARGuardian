"""Regression tests for independent robust outlier handling in GOFF offsets."""

from __future__ import annotations

import h5py
import numpy as np
import pytest

import synth

import goff_reader as goff
from goff_reader import export_clipped, layer_groups, mad_inlier_mask, read_goff, robust_sigma


def test_layer_groups_ignore_processing_metadata_parameter_nodes():
    datasets = {
        "science/LSAR/GOFF/grids/frequencyA/pixelOffsets/HH/layer2/slantRangeOffset": (
            (11, 11), np.dtype("float32")
        ),
        "science/LSAR/GOFF/metadata/processingInformation/parameters/"
        "pixelOffsets/frequencyA/layer2/algorithmVersion": ((), np.dtype("int32")),
    }

    assert layer_groups(datasets) == [
        "science/LSAR/GOFF/grids/frequencyA/pixelOffsets/HH/layer2"
    ]

def test_azimuth_only_outlier_is_rejected_without_borrowing_range_statistics():
    rng = np.random.default_rng(20261022)
    range_m = rng.normal(0.0, 0.01, size=(20, 20))
    azimuth_m = rng.normal(0.0, 0.01, size=(20, 20))
    target = (3, 7)
    range_m[target] = 0.0
    azimuth_m[target] = 12.0
    valid = np.ones(range_m.shape, dtype=bool)

    range_keep = mad_inlier_mask(range_m, valid, 4.0)
    azimuth_keep = mad_inlier_mask(azimuth_m, valid, 4.0)
    joint_keep = range_keep & azimuth_keep

    assert robust_sigma(range_m[valid]) > 0
    assert range_keep[target], "range component itself is not an outlier"
    assert not azimuth_keep[target], "azimuth's own MAD must identify the injected outlier"
    assert not joint_keep[target], "the 2-D offset vector must not retain an axis outlier"


def test_range_only_outlier_is_rejected_even_when_azimuth_is_valid():
    rng = np.random.default_rng(20261023)
    range_m = rng.normal(0.0, 0.02, size=(18, 18))
    azimuth_m = rng.normal(0.0, 0.02, size=(18, 18))
    target = (9, 11)
    range_m[target] = -15.0
    azimuth_m[target] = 0.0
    valid = np.ones(range_m.shape, dtype=bool)

    range_keep = mad_inlier_mask(range_m, valid, 4.0)
    azimuth_keep = mad_inlier_mask(azimuth_m, valid, 4.0)

    assert not range_keep[target]
    assert azimuth_keep[target]
    assert not (range_keep & azimuth_keep)[target]


def test_mad_outlier_mask_preserves_upstream_invalid_cells():
    values = np.linspace(-1.0, 1.0, 121).reshape(11, 11)
    base = np.ones(values.shape, dtype=bool)
    base[2, 3] = False
    values[4, 5] = 1e6

    keep = mad_inlier_mask(values, base, 3.0)

    assert not keep[2, 3], "the outlier gate must not re-enable low-quality pixels"
    assert not keep[4, 5], "a large robust outlier must be rejected"


def test_disabled_or_under_supported_culling_only_preserves_finite_base_pixels():
    values = np.array([[0.0, 1.0, 1000.0]], dtype=float)
    base = np.array([[True, True, False]])
    assert np.array_equal(mad_inlier_mask(values, base, 0.0), base)
    assert np.array_equal(mad_inlier_mask(values, base, 4.0), base)


def test_mad_outlier_mask_rejects_shape_mismatch():
    with pytest.raises(ValueError, match="identical shapes"):
        mad_inlier_mask(np.zeros((3, 3)), np.ones((3, 2), dtype=bool), 4.0)


def _write_minimal_goff(path, *, missing_snr=False):
    shape = (11, 11)
    base = "science/LSAR/GOFF/grids/frequencyA/pixelOffsets/HH/layer2"
    with h5py.File(path, "w") as h5:
        h5.create_dataset(f"{base}/slantRangeOffset", data=np.zeros(shape))
        h5.create_dataset(f"{base}/alongTrackOffset", data=np.zeros(shape))
        h5.create_dataset(f"{base}/correlationSurfacePeak", data=np.ones(shape))
        if not missing_snr:
            h5.create_dataset(f"{base}/snr", data=np.ones(shape) * 5.0)
        h5.create_dataset(f"{base}/xCoordinates", data=np.arange(shape[1]) + 500000.0)
        h5.create_dataset(f"{base}/yCoordinates", data=np.arange(shape[0]) + 3200000.0)
        projection = h5.create_dataset(f"{base}/projection", data=np.array([32645]))
        projection.attrs["epsg_code"] = 32645
    return path


def test_read_goff_uses_active_aoi_ring_for_measurement_mask(tmp_path, monkeypatch):
    shape = (11, 11)
    base = "science/LSAR/GOFF/grids/frequencyA/pixelOffsets/HH/layer2"
    path = tmp_path / "complete_geographic_goff.h5"
    with h5py.File(path, "w") as h5:
        h5.create_dataset(f"{base}/slantRangeOffset", data=np.zeros(shape))
        h5.create_dataset(f"{base}/alongTrackOffset", data=np.zeros(shape))
        h5.create_dataset(f"{base}/correlationSurfacePeak", data=np.ones(shape))
        h5.create_dataset(f"{base}/snr", data=np.ones(shape) * 5.0)
        h5.create_dataset(f"{base}/xCoordinates", data=np.linspace(84.8, 85.3, shape[1]))
        h5.create_dataset(f"{base}/yCoordinates", data=np.linspace(28.1, 28.6, shape[0]))
        projection = h5.create_dataset(f"{base}/projection", data=np.array([4326]))
        projection.attrs["epsg_code"] = 4326

    ring = [(84.9, 28.2), (85.2, 28.2), (85.2, 28.5), (84.9, 28.5)]
    monkeypatch.setitem(goff._g.AOIS, "regression", ring)
    monkeypatch.setattr(goff._g, "AOI_NAME", "regression")
    monkeypatch.setattr(goff._g, "AOI_RING", ring)

    result = read_goff(path, layer="layer2", clip_aoi=True, auto_ref=False, deramp=False)
    layer = result["layers"]["HH/layer2"]

    assert result["aoi_name"] == "regression"
    assert layer["aoi_mask_applied"] is True
    assert layer["aoi_px"] > 0
    assert int(layer["valid"].sum()) > 0

def test_read_goff_refuses_to_silently_skip_a_missing_quality_gate(tmp_path):
    path = _write_minimal_goff(tmp_path / "missing_snr.h5", missing_snr=True)
    with pytest.raises(RuntimeError, match="missing required GOFF measurement/quality"):
        read_goff(path, layer="layer2", deramp=False)


def test_read_goff_refuses_an_unknown_explicit_layer(tmp_path):
    path = _write_minimal_goff(tmp_path / "layer2_only.h5")
    with pytest.raises(RuntimeError, match="Requested GOFF layer 'layer3' not found"):
        read_goff(path, layer="layer3", deramp=False)


def test_coherent_high_offset_patch_survives_for_spatial_review(tmp_path):
    """A real compact offset anomaly may be a global MAD outlier, but remains map-reviewable."""
    import rasterio

    path = tmp_path / synth.granule_name(product="GOFF")
    xs, ys = synth.make_grid(120)
    shape = (len(ys), len(xs))
    rng = np.random.default_rng(20261029)
    range_m = rng.normal(0.0, 0.01, size=shape)
    azimuth_m = rng.normal(0.0, 0.015, size=shape)
    corr = np.full(shape, 0.9, dtype=np.float32)
    snr = np.full(shape, 10.0, dtype=np.float32)
    range_variance = np.full(shape, 0.0004, dtype=np.float32)
    azimuth_variance = np.full(shape, 0.0009, dtype=np.float32)
    cross_offset_variance = np.full(shape, 0.0001, dtype=np.float32)
    lon, lat = synth.LHENDE_INSIDE
    col = int(np.argmin(np.abs(xs - lon)))
    row = int(np.argmin(np.abs(ys - lat)))
    patch = np.zeros(shape, dtype=bool)
    patch[row - 1:row + 2, col - 1:col + 2] = True
    range_m[patch] += 0.8
    azimuth_m[patch] += 0.5

    base = f"{synth.GOFF_GRID}/HH/layer2"
    with h5py.File(path, "w") as h5:
        group = h5.create_group(base)
        group.create_dataset("slantRangeOffset", data=range_m)
        group.create_dataset("alongTrackOffset", data=azimuth_m)
        group.create_dataset("correlationSurfacePeak", data=corr)
        group.create_dataset("snr", data=snr)
        group.create_dataset("slantRangeOffsetVariance", data=range_variance)
        group.create_dataset("alongTrackOffsetVariance", data=azimuth_variance)
        group.create_dataset("crossOffsetVariance", data=cross_offset_variance)
        group.create_dataset("xCoordinates", data=xs)
        group.create_dataset("yCoordinates", data=ys)
        projection = group.create_dataset("projection", data=np.uint32(4326))
        projection.attrs["epsg_code"] = 4326

    result = read_goff(
        path, layer="layer2", clip_aoi=True, auto_ref=False, deramp=False, cull_sigma=4.0
    )
    layer = result["layers"]["HH/layer2"]
    assert layer["screen_valid"][row, col]
    assert layer["range_var"][row, col] == range_variance[row, col]
    assert layer["azimuth_var"][row, col] == azimuth_variance[row, col]
    assert layer["cross_offset_variance"][row, col] == cross_offset_variance[row, col]
    assert layer["outlier_mask"][row, col], "the truth patch must exercise the global MAD gate"
    assert not layer["valid"][row, col], "conservative pairwise stats still exclude the outlier"
    rejected = sum(int(n) for n in layer["gate_rejections"].values())
    assert layer["source_total_pixels"] - int(layer["valid"].sum()) == rejected
    assert layer["gates"]["strict robust-inlier pixels"] == int(layer["valid"].sum())

    exported = export_clipped(result, tmp_path / "exports")
    assert exported and exported[0]["export_status"] == "WRITTEN"
    raster_path = tmp_path / "exports" / exported[0]["export_file"]
    with rasterio.open(raster_path) as src:
        assert src.count == 7
        assert src.descriptions[4] == "slant-range offset variance (mm^2)"
        assert src.descriptions[5] == "along-track offset variance (mm^2)"
        assert src.descriptions[6].startswith("cross offset variance (mm^2")
        rr, cc = src.index(lon, lat)
        assert np.isclose(src.read(5)[rr, cc], 400.0)
        assert np.isclose(src.read(6)[rr, cc], 900.0)
        assert np.isclose(src.read(7)[rr, cc], 100.0)
        assert np.isfinite(src.read(1)[rr, cc]), "screen raster must retain the high range offset"
        assert np.isfinite(src.read(2)[rr, cc]), "screen raster must retain the paired azimuth offset"
        assert src.read(4)[rr, cc] == 0.0, "band 4 must identify a global-MAD outlier"
        import json
        tags = src.tags()
        assert tags["outliers_preserved_for_screening"] == "true"
        assert int(tags["source_total_pixels"]) == shape[0] * shape[1]
        gate_counts = json.loads(tags["gate_counts_json"])
        gate_rejections = json.loads(tags["gate_rejections_json"])
        assert gate_counts["strict robust-inlier pixels"] == int(layer["valid"].sum())
        assert sum(int(n) for n in gate_rejections.values()) == shape[0] * shape[1] - int(layer["valid"].sum())
