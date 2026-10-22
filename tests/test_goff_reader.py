"""Regression tests for independent robust outlier handling in GOFF offsets."""

from __future__ import annotations

import h5py
import numpy as np
import pytest

import goff_reader as goff
from goff_reader import layer_groups, mad_inlier_mask, read_goff, robust_sigma


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
