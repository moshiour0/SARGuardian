"""Regression tests for independent robust outlier handling in GOFF offsets."""

from __future__ import annotations

import numpy as np
import pytest

from goff_reader import mad_inlier_mask, robust_sigma


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
