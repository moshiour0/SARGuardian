"""
Exports must land on one grid, or they cannot be differenced.

The fixture uses an explicit WGS84 CRS and cell-centre coordinates aligned to
the common AOI lattice. Missing CRS must be tested separately, not implied by
passing epsg=None for projected coordinates.
"""
from __future__ import annotations

import numpy as np
import pytest

from gunw_reader import aoi_grid, place_on_grid

RES = 0.0008
RING = [(85.40, 28.20), (85.80, 28.20), (85.80, 28.60), (85.40, 28.60)]


def _centres(x_start, y_start, nx, ny):
    return x_start + RES * np.arange(nx), y_start - RES * np.arange(ny)


def test_two_products_covering_different_ground_get_the_same_grid():
    a = aoi_grid(*_centres(85.5588, 28.4812, 3, 3), RING, 4326)
    b = aoi_grid(*_centres(85.5604, 28.4804, 4, 4), RING, 4326)
    assert a is not None and b is not None
    assert (a["height"], a["width"]) == (b["height"], b["width"])
    assert a["transform"] == b["transform"]


def test_the_same_map_coordinate_lands_in_the_same_cell():
    xs_a, ys_a = _centres(85.5588, 28.4812, 3, 3)
    xs_b, ys_b = _centres(85.5604, 28.4804, 4, 4)
    ga = aoi_grid(xs_a, ys_a, RING, 4326)
    gb = aoi_grid(xs_b, ys_b, RING, 4326)
    assert ga is not None and gb is not None

    A = np.full((len(ys_a), len(xs_a)), np.nan)
    B = np.full((len(ys_b), len(xs_b)), np.nan)
    A[int(np.argmin(abs(ys_a - 28.4804))), int(np.argmin(abs(xs_a - 85.5604)))] = 42.0
    B[int(np.argmin(abs(ys_b - 28.4804))), int(np.argmin(abs(xs_b - 85.5604)))] = 42.0

    oa, ob = place_on_grid(A, ga), place_on_grid(B, gb)
    assert oa.shape == ob.shape
    assert np.array_equal(np.argwhere(oa == 42.0), np.argwhere(ob == 42.0))


def test_the_two_rasters_can_simply_be_subtracted():
    ga = aoi_grid(*_centres(85.5588, 28.4812, 3, 3), RING, 4326)
    gb = aoi_grid(*_centres(85.5604, 28.4804, 4, 4), RING, 4326)
    assert ga is not None and gb is not None
    oa = place_on_grid(np.ones((3, 3)) * 5.0, ga)
    ob = place_on_grid(np.ones((4, 4)) * 3.0, gb)
    d = ob - oa
    assert d.shape == oa.shape
    assert np.nanmin(d) == -2.0


def test_a_product_off_the_lattice_is_refused_not_resampled():
    xs, ys = _centres(85.5588 + RES / 2, 28.4812, 3, 3)
    assert aoi_grid(xs, ys, RING, 4326) is None


def test_a_degenerate_product_grid_is_refused():
    assert aoi_grid(None, None, RING, 4326) is None
    assert aoi_grid(np.array([85.4]), np.array([28.6]), RING, None) is None


def test_a_valid_grid_requires_an_explicit_crs():
    xs, ys = _centres(85.5588, 28.4812, 3, 3)
    with pytest.raises(ValueError, match="without source EPSG/CRS"):
        aoi_grid(xs, ys, RING, None)
