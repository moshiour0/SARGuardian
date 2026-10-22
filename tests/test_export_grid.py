"""
Exports must land on one grid, or they cannot be differenced.

The fixture uses an explicit WGS84 CRS and cell-centre coordinates aligned to
the common AOI lattice. Missing CRS must be tested separately, not implied by
passing epsg=None for projected coordinates.
"""
from __future__ import annotations

import numpy as np
import pytest

from gunw_reader import aoi_grid, export_clipped, place_on_grid

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


def _projected_source_grid():
    from pyproj import Transformer

    tf = Transformer.from_crs(4326, 32645, always_xy=True)
    cx, cy = tf.transform(85.52983, 28.27799)
    # UTM posting centres are aligned to the absolute 80 m lattice; this lets
    # the exporter use the deterministic common AOI grid without resampling.
    x_start = np.floor((cx - 8_000.0 - RES) / 80.0) * 80.0 + 40.0
    y_start = np.floor((cy + 8_000.0 - RES) / 80.0) * 80.0 + 40.0
    xs = x_start + 80.0 * np.arange(201)
    ys = y_start - 80.0 * np.arange(201)
    return xs, ys


def test_all_invalid_gunw_product_writes_an_explicit_nodata_aoi_raster(tmp_path):
    import rasterio

    xs, ys = _projected_source_grid()
    shape = (len(ys), len(xs))
    result = {
        "file": "NISAR_L2_PR_GUNW_024_098_A_016_025_4000_SH_"
                "20260702T233920_20260702T233955_20260714T233920_20260714T233954_P05023_F_F_J_001.h5",
        "reference_date": "20260702",
        "secondary_date": "20260714",
        "valid": np.zeros(shape, dtype=bool),
        "displacement_mm": np.full(shape, np.nan, dtype=np.float32),
        "coherence": np.full(shape, np.nan, dtype=np.float32),
        "xs": xs, "ys": ys, "epsg": 32645,
        "aoi_name": "source", "aoi_mask_applied": True,
        "wavelength_m": 0.2384, "ionosphere_correction_applied": True,
    }
    result_meta = export_clipped(result, tmp_path)
    assert result_meta["export_status"] == "NO_VALID_PIXELS"
    assert result_meta["export_valid_px"] == 0
    assert result_meta["export_file"]

    with rasterio.open(tmp_path / result_meta["export_file"]) as src:
        assert src.count == 2
        assert src.crs.to_epsg() == 32645
        assert src.tags()["aoi_name"] == "source"
        assert src.tags()["aoi_mask_applied"] == "true"
        assert src.tags()["export_status"] == "NO_VALID_PIXELS"
        assert int(src.tags()["valid_pixel_count"]) == 0
        assert np.isnan(src.read(1)).all()
        assert np.isnan(src.read(2)).all()


def test_gunw_export_refuses_to_claim_aoi_provenance_without_a_mask(tmp_path):
    xs, ys = _projected_source_grid()
    shape = (len(ys), len(xs))
    result = {
        "file": "NISAR_L2_PR_GUNW_024_098_A_016_025_4000_SH_"
                "20260702T233920_20260702T233955_20260714T233920_20260714T233954_P05023_F_F_J_001.h5",
        "reference_date": "20260702", "secondary_date": "20260714",
        "valid": np.ones(shape, dtype=bool),
        "displacement_mm": np.ones(shape, dtype=np.float32),
        "coherence": np.ones(shape, dtype=np.float32),
        "xs": xs, "ys": ys, "epsg": 32645,
        "aoi_name": "source", "aoi_mask_applied": False,
        "wavelength_m": 0.2384, "ionosphere_correction_applied": True,
    }
    with pytest.raises(RuntimeError, match="did not apply the AOI mask"):
        export_clipped(result, tmp_path)
    assert not list(tmp_path.glob("*.tif"))
