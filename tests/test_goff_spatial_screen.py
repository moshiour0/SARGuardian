from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from src.goff_spatial_screen import (
    ScreenError,
    box_sum,
    interpolate_regular_volume,
    parse_identity,
    score_pair_arrays,
    validate_series,
)


def granule_name(
    subswath: str = "024",
    crid: str = "P05023",
    track: str = "098",
) -> str:
    return (
        "NISAR_L2_PR_GOFF_023_"
        f"{track}_A_016_{subswath}_4000_SH_"
        "20260620T233921_20260620T233956_"
        "20260702T233920_20260702T233955_"
        f"{crid}_N_F_J_001.h5"
    )


def test_granule_parser_extracts_first_and_third_acquisition_timestamps() -> None:
    identity = parse_identity(Path(granule_name()))
    assert identity.track == "098"
    assert identity.direction == "A"
    assert identity.frame == "016"
    assert identity.subswath == "024"
    assert identity.bandwidth == "4000"
    assert identity.polarization == "SH"
    assert identity.crid == "P05023"
    assert identity.reference == "20260620T233921"
    assert identity.secondary == "20260702T233920"


def test_series_identity_allows_swath_change_but_not_release_or_track_change() -> None:
    a = parse_identity(Path(granule_name("024", "P05023", "098")))
    b = parse_identity(Path(granule_name("028", "P05023", "098").replace(
        "20260620T233921_20260620T233956_20260702T233920_20260702T233955",
        "20260726T233919_20260726T233954_20260819T233918_20260819T233953",
    )))
    assert a.series_key == b.series_key
    x = parse_identity(Path(granule_name("028", "X05026", "098")))
    assert x.series_key != a.series_key
    y = parse_identity(Path(granule_name("028", "P05023", "170")))
    assert y.series_key != a.series_key


def test_box_sum_matches_explicit_constant_zero_padded_window() -> None:
    values = np.arange(35, dtype=np.float32).reshape(5, 7)
    got = box_sum(values, 3)
    padded = np.pad(values.astype(float), 1)
    expected = np.zeros_like(values, dtype=float)
    for row in range(values.shape[0]):
        for col in range(values.shape[1]):
            expected[row, col] = padded[row:row + 3, col:col + 3].sum()
    assert got.shape == values.shape
    np.testing.assert_allclose(got, expected, rtol=0, atol=1e-5)


def test_trilinear_geometry_interpolation_is_bounded_and_exact_for_constant_volume() -> None:
    volume = np.full((2, 2, 2), 7.25, dtype=np.float32)
    elevation = np.array([[50.0, 50.0], [50.0, 50.0]], dtype=np.float32)
    values = interpolate_regular_volume(
        volume,
        np.array([0.0, 100.0]),
        np.array([0.0, 1.0]),
        np.array([0.0, 1.0]),
        elevation,
        np.array([0.25, 0.75]),
        np.array([0.75, 0.25]),
    )
    np.testing.assert_allclose(values, 7.25, rtol=0, atol=1e-6)

    outside = interpolate_regular_volume(
        volume,
        np.array([0.0, 100.0]),
        np.array([0.0, 1.0]),
        np.array([0.0, 1.0]),
        elevation,
        np.array([0.25, 2.0]),
        np.array([0.75, 0.25]),
    )
    assert np.isnan(outside[:, 1]).all()


def _fixture_arrays(seed: int = 4) -> tuple[np.ndarray, ...]:
    rng = np.random.default_rng(seed)
    shape = (151, 151)
    east = rng.normal(0.0, 0.015, shape).astype(np.float32)
    north = rng.normal(0.0, 0.015, shape).astype(np.float32)
    snr = np.full(shape, 12.0, dtype=np.float32)
    correlation = np.full(shape, 0.85, dtype=np.float32)
    aoi = np.ones(shape, dtype=bool)
    east[65:76, 65:76] += 1.7
    north[65:76, 65:76] += 0.2
    return east, north, snr, correlation, aoi


def test_target_free_screen_ranks_a_synthetic_persistent_coherent_patch() -> None:
    east, north, snr, correlation, aoi = _fixture_arrays()
    result = score_pair_arrays(east, north, snr, correlation, aoi)
    assert result["candidate_pixels"] >= 10
    assert result["candidate"][70, 70]
    assert result["valid_aoi_pixels"] == aoi.size


def test_noise_only_and_low_quality_patch_do_not_become_candidates() -> None:
    rng = np.random.default_rng(9)
    shape = (151, 151)
    east = rng.normal(0.0, 0.02, shape).astype(np.float32)
    north = rng.normal(0.0, 0.02, shape).astype(np.float32)
    snr = np.full(shape, 10.0, dtype=np.float32)
    correlation = np.full(shape, 0.9, dtype=np.float32)
    aoi = np.ones(shape, dtype=bool)
    result = score_pair_arrays(east, north, snr, correlation, aoi)
    assert result["candidate_pixels"] == 0

    east[65:76, 65:76] += 1.8
    snr[65:76, 65:76] = 1.0
    correlation[65:76, 65:76] = 0.1
    rejected = score_pair_arrays(east, north, snr, correlation, aoi)
    assert rejected["candidate_pixels"] == 0


def test_aoi_gate_prevents_outside_patch_from_becoming_a_target_free_candidate() -> None:
    east, north, snr, correlation, _ = _fixture_arrays()
    aoi = np.zeros(east.shape, dtype=bool)
    aoi[:40, :40] = True
    result = score_pair_arrays(east, north, snr, correlation, aoi)
    assert result["candidate_pixels"] == 0
    assert not result["candidate"][70, 70]


def test_series_gate_refuses_mixed_crids_before_temporal_combination() -> None:
    # validate_series is tested through identity metadata before any inversion.
    # Constructing the raw rasters is deliberately unnecessary for this gate.
    from src.goff_spatial_screen import RawPair
    from rasterio.crs import CRS
    from rasterio.transform import Affine

    identity = parse_identity(Path(granule_name()))
    mismatch = parse_identity(Path(granule_name("028", "X05026").replace(
        "20260620T233921_20260620T233956_20260702T233920_20260702T233955",
        "20260726T233919_20260726T233954_20260819T233918_20260819T233953",
    )))
    arr = np.zeros((3, 3), dtype=np.float32)
    mask = np.ones((3, 3), dtype=bool)
    transform = Affine(80, 0, 0, 0, -80, 240)
    def pair(i):
        return RawPair(
            identity=i, east_offsets_m=arr, north_offsets_m=arr, snr=arr + 10,
            correlation=arr + 0.8, xs=np.array([40, 120, 200.]), ys=np.array([200, 120, 40.]),
            crs=CRS.from_epsg(32645), transform=transform, grid_mask=mask,
            source_sha256="0" * 64, source_bytes=1, geometry_note="synthetic",
        )
    with pytest.raises(ScreenError, match="MIXED_SERIES_IDENTITY"):
        validate_series([pair(identity), pair(mismatch)])
