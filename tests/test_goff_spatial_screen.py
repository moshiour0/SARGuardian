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
    seconds_variant = granule_name("028", "P05023", "098").replace(
        "20260620T233921_20260620T233956_20260702T233920_20260702T233955",
        "20260620T233922_20260620T233957_20260702T233921_20260702T233956",
    )
    same_interval = parse_identity(Path(seconds_variant))
    assert same_interval.series_key == a.series_key
    assert same_interval.pair_key == a.pair_key
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

    seconds_variant = parse_identity(Path(granule_name("028", "P05023").replace(
        "20260620T233921_20260620T233956_20260702T233920_20260702T233955",
        "20260620T233922_20260620T233957_20260702T233921_20260702T233956",
    )))
    with pytest.raises(ScreenError, match="DUPLICATE_ACQUISITION_INTERVALS"):
        validate_series([pair(identity), pair(seconds_variant)])


def test_raw_hdf5_to_ranked_candidate_outputs_end_to_end(tmp_path: Path) -> None:
    import h5py
    import rasterio
    from pyproj import Transformer
    from rasterio.transform import Affine

    from src.goff_spatial_screen import run_screen

    height, width = 64, 64
    xs = 500000.0 + np.arange(width, dtype=np.float64) * 80.0
    ys = 3150000.0 - np.arange(height, dtype=np.float64) * 80.0
    transform = Affine(80.0, 0.0, xs[0] - 40.0, 0.0, -80.0, ys[0] + 40.0)
    dem_path = tmp_path / "dem.tif"
    with rasterio.open(
        dem_path, "w", driver="GTiff", width=width, height=height, count=1,
        dtype="float32", crs="EPSG:32645", transform=transform, nodata=-9999.0,
    ) as dataset:
        dataset.write(np.full((height, width), 4500.0, dtype=np.float32), 1)

    inverse = Transformer.from_crs("EPSG:32645", "EPSG:4326", always_xy=True)
    west, south = inverse.transform(xs[8], ys[52])
    east, north = inverse.transform(xs[55], ys[8])
    ring_xy = [
        inverse.transform(xs[8], ys[8]),
        inverse.transform(xs[55], ys[8]),
        inverse.transform(xs[55], ys[52]),
        inverse.transform(xs[8], ys[52]),
        inverse.transform(xs[8], ys[8]),
    ]
    aoi_wkt = "POLYGON ((" + ", ".join(f"{lon:.8f} {lat:.8f}" for lon, lat in ring_xy) + "))"

    def write_goff(path: Path, granule: str, reference: str, secondary: str, seed: int) -> None:
        rng = np.random.default_rng(seed)
        at = rng.normal(0.0, 0.01, (height, width)).astype(np.float32)
        sr = rng.normal(0.0, 0.005, (height, width)).astype(np.float32)
        at[28:34, 28:34] += 2.0
        snr = np.full((height, width), 12.0, dtype=np.float32)
        corr = np.full((height, width), 0.9, dtype=np.float32)
        parent = "science/LSAR/GOFF/grids/frequencyA/pixelOffsets/HH/layer2"
        radar = "science/LSAR/GOFF/metadata/radarGrid"
        with h5py.File(path, "w") as h5:
            group = h5.require_group(parent)
            group.create_dataset("alongTrackOffset", data=at)
            group.create_dataset("slantRangeOffset", data=sr)
            group.create_dataset("snr", data=snr)
            group.create_dataset("correlationSurfacePeak", data=corr)
            group.create_dataset("xCoordinates", data=xs)
            group.create_dataset("yCoordinates", data=ys)
            projection = group.create_dataset("projection", data=np.int32(32645))
            projection.attrs["epsg_code"] = 32645

            metadata = h5.require_group(radar)
            hgt = np.array([0.0, 9000.0])
            gy = np.linspace(ys[-1], ys[0], 3)  # radarGrid Y is ascending
            gx = np.linspace(xs[0], xs[-1], 3)
            metadata.create_dataset("heightAboveEllipsoid", data=hgt)
            metadata.create_dataset("xCoordinates", data=gx)
            metadata.create_dataset("yCoordinates", data=gy)
            shape = (2, 3, 3)
            metadata.create_dataset("incidenceAngle", data=np.full(shape, 40.0))
            metadata.create_dataset("losUnitVectorX", data=np.full(shape, -0.5))
            metadata.create_dataset("losUnitVectorY", data=np.full(shape, 0.2))
            metadata.create_dataset("alongTrackUnitVectorX", data=np.full(shape, 0.8))
            metadata.create_dataset("alongTrackUnitVectorY", data=np.full(shape, 0.6))

    products = []
    specifications = [
        ("023", "024", "20260620T000000", "20260702T000000", 21),
        ("026", "028", "20260726T000000", "20260819T000000", 22),
    ]
    for granule, subswath, reference, secondary, seed in specifications:
        name = (
            f"NISAR_L2_PR_GOFF_{granule}_098_A_016_{subswath}_4000_SH_"
            f"{reference}_{reference[:9]}000030_{secondary}_{secondary[:9]}000030_"
            "P05023_N_F_J_001.h5"
        )
        path = tmp_path / name
        write_goff(path, granule, reference, secondary, seed)
        products.append(path)

    output = tmp_path / "screen-output"
    summary = run_screen(
        products, dem_path, aoi_wkt, output,
        local_window_px=15, direction_window_px=3, z_min=2.0,
        anomaly_min_m=0.3, coherence_min=0.6,
        minimum_pairs=2, minimum_component_pixels=5,
    )
    assert summary["status"] == "EXPLORATORY_RANKING_ONLY"
    assert summary["pair_count"] == 2
    assert summary["aoi_grid_pixels"] > 1000
    for filename in [
        "pair_persistence_count.tif", "mean_anomaly_m.tif",
        "mean_directional_coherence.tif", "persistent_candidate_mask.tif",
        "candidates.geojson", "candidates.csv", "screen-summary.json",
    ]:
        assert (output / filename).is_file()
    with rasterio.open(output / "pair_persistence_count.tif") as dataset:
        assert dataset.crs.to_epsg() == 32645
        assert dataset.shape == (height, width)
    import json
    geojson = json.loads((output / "candidates.geojson").read_text(encoding="utf-8"))
    for feature in geojson["features"]:
        coords = feature["geometry"]["coordinates"][0]
        assert all(-180 <= point[0] <= 180 and -90 <= point[1] <= 90 for point in coords)
