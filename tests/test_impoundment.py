"""
The ranking has to favour what it says it favours.

Why this exists
===============
`--rank-by efficiency` is documented as favouring "sites a small landslide
could dam", and the whole Blatten argument rests on it: the observed lake there
was about 10 m deep, so a screening tool that cannot rank small-blockage
sensitivity has not screened for the thing that happened.

It ranked on volume-per-metre alone, evaluated at the smallest height each site
happened to respond to. Those heights differ between sites and impounded volume
grows faster than linearly with depth, so dividing by a LARGER threshold can
still give the bigger ratio. In the committed source-zone run the reach 790 m
from the failure point - which holds nothing at all until a 100 m blockage -
came SECOND of twelve, above four reaches that start filling at 10 m.

That is the opposite of the stated intent, and the module had no tests, so
nothing contradicted the docstring.

The rest of these cover the hydrology the volumes rest on, on a DEM small
enough to state by hand - the module otherwise needs several thousand external
elevation queries to run at all, which is why none of it was ever exercised.
"""

from __future__ import annotations

import json
import math
import urllib.parse

import impoundment
import numpy as np
import pytest

from impoundment import (AOIS, API_ENDPOINT, dem_from_api, dem_from_file,
                         dem_provenance, d8_receivers, elevation_grid_sha256,
                         fill_depressions, flow_accumulation, rank_sites,
                         file_sha256, write_geojson, write_provenance)


def site(lat, lon, by_height, name=""):
    """A site dict shaped exactly as analyse() builds one."""
    bh = {float(h): {"volume_Mm3": v, "area_km2": v / 10, "backwater_km": 1.0,
                     "truncated": False} for h, v in by_height.items()}
    lo = min(bh)
    return {"lat": lat, "lon": lon, "name": name, "by_height": bh,
            "sill_m": 3000.0, "upstream_cells": 100, "upstream_area_km2": 1.0,
            "max_volume_Mm3": max(v["volume_Mm3"] for v in bh.values()),
            "threshold_height_m": lo,
            "efficiency_Mm3_per_m": bh[lo]["volume_Mm3"] / lo,
            "volume_at_min_height_Mm3": bh[lo]["volume_Mm3"],
            "truncated": False}


# ---------------------------------------------------------------------------
# The ranking
# ---------------------------------------------------------------------------
def test_a_site_needing_a_hundred_metres_does_not_outrank_one_that_fills_at_ten():
    """
    The exact case from the committed source-zone run, in the numbers it was
    published with: site 2 responds only at 100 m (3.23 Mm3) and site 9
    responds at 10 m (0.21 Mm3).

        volume-per-metre    3.23/100 = 0.0323   vs   0.21/10 = 0.0209

    so the old metric put the 100 m site first - while the column the README
    itself calls "the one that actually separates them" is how many reaches
    respond at 10 m.
    """
    needs_100 = site(28.2944, 85.5308, {100: 3.23, 150: 10.14}, "needs 100 m")
    fills_at_10 = site(28.3495, 85.5474, {10: 0.21, 25: 0.65, 50: 2.37,
                                          100: 7.31, 150: 14.28}, "fills at 10 m")

    assert needs_100["efficiency_Mm3_per_m"] > fills_at_10["efficiency_Mm3_per_m"], (
        "fixture must reproduce the inversion, or it is not testing anything")

    ranked = rank_sites([needs_100, fills_at_10], "efficiency",
                        top_n=12, min_separation_km=0.0)
    assert [r["name"] for r in ranked] == ["fills at 10 m", "needs 100 m"]


def test_every_small_threshold_site_outranks_every_large_threshold_site():
    """Not just the pair above - the property, across a mixed table."""
    small = [site(28.0 + i * 0.1, 85.0, {10: 0.1 + i * 0.05, 150: 20.0}, f"s{i}")
             for i in range(3)]
    large = [site(29.0 + i * 0.1, 85.0, {150: 60.0 + i}, f"L{i}") for i in range(3)]

    ranked = rank_sites(large + small, "efficiency", top_n=12, min_separation_km=0.0)
    thresholds = [r["threshold_height_m"] for r in ranked]
    assert thresholds == sorted(thresholds), (
        f"threshold heights must be non-decreasing down the table, got {thresholds}")
    assert all(r["name"].startswith("s") for r in ranked[:3])


def test_within_one_threshold_the_more_efficient_site_wins():
    """Efficiency still orders the table - it just no longer crosses thresholds."""
    a = site(28.0, 85.0, {25: 3.0}, "efficient")
    b = site(28.5, 85.0, {25: 0.5}, "less efficient")
    ranked = rank_sites([b, a], "efficiency", top_n=12, min_separation_km=0.0)
    assert [r["name"] for r in ranked] == ["efficient", "less efficient"]


def test_volume_ranking_is_unchanged_and_still_available():
    """--rank-by volume answers a different question and must keep answering it."""
    big = site(28.0, 85.0, {150: 80.0}, "big")
    small = site(28.5, 85.0, {10: 0.3, 150: 20.0}, "small but responsive")
    ranked = rank_sites([small, big], "volume", top_n=12, min_separation_km=0.0)
    assert ranked[0]["name"] == "big"


def test_near_duplicate_reaches_are_suppressed():
    """
    Consecutive channel cells describe the same pool with the dam nudged a few
    tens of metres. Without suppression one valley fills the whole table.
    """
    close = [site(28.0 + i * 0.001, 85.0, {10: 1.0 - i * 0.01}, f"c{i}")
             for i in range(5)]
    far = site(28.5, 85.0, {10: 0.5}, "far")
    ranked = rank_sites(close + [far], "efficiency", top_n=12, min_separation_km=1.5)
    assert len(ranked) == 2
    assert {r["name"] for r in ranked} == {"c0", "far"}


def test_suppression_keeps_the_better_site_of_a_cluster():
    """It must drop the neighbours of the winner, not the winner."""
    weak = site(28.000, 85.0, {25: 5.0}, "weak")
    strong = site(28.002, 85.0, {10: 1.0}, "strong")
    ranked = rank_sites([weak, strong], "efficiency", top_n=12, min_separation_km=1.5)
    assert [r["name"] for r in ranked] == ["strong"]


def test_top_n_is_respected():
    many = [site(28.0 + i * 0.1, 85.0, {10: 1.0}, f"s{i}") for i in range(20)]
    assert len(rank_sites(many, "efficiency", top_n=12, min_separation_km=0.0)) == 12


# ---------------------------------------------------------------------------
# The hydrology the volumes rest on
# ---------------------------------------------------------------------------
def _valley(ny=24, nx=24, wall=3.0, drop=2.0):
    """A V-shaped valley draining toward increasing row."""
    r = np.arange(ny)[:, None]
    c = np.arange(nx)[None, :]
    return (ny - 1 - r) * drop + np.abs(c - nx // 2) * wall


def test_filling_never_lowers_the_ground():
    z = _valley()
    z[10, 12] -= 25.0                       # a pit
    filled = fill_depressions(z)
    assert np.all(filled >= z - 1e-9)
    assert filled[10, 12] > z[10, 12], "the pit must be filled for routing"


def test_filling_leaves_terrain_without_pits_alone():
    z = _valley()
    assert np.allclose(fill_depressions(z), z, atol=1e-9)


def test_flow_accumulates_downstream_and_conserves_cells():
    z = _valley()
    filled = fill_depressions(z)
    recv = d8_receivers(filled, dx=30.0, dy=30.0)
    acc = flow_accumulation(filled, recv)

    assert acc.min() >= 1, "every cell drains at least itself"
    # Accumulation must grow toward the outlet, never shrink going downstream.
    ny, nx = z.shape
    for idx, r in enumerate(recv):
        if r >= 0:
            assert acc.ravel()[r] >= acc.ravel()[idx]
    assert acc.max() > 0.5 * z.size, "a single valley should concentrate flow"


def test_a_cell_on_the_edge_drains_off_the_grid():
    """Boundary cells have no lower neighbour inside the grid and terminate."""
    z = _valley()
    filled = fill_depressions(z)
    recv = d8_receivers(filled, dx=30.0, dy=30.0)
    ny, nx = z.shape
    outlet_row = recv.reshape(ny, nx)[-1, :]
    assert (outlet_row < 0).any(), "the downstream edge must contain a terminus"


def test_diagonal_steps_are_not_treated_as_cardinal():
    """
    A diagonal move is sqrt(2) cells long, so a diagonal drop must be steeper
    than a cardinal one to win. Getting this wrong biases every channel toward
    the diagonals.
    """
    z = np.array([[10.0, 9.4, 20.0],
                  [20.0, 20.0, 20.0],
                  [20.0, 20.0, 20.0]])
    # From (0,0): cardinal east drops 0.6 over 1 cell; there is no better move.
    recv = d8_receivers(z, dx=1.0, dy=1.0)
    assert recv[0] == 1, "steepest descent per unit LENGTH, not per step"


def test_diagonal_distance_uses_rectangular_cell_geometry_exactly():
    # At (0,0), the eastward slope is 0.305 m/m. The diagonal drop is 1.09 m:
    # its slope is 1.09/sqrt(2^2+3^2) = 0.3023, so east must win. The old
    # sqrt(2)*mean(dx,dy) approximation made the diagonal slope ~0.3085 and
    # selected the wrong receiver.
    z = np.array([
        [10.0, 9.39, 30.0],
        [30.0, 8.91, 30.0],
        [30.0, 30.0, 30.0],
    ])
    receivers = d8_receivers(z, dx=2.0, dy=3.0)
    assert receivers[0] == 1


def test_api_dem_cache_reuses_exact_grid_without_a_second_network_call(
    tmp_path, monkeypatch
):
    import impoundment

    calls = []
    requested_points = []

    class Response:
        def __init__(self, payload):
            self.payload = payload
        def read(self):
            return self.payload

    def fake_urlopen(url, timeout=0):
        calls.append((url, timeout))
        query = urllib.parse.urlsplit(url).query
        locations = urllib.parse.parse_qs(query)["locations"][0]
        points = urllib.parse.unquote(locations).split("|")
        requested_points.extend(
            (float(point.split(",")[0]), float(point.split(",")[1]))
            for point in points
        )
        payload = json.dumps({
            "results": [{"elevation": 100.0 + i} for i in range(len(points))]
        }).encode("utf-8")
        return Response(payload)

    monkeypatch.setattr(impoundment.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(impoundment.time, "sleep", lambda *_: None)
    bbox = AOIS["source"]
    cache = tmp_path / "dem" / "source_srtm.npz"

    first_z, first_lats, first_lons = dem_from_api(bbox, target_cells=12, cache_path=cache)
    first_sha = elevation_grid_sha256(
        first_z, first_lats, first_lons, bbox, 12, endpoint=API_ENDPOINT
    )
    assert cache.is_file()
    assert len(calls) == 1
    expected_points = [
        (float(lat), float(lon))
        for lat in first_lats for lon in first_lons
    ]
    assert requested_points == expected_points, (
        "Cached axes must identify the exact coordinates sent to the API, not the unrounded linspace."
    )
    assert all(bbox[1] <= lat <= bbox[3] and bbox[0] <= lon <= bbox[2]
               for lat, lon in requested_points)

    def no_network(*args, **kwargs):
        raise AssertionError("a valid existing cache must not trigger a network request")

    monkeypatch.setattr(impoundment.urllib.request, "urlopen", no_network)
    second_z, second_lats, second_lons = dem_from_api(bbox, target_cells=12, cache_path=cache)
    second_sha = elevation_grid_sha256(
        second_z, second_lats, second_lons, bbox, 12, endpoint=API_ENDPOINT
    )

    assert first_sha == second_sha
    assert np.array_equal(first_z, second_z)
    assert np.array_equal(first_lats, second_lats)
    assert np.array_equal(first_lons, second_lons)
    assert len(calls) == 1


def test_api_dem_cache_fails_closed_when_grid_bytes_change(tmp_path, monkeypatch):
    import impoundment

    class Response:
        def read(self):
            return json.dumps({"results": [{"elevation": 100.0 + i} for i in range(12)]}).encode("utf-8")

    monkeypatch.setattr(impoundment.urllib.request, "urlopen", lambda *a, **k: Response())
    monkeypatch.setattr(impoundment.time, "sleep", lambda *_: None)
    bbox = AOIS["source"]
    cache = tmp_path / "source.npz"
    z, lats, lons = dem_from_api(bbox, target_cells=12, cache_path=cache)
    with np.load(cache, allow_pickle=False) as data:
        metadata_json = data["metadata_json"].copy()
    z[0, 0] += 999.0
    np.savez_compressed(cache, z=z, lats=lats, lons=lons, metadata_json=metadata_json)

    def no_network(*args, **kwargs):
        raise AssertionError("corrupt cache must not fall back to live data")

    monkeypatch.setattr(impoundment.urllib.request, "urlopen", no_network)
    with pytest.raises(ValueError, match="fingerprint mismatch"):
        dem_from_api(bbox, target_cells=12, cache_path=cache)


def test_api_dem_cache_refuses_a_different_bbox(tmp_path, monkeypatch):
    import impoundment

    class Response:
        def read(self):
            return json.dumps({"results": [{"elevation": 100.0 + i} for i in range(12)]}).encode("utf-8")

    monkeypatch.setattr(impoundment.urllib.request, "urlopen", lambda *a, **k: Response())
    monkeypatch.setattr(impoundment.time, "sleep", lambda *_: None)
    cache = tmp_path / "source.npz"
    dem_from_api(AOIS["source"], target_cells=12, cache_path=cache)
    with pytest.raises(ValueError, match="bbox does not match"):
        dem_from_api(AOIS["langtang"], target_cells=12, cache_path=cache)


def test_api_dem_rejects_partial_results_instead_of_misaligning_elevations(
    monkeypatch
):
    import impoundment

    class Response:
        def read(self):
            # For 12 requested points, returning one value would silently
            # scramble/rescale the elevation grid without a count check.
            return json.dumps({"results": [{"elevation": 123.0}]}).encode()

    monkeypatch.setattr(impoundment.urllib.request, "urlopen", lambda *a, **k: Response())
    monkeypatch.setattr(impoundment.time, "sleep", lambda *_: None)
    with pytest.raises(RuntimeError, match="returned 1 results for 12 requested cells"):
        dem_from_api(AOIS["source"], target_cells=12)


def test_local_dem_is_cropped_to_aoi_bbox_and_keeps_lon_lat_axes(tmp_path):
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin

    path = tmp_path / "larger_dem.tif"
    data = np.arange(30 * 30, dtype=np.float32).reshape(30, 30)
    transform = from_origin(85.44, 28.36, 0.005, 0.005)
    with rasterio.open(
        path, "w", driver="GTiff", height=30, width=30, count=1,
        dtype="float32", crs="EPSG:4326", transform=transform, nodata=-9999,
    ) as ds:
        ds.write(data, 1)

    bbox = AOIS["source"]
    z, lats, lons = dem_from_file(path, bbox=bbox)
    assert z.shape == (lats.size, lons.size)
    assert z.shape[0] >= 2 and z.shape[1] >= 2
    assert bbox[1] <= lats.min() < lats.max() <= bbox[3]
    assert bbox[0] <= lons.min() < lons.max() <= bbox[2]


def test_local_dem_rejects_projected_crs_instead_of_treating_metres_as_degrees(tmp_path):
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin

    path = tmp_path / "projected_dem.tif"
    with rasterio.open(
        path, "w", driver="GTiff", height=4, width=4, count=1,
        dtype="float32", crs="EPSG:32645",
        transform=from_origin(300000, 3150000, 30, 30),
    ) as ds:
        ds.write(np.arange(16, dtype=np.float32).reshape(4, 4), 1)

    with pytest.raises(ValueError, match="requires a WGS84/EPSG:4326"):
        dem_from_file(path, bbox=AOIS["source"])


def test_dem_provenance_identifies_local_input_without_claiming_api_source(tmp_path):
    z = np.arange(12, dtype=float).reshape(3, 4)
    lats = np.array([28.3, 28.2, 28.1])
    lons = np.array([85.4, 85.5, 85.6, 85.7])
    bbox = [85.4, 28.1, 85.7, 28.3]
    provenance = dem_provenance(z, lats, lons, bbox=bbox, source="local GeoTIFF",
                                target_cells=z.size)
    assert provenance["endpoint"] is None
    assert provenance["source"] == "local GeoTIFF"
    assert len(provenance["grid_sha256"]) == 64


def test_geojson_output_and_provenance_sidecar_are_atomic_and_bound_to_the_grid(tmp_path):
    output = tmp_path / "nested" / "dams.geojson"
    row = site(28.30, 85.50, {10: 0.2, 25: 0.7}, "fixture site")
    write_geojson([row], output)
    parsed = json.loads(output.read_text(encoding="utf-8"))
    assert parsed["type"] == "FeatureCollection"
    assert len(parsed["features"]) == 1
    assert parsed["features"][0]["geometry"]["coordinates"] == [85.50, 28.30]
    assert output.read_text(encoding="utf-8").endswith("\n")

    z = np.arange(12, dtype=float).reshape(3, 4)
    lats = np.array([28.3, 28.2, 28.1])
    lons = np.array([85.4, 85.5, 85.6, 85.7])
    bbox = [85.4, 28.1, 85.7, 28.3]
    provenance = dem_provenance(
        z, lats, lons, bbox=bbox, source="local GeoTIFF", target_cells=z.size
    )
    sidecar = write_provenance(
        output, provenance, {"aoi": "fixture", "heights_m": [10, 25]}, [row]
    )
    saved = json.loads(sidecar.read_text(encoding="utf-8"))
    assert saved["grid_sha256"] == provenance["grid_sha256"]
    assert saved["geojson_sha256"] == file_sha256(output)
    assert len(saved["result_rows_sha256"]) == 64
    assert saved["result_features"] == 1


def test_geojson_refuses_nonfinite_site_coordinates_without_creating_output(tmp_path):
    output = tmp_path / "invalid.geojson"
    bad = site(float("nan"), 85.5, {10: 0.2}, "bad coordinate")
    with pytest.raises(ValueError, match="invalid WGS84 site coordinates"):
        write_geojson([bad], output)
    assert not output.exists()


def test_api_geojson_cli_requires_an_explicit_reusable_dem_cache(
    tmp_path, monkeypatch, capsys
):
    import impoundment
    import sys

    monkeypatch.setattr(sys, "argv", [
        "impoundment.py", "--api", "--aoi", "source",
        "--geojson", str(tmp_path / "dams.geojson"),
    ])
    with pytest.raises(SystemExit) as exc:
        impoundment.main()
    assert exc.value.code == 2
    assert "--api-cache" in capsys.readouterr().err


def test_d8_receiver_rejects_invalid_postings_and_degenerate_grids():
    with pytest.raises(ValueError, match="finite positive metre distances"):
        d8_receivers(np.ones((3, 3)), dx=0.0, dy=1.0)
    with pytest.raises(ValueError, match="at least two rows and columns"):
        d8_receivers(np.ones((1, 3)), dx=1.0, dy=1.0)


def test_impoundment_analysis_rejects_mismatched_axes_and_duplicate_heights():
    z = np.arange(12, dtype=float).reshape(3, 4)
    lats = np.array([28.3, 28.2, 28.1])
    lons = np.array([85.4, 85.5, 85.6, 85.7])
    with pytest.raises(ValueError, match="does not match"):
        impoundment.analyse(z, lats[:-1], lons, [10, 25], 1, 2)
    with pytest.raises(ValueError, match=r"unique, finite (?:and )?positive"):
        impoundment.analyse(z, lats, lons, [10, 10], 1, 2)


def test_rank_sites_rejects_invalid_parameters_and_nonfinite_coordinates():
    sample = site(28.30, 85.50, {10: 0.2, 25: 0.7}, "fixture")
    with pytest.raises(ValueError, match="top_n"):
        rank_sites([sample], "volume", 0, 1.5)
    with pytest.raises(ValueError, match="min_separation"):
        rank_sites([sample], "volume", 1, -1.0)
    sample["lat"] = float("nan")
    with pytest.raises(ValueError, match="Invalid WGS84"):
        rank_sites([sample], "volume", 1, 1.5)


def test_api_coordinate_grid_refuses_resolution_finer_than_request_precision():
    from impoundment import _api_coordinate_axes

    # Four-by-four cells inside a single 0.00001-degree interval would collapse
    # onto duplicated coordinate requests after five-decimal serialization.
    with pytest.raises(ValueError, match="finer than the API"):
        _api_coordinate_axes((85.0, 28.0, 85.00001, 28.00001), 4, 4)


def test_elevation_grid_rejects_infinite_samples_but_keeps_nan_nodata():
    from impoundment import _validate_elevation_grid

    z = np.array([[100.0, np.nan], [101.0, 102.0]])
    lats = np.array([28.1, 28.0])
    lons = np.array([85.1, 85.2])
    _validate_elevation_grid(z, lats, lons, (85.1, 28.0, 85.2, 28.1))

    z[0, 0] = np.inf
    with pytest.raises(ValueError, match="infinite elevations"):
        _validate_elevation_grid(z, lats, lons, (85.1, 28.0, 85.2, 28.1))



def test_analyse_runs_end_to_end_with_a_valid_rectangular_dem(capsys):
    """A live DEM must pass grid validation, logging, routing, and ranking."""
    import impoundment

    rows, cols = np.mgrid[0:7, 0:8]
    z = (rows * 3.0 + cols * 2.0).astype(float)
    lats = np.linspace(28.30, 28.24, z.shape[0])
    lons = np.linspace(85.48, 85.55, z.shape[1])
    result = impoundment.analyse(
        z, lats, lons, heights=[10, 25], min_accum_cells=1, top_n=3,
        min_separation_km=0.0, rank_by="efficiency",
    )
    captured = capsys.readouterr()
    assert "Grid 7 x 8" in captured.out or "Grid 7 x 8" in captured.err
    assert isinstance(result, list)
    assert result, "the sloping fixture should produce at least one candidate pool"
    for site in result:
        assert np.isfinite(site["lat"]) and np.isfinite(site["lon"])
        assert site["threshold_height_m"] in {10.0, 25.0}
        assert np.isfinite(site["max_volume_Mm3"])
