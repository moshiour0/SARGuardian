from __future__ import annotations

import pytest

from aoi_wkt import parse_polygon_wkt


def test_closed_and_open_wkt_are_normalized_to_closed_ring():
    open_ring = "POLYGON ((85 28, 86 28, 86 29, 85 29))"
    ring = parse_polygon_wkt(open_ring)
    assert ring[0] == ring[-1]
    assert len(ring) == 5
    assert ring[1] == (86.0, 28.0)


def test_srid_4326_is_accepted():
    ring = parse_polygon_wkt("SRID=4326;POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))")
    assert len(ring) == 5


@pytest.mark.parametrize(
    "wkt, code",
    [
        ("", "AOI_WKT_EMPTY"),
        ("POINT (85 28)", "AOI_WKT_MUST_BE_SINGLE_RING_POLYGON_WGS84_NO_HOLES"),
        ("MULTIPOLYGON (((0 0, 1 0, 1 1, 0 0)))", "AOI_WKT_MUST_BE_SINGLE_RING_POLYGON_WGS84_NO_HOLES"),
        ("SRID=3857;POLYGON ((0 0, 1 0, 1 1, 0 0))", "AOI_WKT_SRID_MUST_BE_4326"),
        ("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0), (0.2 0.2, 0.3 0.2, 0.3 0.3, 0.2 0.2))", "AOI_WKT_MUST_BE_SINGLE_RING_POLYGON_WGS84_NO_HOLES"),
        ("POLYGON ((0 0, 1 1, 0 1, 1 0, 0 0))", "AOI_WKT_SELF_INTERSECTING_RING"),
        ("POLYGON ((0 0, 200 0, 1 1, 0 0))", "AOI_WKT_COORDINATE_OUTSIDE_WGS84_BOUNDS"),
        ("POLYGON ((0 0 10, 1 0 10, 1 1 10, 0 0 10))", "AOI_WKT_COORDINATES_MUST_BE_LON_LAT_PAIRS"),
    ],
)
def test_invalid_wkt_is_rejected(wkt, code):
    with pytest.raises(ValueError, match=code):
        parse_polygon_wkt(wkt)


def test_degenerate_area_is_rejected():
    with pytest.raises(ValueError, match="AOI_WKT_POLYGON_AREA_IS_ZERO"):
        parse_polygon_wkt("POLYGON ((0 0, 1 1, 2 2, 0 0))")
