"""Strict parser for user-supplied WGS84 polygon AOIs.

Accepted format: POLYGON ((lon lat, lon lat, ...)) with one exterior ring.
Holes, MULTIPOLYGON, Z/M coordinates and non-WGS84 SRIDs are rejected rather
than silently changing geometry. Coordinates use longitude, latitude order.
"""
from __future__ import annotations

import math
import re
from typing import List, Tuple

Point = Tuple[float, float]

def parse_polygon_wkt(value: str, *, max_vertices: int = 500) -> List[Point]:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("AOI_WKT_EMPTY")
    wkt = value.strip()
    srid = re.match(r"(?is)^\s*SRID\s*=\s*(\d+)\s*;\s*(.*)$", wkt)
    if srid:
        if int(srid.group(1)) != 4326:
            raise ValueError("AOI_WKT_SRID_MUST_BE_4326")
        wkt = srid.group(2).strip()
    match = re.fullmatch(r"(?is)POLYGON\s*\(\s*\((.*?)\)\s*\)", wkt)
    if not match:
        raise ValueError("AOI_WKT_MUST_BE_SINGLE_RING_POLYGON_WGS84_NO_HOLES")
    body = match.group(1).strip()
    chunks = body.split(",")
    if len(chunks) > max_vertices:
        raise ValueError(f"AOI_WKT_TOO_MANY_VERTICES_MAX_{max_vertices}")
    ring: List[Point] = []
    for chunk in chunks:
        fields = chunk.split()
        if len(fields) != 2:
            raise ValueError("AOI_WKT_COORDINATES_MUST_BE_LON_LAT_PAIRS")
        try:
            lon, lat = float(fields[0]), float(fields[1])
        except (TypeError, ValueError, OverflowError):
            raise ValueError("AOI_WKT_COORDINATE_NOT_NUMERIC") from None
        if not math.isfinite(lon) or not math.isfinite(lat):
            raise ValueError("AOI_WKT_COORDINATE_NOT_FINITE")
        if not (-180 <= lon <= 180 and -90 <= lat <= 90):
            raise ValueError("AOI_WKT_COORDINATE_OUTSIDE_WGS84_BOUNDS")
        point = (lon, lat)
        if ring and point == ring[-1]:
            raise ValueError("AOI_WKT_HAS_DUPLICATE_CONSECUTIVE_VERTEX")
        ring.append(point)
    if len(ring) < 3:
        raise ValueError("AOI_WKT_NEEDS_AT_LEAST_THREE_DISTINCT_VERTICES")
    if ring[0] != ring[-1]:
        ring.append(ring[0])
    if len(set(ring[:-1])) < 3:
        raise ValueError("AOI_WKT_NEEDS_AT_LEAST_THREE_DISTINCT_VERTICES")
    twice_area = sum(
        ring[i][0] * ring[i + 1][1] - ring[i + 1][0] * ring[i][1]
        for i in range(len(ring) - 1)
    )
    if not math.isfinite(twice_area) or abs(twice_area) < 1e-12:
        raise ValueError("AOI_WKT_POLYGON_AREA_IS_ZERO")
    # Reject self-intersecting rings. Adjacent edges may share their endpoint.
    def orient(a: Point, b: Point, c: Point) -> float:
        return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
    def on_segment(a: Point, b: Point, p: Point) -> bool:
        return (min(a[0], b[0])-1e-12 <= p[0] <= max(a[0], b[0])+1e-12
                and min(a[1], b[1])-1e-12 <= p[1] <= max(a[1], b[1])+1e-12)
    edges = list(zip(ring[:-1], ring[1:]))
    n = len(edges)
    for i, (a, b) in enumerate(edges):
        for j in range(i + 1, n):
            if j == i + 1 or (i == 0 and j == n - 1):
                continue
            c, d = edges[j]
            o1, o2, o3, o4 = orient(a,b,c), orient(a,b,d), orient(c,d,a), orient(c,d,b)
            crosses = ((o1 > 1e-12 and o2 < -1e-12) or
                       (o1 < -1e-12 and o2 > 1e-12)) and (
                       (o3 > 1e-12 and o4 < -1e-12) or
                       (o3 < -1e-12 and o4 > 1e-12))
            touches = ((abs(o1) <= 1e-12 and on_segment(a,b,c)) or
                       (abs(o2) <= 1e-12 and on_segment(a,b,d)) or
                       (abs(o3) <= 1e-12 and on_segment(c,d,a)) or
                       (abs(o4) <= 1e-12 and on_segment(c,d,b)))
            if crosses or touches:
                raise ValueError("AOI_WKT_SELF_INTERSECTING_RING")
    return ring
