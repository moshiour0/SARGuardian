"""
impoundment.py
--------------
Find places where a landslide could dam a river, and how much water it would
hold.

Why this exists
===============
The 26 August 2026 Nepal cascade was not a single failure - it was ice/rock
detachment, then channel blockage, then breach, then a surge that raised the
Trishuli by metres in half an hour. The blockage step is the multiplier: it
converts a local slope failure into a downstream flood.

Blockage potential is a property of the *terrain*, so it can be mapped in
advance from a DEM alone, with no SAR at all. Cross that map with InSAR
deformation on the flanking slopes and you have a genuine two-factor alert:
a site that is both dammable and moving.

Method
======
1. Fill depressions (priority-flood) for routing purposes only.
2. D8 flow directions and flow accumulation.
3. Channel = cells whose upstream contributing area exceeds a threshold.
4. For each channel cell, treat it as a dam site. For a range of dam heights,
   flood the *upstream contributing area* to that water level and integrate
   the depth.
5. Rank by impounded volume; flag sites whose pool touches the AOI edge
   (truncated, so the volume is a lower bound).

DEM sources
===========
--dem FILE.tif    a local GeoTIFF (NASADEM, SRTM, Copernicus GLO-30). Proper.
--api             OpenTopoData SRTM 30 m over the AOI bbox. No download, no
                  credentials, rate-limited to ~100 points/second, so keep the
                  grid modest. Good enough to prototype; not production.

Limits worth stating in a paper
===============================
- A pool is only closed if it does not spill over a col into a neighbouring
  catchment. This tool does not test for cols. Volumes are therefore upper
  bounds for any site near a drainage divide.
- Dam height is imposed, not predicted. Whether a slope can actually deliver
  that volume of debris is a separate question - couple this to a runout model
  before calling a number a forecast.
- SRTM dates from 2000. In terrain reshaped since (the 2015 Gorkha avalanche
  through Langtang, for instance) the geometry is out of date.

Usage
-----
    python src/impoundment.py --api --aoi langtang --heights 50 100 150
    python src/impoundment.py --dem nasadem.tif --top 15 --geojson outputs/dams.geojson
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import logging
import math
import sys
import time
import urllib.parse
from datetime import datetime, timezone
import urllib.request
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger("impoundment")

AOIS = {
    # The confirmed 26 Aug 2026 source zone. Bounding box of SOURCE_RING in
    # gunw_reader; the failure point is at 28.28771 N, 85.52809 E.
    "source": (85.46452, 28.24529, 85.55619, 28.35290),
    "langtang": (85.45911, 28.24474, 85.56485, 28.32969),
    "lhende": (85.44, 28.34, 85.62, 28.47),
    # A second range, a second continent, and an outcome already on the record.
    # The Birch Glacier collapsed into the Loetschental on 28 May 2025, and the
    # deposit obstructed the Lonza and impounded a lake - a documented landslide
    # dam. Everything else here is calibrated against that event, so running the
    # susceptibility map over it is the one available test of whether this
    # travels beyond the Himalaya. Blatten village sits at 46.42028 N 7.81972 E.
    "blatten": (7.75, 46.38, 7.90, 46.46),
}

# D8 neighbour offsets (row, col) and their relative distances
NB = [(-1, 0, 1.0), (1, 0, 1.0), (0, -1, 1.0), (0, 1, 1.0),
      (-1, -1, math.sqrt(2)), (-1, 1, math.sqrt(2)),
      (1, -1, math.sqrt(2)), (1, 1, math.sqrt(2))]


# ---------------------------------------------------------------------------
# DEM loading
# ---------------------------------------------------------------------------
API_ENDPOINT = "https://api.opentopodata.org/v1/srtm30m"
API_CACHE_SCHEMA = 1


def api_grid_shape(bbox, target_cells: int = 12000) -> tuple[int, int, float]:
    """Return the deterministic north-to-south grid dimensions for an API AOI."""
    if len(bbox) != 4:
        raise ValueError("AOI bbox must be [west, south, east, north].")
    lon0, lat0, lon1, lat1 = map(float, bbox)
    if not (all(math.isfinite(v) for v in (lon0, lat0, lon1, lat1))
            and -180 <= lon0 < lon1 <= 180 and -90 <= lat0 < lat1 <= 90):
        raise ValueError(f"Invalid WGS84 bbox: {bbox!r}")
    if int(target_cells) != target_cells or target_cells < 4:
        raise ValueError("--cells must be an integer >= 4 so both grid axes contain at least two samples.")
    midlat = (lat0 + lat1) / 2
    w_km = (lon1 - lon0) * 111.32 * math.cos(math.radians(midlat))
    h_km = (lat1 - lat0) * 110.57
    aspect = w_km / h_km
    if not (math.isfinite(aspect) and aspect > 0):
        raise ValueError("AOI grid aspect ratio is invalid.")
    ny = max(2, int(round(math.sqrt(target_cells / aspect))))
    nx = max(2, int(round(target_cells / ny)))
    return ny, nx, w_km


def elevation_grid_sha256(z, lats, lons, bbox, target_cells, endpoint: str | None = None) -> str:
    """Stable digest of the exact elevation array, axes and request geometry."""
    digest = hashlib.sha256()
    descriptor = {
        "schema_version": API_CACHE_SCHEMA,
        "endpoint": endpoint,
        "bbox": [float(v) for v in bbox],
        "target_cells": int(target_cells),
        "shape": list(np.asarray(z).shape),
    }
    digest.update(json.dumps(descriptor, sort_keys=True, separators=(",", ":")).encode("utf-8"))
    for value in (z, lats, lons):
        array = np.ascontiguousarray(np.asarray(value, dtype="<f8"))
        digest.update(json.dumps(list(array.shape), separators=(",", ":")).encode("ascii"))
        digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def _validate_elevation_grid(z, lats, lons, bbox) -> None:
    z = np.asarray(z)
    lats = np.asarray(lats)
    lons = np.asarray(lons)
    if z.ndim != 2 or lats.ndim != 1 or lons.ndim != 1:
        raise ValueError("DEM cache must contain a 2-D elevation grid and 1-D axes.")
    if z.shape != (lats.size, lons.size) or min(z.shape) < 2:
        raise ValueError(f"DEM shape {z.shape} does not match coordinate axes ({lats.size}, {lons.size}).")
    if not (np.all(np.isfinite(lats)) and np.all(np.isfinite(lons)) and np.isfinite(z).any()):
        raise ValueError("DEM grid has invalid coordinate axes or no finite elevation samples.")
    if not (np.all(np.diff(lats) < 0) and np.all(np.diff(lons) > 0)):
        raise ValueError("DEM axes must be strictly north-to-south in latitude and west-to-east in longitude.")
    west, south, east, north = map(float, bbox)
    eps = 1e-9
    if (lons.min() < west - eps or lons.max() > east + eps
            or lats.min() < south - eps or lats.max() > north + eps):
        raise ValueError("DEM coordinate axes extend outside the requested WGS84 bbox.")


def _load_api_cache(path: Path, bbox, target_cells: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    try:
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata_json"].item()))
            z = np.asarray(data["z"], dtype=np.float64)
            lats = np.asarray(data["lats"], dtype=np.float64)
            lons = np.asarray(data["lons"], dtype=np.float64)
    except Exception as exc:
        raise ValueError(f"DEM cache {path} is unreadable; refusing to fall back silently to live API data.") from exc

    if metadata.get("schema_version") != API_CACHE_SCHEMA:
        raise ValueError(f"DEM cache schema mismatch in {path}.")
    if metadata.get("endpoint") != API_ENDPOINT:
        raise ValueError(f"DEM cache endpoint mismatch in {path}.")
    if metadata.get("bbox") != [float(v) for v in bbox]:
        raise ValueError(f"DEM cache AOI bbox does not match the requested AOI in {path}.")
    if metadata.get("target_cells") != int(target_cells):
        raise ValueError(f"DEM cache target-cell request does not match --cells in {path}.")
    ny, nx, _ = api_grid_shape(bbox, target_cells)
    if z.shape != (ny, nx):
        raise ValueError(f"DEM cache shape mismatch: expected {(ny, nx)}, got {z.shape}.")
    _validate_elevation_grid(z, lats, lons, bbox)
    fingerprint = elevation_grid_sha256(z, lats, lons, bbox, target_cells, endpoint=API_ENDPOINT)
    if metadata.get("grid_sha256") != fingerprint:
        raise ValueError(f"DEM cache fingerprint mismatch in {path}; the cached grid may be incomplete or modified.")
    logger.info("Loaded fixed DEM grid %s sha256=%s", path, fingerprint)
    return z, lats, lons


def _save_api_cache(path: Path, z, lats, lons, bbox, target_cells: int) -> None:
    _validate_elevation_grid(z, lats, lons, bbox)
    metadata = {
        "schema_version": API_CACHE_SCHEMA,
        "source": "OpenTopoData SRTM 30m",
        "endpoint": API_ENDPOINT,
        "bbox": [float(v) for v in bbox],
        "target_cells": int(target_cells),
        "shape": list(np.asarray(z).shape),
        "nodata_cells": int(np.count_nonzero(~np.isfinite(z))),
        "grid_sha256": elevation_grid_sha256(z, lats, lons, bbox, target_cells, endpoint=API_ENDPOINT),
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("wb") as stream:
        np.savez_compressed(
            stream, z=np.asarray(z, dtype=np.float64),
            lats=np.asarray(lats, dtype=np.float64),
            lons=np.asarray(lons, dtype=np.float64),
            metadata_json=np.array(json.dumps(metadata, sort_keys=True, separators=(",", ":"))),
        )
    temp.replace(path)
    logger.info("Saved fixed DEM grid %s sha256=%s", path, metadata["grid_sha256"])


def dem_from_api(bbox, target_cells: int = 12000, cache_path: Path | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sample SRTM 30 m through OpenTopoData, optionally pinning the exact grid.

    When cache_path exists it is verified and reused without a network request.
    A malformed/mismatched cache is an error, never an implicit cache miss.
    """
    ny, nx, w_km = api_grid_shape(bbox, target_cells)
    if cache_path is not None and cache_path.exists():
        return _load_api_cache(cache_path, bbox, target_cells)

    lon0, lat0, lon1, lat1 = map(float, bbox)
    h_km = (lat1 - lat0) * 110.57
    lats = np.linspace(lat1, lat0, ny)
    lons = np.linspace(lon0, lon1, nx)
    pts = [(la, lo) for la in lats for lo in lons]
    n_calls = (len(pts) + 99) // 100
    logger.info("Sampling %d x %d = %d cells (~%.0f m) via OpenTopoData: %d calls, ~%.0f s",
                ny, nx, len(pts), w_km * 1000 / nx, n_calls, n_calls * 1.15)

    elevations = []
    for k in range(0, len(pts), 100):
        chunk = pts[k:k + 100]
        locations = "|".join(f"{lat:.5f},{lon:.5f}" for lat, lon in chunk)
        url = API_ENDPOINT + "?locations=" + urllib.parse.quote(locations)
        for attempt in range(4):
            try:
                response = urllib.request.urlopen(url, timeout=120)
                payload = json.loads(response.read().decode("utf-8"))
                results = payload.get("results")
                if not isinstance(results, list) or len(results) != len(chunk):
                    raise ValueError(f"OpenTopoData returned {len(results) if isinstance(results, list) else 'no'} results for {len(chunk)} requested cells.")
                elevations.extend(item.get("elevation") for item in results)
                break
            except Exception as exc:
                if attempt == 3:
                    raise RuntimeError(f"OpenTopoData failed for request batch {k//100 + 1}/{n_calls}: {type(exc).__name__}: {exc}") from exc
                logger.warning("retrying OpenTopoData batch %d/%d after %s", k//100 + 1, n_calls, type(exc).__name__)
                time.sleep(3)
        if k + len(chunk) < len(pts):
            time.sleep(1.1)

    z = np.array([np.nan if value is None else value for value in elevations], dtype=np.float64).reshape(ny, nx)
    _validate_elevation_grid(z, lats, lons, bbox)
    if cache_path is not None:
        _save_api_cache(cache_path, z, lats, lons, bbox, target_cells)
    return z, lats, lons


def dem_from_file(path: Path, bbox=None):
    """Load only an axis-aligned north-up EPSG:4326 GeoTIFF for this grid solver.

    The solver computes metric postings from degree axes. Treating a projected
    UTM GeoTIFF as longitude/latitude would silently corrupt flow slopes,
    cell areas, and output site coordinates, so unsupported CRS/affines fail
    closed instead of being interpreted by assumption.
    """
    try:
        import rasterio
    except ImportError:
        logger.error("rasterio required for --dem: pip install rasterio")
        raise SystemExit(1)
    from rasterio.windows import Window

    with rasterio.open(path) as src:
        if src.crs is None or src.crs.to_epsg() != 4326:
            raise ValueError("The impoundment solver currently requires a WGS84/EPSG:4326 DEM; reproject the GeoTIFF before analysis.")
        transform = src.transform
        if (abs(transform.b) > 1e-10 or abs(transform.d) > 1e-10
                or transform.a <= 0 or transform.e >= 0):
            raise ValueError("DEM must be an unrotated north-up EPSG:4326 grid (positive x posting, negative y posting).")
        h, w = src.height, src.width
        if min(h, w) < 2:
            raise ValueError("DEM must contain at least two rows and columns.")
        full_xs = np.array([src.xy(0, col)[0] for col in range(w)], dtype=float)
        full_ys = np.array([src.xy(row, 0)[1] for row in range(h)], dtype=float)
        if bbox is None:
            west, south, east, north = (
                float(full_xs.min()), float(full_ys.min()),
                float(full_xs.max()), float(full_ys.max())
            )
        else:
            west, south, east, north = map(float, bbox)
        cols = np.flatnonzero((full_xs >= west) & (full_xs <= east))
        rows = np.flatnonzero((full_ys >= south) & (full_ys <= north))
        if cols.size < 2 or rows.size < 2:
            raise ValueError("DEM has fewer than two pixel centres in each axis inside the selected AOI bbox.")
        # Crop by pixel centres in their native geographic grid, not by a
        # rounded window that could reintroduce samples outside the AOI bbox.
        r0, c0 = int(rows[0]), int(cols[0])
        rh, cw = int(rows[-1] - rows[0] + 1), int(cols[-1] - cols[0] + 1)
        arr = src.read(1, window=Window(c0, r0, cw, rh)).astype(float)
        if src.nodata is not None:
            arr[arr == src.nodata] = np.nan
        xs, ys = full_xs[cols], full_ys[rows]
    _validate_elevation_grid(arr, ys, xs, (west, south, east, north))
    return arr, ys, xs


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dem_provenance(z, lats, lons, *, bbox, source: str,
                   target_cells: int, source_file: Path | None = None,
                   cache_path: Path | None = None) -> dict:
    """Machine-readable provenance for the exact elevation grid used."""
    _validate_elevation_grid(z, lats, lons, bbox)
    data = {
        "schema_version": 1,
        "source": source,
        "endpoint": API_ENDPOINT if source == "OpenTopoData SRTM 30m" else None,
        "bbox_wgs84": [float(v) for v in bbox],
        "shape_rows_cols": list(np.asarray(z).shape),
        "latitude_north_south": [float(np.max(lats)), float(np.min(lats))],
        "longitude_west_east": [float(np.min(lons)), float(np.max(lons))],
        "latitude_posting_deg": float(np.median(-np.diff(lats))),
        "longitude_posting_deg": float(np.median(np.diff(lons))),
        "elevation_nodata_cells": int(np.count_nonzero(~np.isfinite(z))),
        "elevation_min_m": float(np.nanmin(z)),
        "elevation_max_m": float(np.nanmax(z)),
        "grid_sha256": elevation_grid_sha256(
            z, lats, lons, bbox, target_cells,
            endpoint=API_ENDPOINT if source == "OpenTopoData SRTM 30m" else None,
        ),
    }
    if source_file is not None:
        data["source_file"] = str(source_file)
        data["source_file_sha256"] = file_sha256(source_file)
    if cache_path is not None:
        data["api_cache_file"] = str(cache_path)
    return data


def write_provenance(out: Path, provenance: dict, parameters: dict, rows: list[dict]) -> Path:
    """Write an atomic sidecar binding a GeoJSON result to input grid and parameters."""
    rows_payload = json.dumps(rows, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    result = dict(provenance)
    result["analysis_parameters"] = parameters
    result["result_features"] = len(rows)
    result["result_rows_sha256"] = hashlib.sha256(rows_payload).hexdigest()
    if out.is_file():
        result["geojson_sha256"] = file_sha256(out)
    path = out.with_suffix(out.suffix + ".provenance.json")
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)
    return path


# ---------------------------------------------------------------------------
# Hydrology
# ---------------------------------------------------------------------------
def fill_depressions(z: np.ndarray) -> np.ndarray:
    """Priority-flood depression filling (Barnes et al.). Routing only."""
    ny, nx = z.shape
    filled = np.full_like(z, np.inf)
    closed = np.zeros(z.shape, dtype=bool)
    heap: list[tuple[float, int, int]] = []

    for i in range(ny):
        for j in (0, nx - 1):
            if np.isfinite(z[i, j]):
                heapq.heappush(heap, (z[i, j], i, j)); filled[i, j] = z[i, j]; closed[i, j] = True
    for j in range(nx):
        for i in (0, ny - 1):
            if np.isfinite(z[i, j]) and not closed[i, j]:
                heapq.heappush(heap, (z[i, j], i, j)); filled[i, j] = z[i, j]; closed[i, j] = True

    while heap:
        zc, i, j = heapq.heappop(heap)
        for di, dj, _ in NB:
            ni, nj = i + di, j + dj
            if 0 <= ni < ny and 0 <= nj < nx and not closed[ni, nj] and np.isfinite(z[ni, nj]):
                filled[ni, nj] = max(z[ni, nj], zc)
                closed[ni, nj] = True
                heapq.heappush(heap, (filled[ni, nj], ni, nj))
    filled[~np.isfinite(z)] = np.nan
    return filled


def d8_receivers(filled: np.ndarray, dx: float, dy: float) -> np.ndarray:
    """Index of the steepest-descent neighbour for each cell (-1 if none)."""
    ny, nx = filled.shape
    recv = np.full(ny * nx, -1, dtype=np.int64)
    for i in range(ny):
        for j in range(nx):
            if not np.isfinite(filled[i, j]):
                continue
            best, best_idx = 0.0, -1
            for di, dj, dist in NB:
                ni, nj = i + di, j + dj
                if not (0 <= ni < ny and 0 <= nj < nx) or not np.isfinite(filled[ni, nj]):
                    continue
                length = dist * (dx if dj else dy) if di == 0 or dj == 0 else dist * (dx + dy) / 2
                slope = (filled[i, j] - filled[ni, nj]) / length
                if slope > best:
                    best, best_idx = slope, ni * nx + nj
            recv[i * nx + j] = best_idx
    return recv


def flow_accumulation(filled: np.ndarray, recv: np.ndarray) -> np.ndarray:
    """Number of upstream cells draining through each cell (including itself)."""
    ny, nx = filled.shape
    acc = np.ones(ny * nx, dtype=np.float64)
    acc[~np.isfinite(filled).ravel()] = 0.0
    order = np.argsort(-np.nan_to_num(filled.ravel(), nan=-np.inf))
    for idx in order:
        r = recv[idx]
        if r >= 0:
            acc[r] += acc[idx]
    return acc.reshape(ny, nx)


def upstream_of(idx: int, recv: np.ndarray, children: dict[int, list[int]]) -> np.ndarray:
    """All cell indices draining through idx, by BFS up the receiver tree."""
    out, stack = [idx], [idx]
    while stack:
        c = stack.pop()
        for ch in children.get(c, ()):
            out.append(ch); stack.append(ch)
    return np.array(out, dtype=np.int64)


# ---------------------------------------------------------------------------
def analyse(z, lats, lons, heights, min_accum_cells, top_n,
            min_separation_km=1.5, rank_by="volume"):
    ny, nx = z.shape
    midlat = float(np.mean(lats))
    dy = abs(float(lats[1] - lats[0])) * 110570.0
    dx = abs(float(lons[1] - lons[0])) * 111320.0 * math.cos(math.radians(midlat))
    cell_area = dx * dy
    logger.info("Grid %d x %d, cell %.0f x %.0f m (%.4f km2)", ny, nx, dx, dy, cell_area / 1e6)

    logger.info("Filling depressions...")
    filled = fill_depressions(z)
    logger.info("D8 flow directions...")
    recv = d8_receivers(filled, dx, dy)
    logger.info("Flow accumulation...")
    acc = flow_accumulation(filled, recv)

    children: dict[int, list[int]] = {}
    for idx, r in enumerate(recv):
        if r >= 0:
            children.setdefault(int(r), []).append(idx)

    channel = np.argwhere(acc >= min_accum_cells)
    logger.info("Channel cells (>= %d upstream cells): %d", min_accum_cells, len(channel))
    if len(channel) == 0:
        logger.error("No channel found - lower --min-accum.")
        return []

    zf = z.ravel()
    results = []
    for (i, j) in channel:
        idx = int(i) * nx + int(j)
        sill = z[i, j]
        if not np.isfinite(sill):
            continue
        up = upstream_of(idx, recv, children)
        zup = zf[up]
        good = np.isfinite(zup)
        up, zup = up[good], zup[good]
        if len(up) < 3:
            continue
        for H in heights:
            wl = sill + H
            flooded = zup < wl
            if flooded.sum() < 2:
                continue
            vol = float(np.sum(wl - zup[flooded]) * cell_area)
            area = float(flooded.sum() * cell_area)
            fi = up[flooded] // nx
            fj = up[flooded] % nx
            edge = bool((fi.min() == 0) or (fj.min() == 0) or
                        (fi.max() == ny - 1) or (fj.max() == nx - 1))
            back = float(math.hypot((fi.max() - fi.min()) * dy,
                                    (fj.max() - fj.min()) * dx) / 1000.0)
            results.append({
                "lat": float(lats[i]), "lon": float(lons[j]),
                "sill_m": float(sill), "dam_height_m": float(H),
                "volume_Mm3": vol / 1e6, "area_km2": area / 1e6,
                "backwater_km": back,
                "upstream_cells": int(len(up)),
                "upstream_area_km2": float(len(up) * cell_area / 1e6),
                "truncated": edge,
            })

    # Group by site so each location carries its full height response, not
    # just the biggest number.
    sites: dict[tuple, dict] = {}
    for r in results:
        key = (r["lat"], r["lon"])
        s_ = sites.setdefault(key, {k: v for k, v in r.items()
                                    if k not in ("dam_height_m", "volume_Mm3",
                                                 "area_km2", "backwater_km", "truncated")})
        s_.setdefault("by_height", {})[r["dam_height_m"]] = {
            "volume_Mm3": r["volume_Mm3"], "area_km2": r["area_km2"],
            "backwater_km": r["backwater_km"], "truncated": r["truncated"]}

    hs = sorted(heights)
    for s_ in sites.values():
        bh = s_["by_height"]
        s_["max_volume_Mm3"] = max(v["volume_Mm3"] for v in bh.values())
        # The smallest blockage this site responds to at all. A reach that needs
        # 100 m of debris before it holds anything is a different proposition
        # from one that starts filling at 10 m, and this is the number that
        # separates them.
        lo = min(bh)
        s_["threshold_height_m"] = lo
        s_["efficiency_Mm3_per_m"] = bh[lo]["volume_Mm3"] / lo
        s_["volume_at_min_height_Mm3"] = bh[lo]["volume_Mm3"]
        s_["truncated"] = any(v["truncated"] for v in bh.values())

    return rank_sites(list(sites.values()), rank_by, top_n, min_separation_km)


def rank_sites(sites: list[dict], rank_by: str, top_n: int,
               min_separation_km: float) -> list[dict]:
    """
    Order the sites and suppress near-duplicates.

    Pure: takes dicts, returns dicts, touches no DEM. The ranking is the part
    of this module that makes a claim - "favours sites a small landslide could
    dam" - so it has to be checkable without several thousand DEM queries.
    """
    if rank_by == "efficiency":
        # Threshold height FIRST, then volume per metre within it.
        #
        # Ranking on volume-per-metre alone did the opposite of what it claims.
        # `lo` differs between sites and impounded volume grows faster than
        # linearly with height, so dividing by a LARGER lo can still give a
        # bigger ratio - and a site that holds nothing until 100 m outranks
        # sites that start filling at 10 m. In the committed source-zone run
        # that is exactly what happened: site 2 needs a 100 m blockage and came
        # second of twelve, above four reaches that respond at 10 m.
        #
        # "Favours sites a small landslide could dam" has to mean the small
        # blockage comes first. Efficiency then orders the sites that share a
        # threshold, which is a like-for-like comparison because lo is equal.
        ranked = sorted(sites,
                        key=lambda r: (r["threshold_height_m"],
                                       -r["efficiency_Mm3_per_m"]))
    else:
        ranked = sorted(sites, key=lambda r: -r["max_volume_Mm3"])

    # Non-maximum suppression. Consecutive channel cells describe the SAME pool
    # with the dam nudged a few tens of metres, so without this one valley reach
    # fills the entire table with near-identical rows.
    keep: list[dict] = []
    for r in ranked:
        if any(math.hypot((r["lat"] - k["lat"]) * 110.57,
                          (r["lon"] - k["lon"]) * 111.32 * math.cos(math.radians(r["lat"])))
               < min_separation_km for k in keep):
            continue
        keep.append(r)
        if len(keep) >= top_n:
            break
    return keep


def report(rows, heights, rank_by):
    if not rows:
        print("\nNo impoundment sites found.")
        return
    hs = sorted(heights)
    basis = "volume per metre of blockage" if rank_by == "efficiency" else "maximum impounded volume"
    print("\nRanked by " + basis + "\n")

    head = f"{'#':>3}  {'LAT':>9}{'LON':>10}{'SILL m':>8}{'UPSTR km2':>10}"
    head += "".join(f"{int(h):>8}m" for h in hs)
    head += f"{'THRESH':>8}{'Mm3/m':>9}  FLAG"
    print(head)
    print("-" * len(head))

    for k, r in enumerate(rows, 1):
        line = (f"{k:>3}  {r['lat']:>9.4f}{r['lon']:>10.4f}"
                f"{r['sill_m']:>8.0f}{r['upstream_area_km2']:>10.1f}")
        for h in hs:
            v = r["by_height"].get(h)
            line += f"{v['volume_Mm3']:>9.1f}" if v else f"{'-':>9}"
        line += f"{r['threshold_height_m']:>7.0f}m{r['efficiency_Mm3_per_m']:>9.2f}  "
        line += "truncated" if r["truncated"] else ""
        print(line)

    print("\n  columns are impounded volume (Mm3) at each dam height")
    print("  Mm3/m = volume at the smallest height THIS SITE responds to, per")
    print("  metre of blockage. The threshold differs between sites, so this")
    print("  column compares like for like only within one threshold - which is")
    print("  why the ranking sorts on threshold first.")
    if rank_by == "efficiency":
        by_thr: dict[float, int] = {}
        for r in rows:
            by_thr[r["threshold_height_m"]] = by_thr.get(r["threshold_height_m"], 0) + 1
        summary = ", ".join(f"{n} at {h:g} m" for h, n in sorted(by_thr.items()))
        print(f"\n  threshold blockage: {summary}")
    print("\nScale reference:")
    print("  2024 Thame GLOF, Nepal            ~  2 Mm3")
    print("  2023 South Lhonak GLOF, Sikkim    ~ 50 Mm3  (~180 dead)")

    small = [r for r in rows if r["volume_at_min_height_Mm3"] >= 2]
    if small:
        print(f"\n  {len(small)} site(s) already impound >2 Mm3 with only a {min(hs):g} m blockage.")
        print("  Those are the realistic ones - a 150 m dam is a rare event, 25 m is not.")
    if any(r["truncated"] for r in rows):
        print("  'truncated' = pool reaches the AOI edge; volume is a LOWER bound.")


def write_geojson(rows, out: Path):
    fc = {"type": "FeatureCollection", "features": [
        {"type": "Feature",
         "geometry": {"type": "Point", "coordinates": [r["lon"], r["lat"]]},
         "properties": {k: v for k, v in r.items() if k not in ("lat", "lon")}}
        for r in rows]}
    out.write_text(json.dumps(fc, indent=1))
    logger.info("Wrote %s (%d sites)", out, len(rows))


def main() -> int:
    ap = argparse.ArgumentParser(description="Landslide-dam impoundment susceptibility")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--api", action="store_true", help="sample SRTM via OpenTopoData")
    src.add_argument("--dem", metavar="FILE.tif", help="local DEM GeoTIFF")
    ap.add_argument("--aoi", choices=sorted(AOIS), default="source",
                    help="terrain AOI. 'lhende' is the valley the 26 Aug 2026 avalanche dammed and breached, so its map is a testable prediction")
    ap.add_argument("--cells", type=int, default=12000, help="approx grid cells for --api")
    ap.add_argument("--api-cache", type=Path, default=None,
                    help="save/reuse the exact API elevation grid (.npz); existing cache is validated, not silently refreshed")
    ap.add_argument("--heights", type=float, nargs="+", default=[25, 50, 100, 150])
    ap.add_argument("--min-accum", type=int, default=40,
                    help="upstream cells needed to count as channel")
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--min-separation", type=float, default=1.5,
                    help="km; suppress sites closer than this to a better one")
    ap.add_argument("--rank-by", choices=("volume", "efficiency"), default="efficiency",
                    help="efficiency = volume per metre of blockage (default; "
                         "favours sites a small landslide could dam)")
    ap.add_argument("--geojson", metavar="OUT.geojson")
    args = ap.parse_args()

    if args.cells < 4:
        ap.error("--cells must be at least 4")
    if args.api_cache is not None and not args.api:
        ap.error("--api-cache is only valid together with --api")
    if args.api and args.geojson and args.api_cache is None:
        ap.error("for reproducible GeoJSON output, pass --api-cache PATH so the exact API elevation grid is retained and can be reused")
    if (not args.heights or any(not math.isfinite(h) or h <= 0 for h in args.heights)
            or len(set(args.heights)) != len(args.heights)):
        ap.error("--heights must be unique, finite, positive blockage heights")
    if args.min_accum < 1:
        ap.error("--min-accum must be a positive integer")
    if args.top < 1:
        ap.error("--top must be at least 1")
    if not math.isfinite(args.min_separation) or args.min_separation < 0:
        ap.error("--min-separation must be a finite non-negative number")

    bbox = AOIS[args.aoi]
    print(f"\nAOI {args.aoi}: {bbox}")
    print(f"Dam heights tested: {', '.join(f'{h:g} m' for h in args.heights)}")

    source_file = None
    if args.api:
        z, lats, lons = dem_from_api(bbox, args.cells, args.api_cache)
        source = "OpenTopoData SRTM 30m"
        target_cells = args.cells
    else:
        source_file = Path(args.dem).resolve()
        z, lats, lons = dem_from_file(source_file, bbox=bbox)
        source = "local GeoTIFF"
        target_cells = int(z.size)
    provenance = dem_provenance(z, lats, lons, bbox=bbox, source=source,
                                target_cells=target_cells, source_file=source_file,
                                cache_path=args.api_cache)
    logger.info("DEM_GRID_SHA256=%s shape=%s nodata_cells=%d",
                provenance["grid_sha256"], provenance["shape_rows_cols"],
                provenance["elevation_nodata_cells"])
    logger.info("DEM_PROVENANCE: %s", json.dumps(provenance, sort_keys=True))

    rows = analyse(z, lats, lons, args.heights, args.min_accum, args.top,
                   args.min_separation, args.rank_by)
    report(rows, args.heights, args.rank_by)

    if args.geojson:
        output = Path(args.geojson)
        write_geojson(rows, output)
        sidecar = write_provenance(
            output, provenance,
            {"aoi": args.aoi, "heights_m": args.heights,
             "min_accum_cells": args.min_accum, "top": args.top,
             "min_separation_km": args.min_separation, "rank_by": args.rank_by},
            rows,
        )
        logger.info("DEM_AND_OUTPUT_PROVENANCE=%s", sidecar)

    print("\nCaveats: pools are not tested for spilling over cols, so volumes are")
    print("upper bounds near divides. Dam height is imposed, not predicted -")
    print("couple to a runout model before calling any number a forecast.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
