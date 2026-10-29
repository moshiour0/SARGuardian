#!/usr/bin/env python3
"""Exploratory, target-free NISAR GOFF spatial screening for a custom WKT AOI.

This module ranks repeatable spatial anomalies; it is NOT a hazard detector.
It reads original GOFF HDF5 so SNR and per-pixel radar geometry are not lost.
The original range/azimuth offsets are projected to horizontal East/North using
the product's metadata grid and a user-supplied DEM. It refuses to combine
different track/frame/bandwidth/polarization/CRID groups or mismatched grids.

The anomaly thresholds are adapted from Meena (2026), DOI 10.5281/zenodo.22730419,
but custom-AOI thresholds are not independently calibrated. Correlation >= 0.3
is an additional SARGuardian quality gate and means this is not an exact re-run
of the author's screen. Output status never says "detected hazard".
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.features import geometry_mask, shapes
from rasterio.transform import Affine
from rasterio.warp import Resampling, reproject, transform_geom
from pyproj import Transformer

sys.path.insert(0, str(Path(__file__).resolve().parent))
from aoi_wkt import parse_polygon_wkt  # noqa: E402


SCHEMA = "sarguardian-goff-spatial-screen/1.0"
LAYER_RE = re.compile(r"/pixelOffsets/([^/]+)/layer2$")
GRANULE_RE = re.compile(
    r"NISAR_L2_(?P<processing>PR|UR)_GOFF_(?P<granule>\d{3})_"
    r"(?P<track>\d{3})_(?P<direction>[AD])_(?P<frame>\d{3})_"
    r"(?P<subswath>\d{3})_(?P<bandwidth>\d{4})_(?P<polarization>[A-Z]{2})_"
    r"(?P<reference>\d{8}T\d{6})_\d{8}T\d{6}_"
    r"(?P<secondary>\d{8}T\d{6})_\d{8}T\d{6}_"
    r"(?P<crid>[PX]\d{5})_"
)
DEFAULTS = {
    "snr_min": 5.0,
    "correlation_min": 0.30,
    "magnitude_clip_m": 10.0,
    "local_window_px": 41,
    "direction_window_px": 5,
    "z_min": 4.0,
    "anomaly_min_m": 0.50,
    "directional_coherence_min": 0.75,
    "minimum_pairs": 2,
    "minimum_component_pixels": 10,
}


class ScreenError(RuntimeError):
    """A fail-closed input or scientific-contract error."""


@dataclass(frozen=True)
class GranuleIdentity:
    path: Path
    granule: str
    track: str
    direction: str
    frame: str
    subswath: str
    bandwidth: str
    polarization: str
    reference: str
    secondary: str
    crid: str
    processing: str

    @property
    def series_key(self) -> tuple[str, ...]:
        # Subswath is deliberately retained as metadata but not used as a
        # series key: the author's five consecutive T098 pairs change swath
        # while remaining on the same frame/track and common geocoded grid.
        return (
            self.processing, self.track, self.direction, self.frame,
            self.bandwidth, self.polarization, self.crid,
        )

    @property
    def pair_key(self) -> tuple[str, str]:
        return self.reference, self.secondary

    @property
    def identity_string(self) -> str:
        return (
            f"{self.direction}{self.track} frame {self.frame} "
            f"bandwidth {self.bandwidth} pol {self.polarization} CRID {self.crid}"
        )


@dataclass
class RawPair:
    identity: GranuleIdentity
    east_offsets_m: np.ndarray
    north_offsets_m: np.ndarray
    snr: np.ndarray
    correlation: np.ndarray
    xs: np.ndarray
    ys: np.ndarray
    crs: CRS
    transform: Affine
    grid_mask: np.ndarray
    source_sha256: str
    source_bytes: int
    geometry_note: str


def parse_identity(path: Path) -> GranuleIdentity:
    name = path.name[:-3] if path.name.lower().endswith(".h5") else path.name
    match = GRANULE_RE.search(name)
    if not match:
        raise ScreenError(
            "GRANULE_IDENTITY_UNRESOLVED: expected a full NISAR L2 GOFF granule name "
            "including track, frame, bandwidth, polarization, dates and CRID: " + path.name
        )
    v = match.groupdict()
    return GranuleIdentity(
        path=path, granule=v["granule"], track=v["track"], direction=v["direction"],
        frame=v["frame"], subswath=v["subswath"], bandwidth=v["bandwidth"],
        polarization=v["polarization"], reference=v["reference"], secondary=v["secondary"],
        crid=v["crid"], processing=v["processing"],
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _epsg_from_group(h5: h5py.File, group_path: str) -> int:
    projection_path = group_path + "/projection"
    if projection_path not in h5:
        raise ScreenError("GOFF_PROJECTION_DATASET_MISSING:" + projection_path)
    node = h5[projection_path]
    raw = node.attrs.get("epsg_code")
    if raw is None:
        raw = node.attrs.get("epsg")
    if raw is None:
        try:
            value = np.asarray(node[()]).ravel()
            if value.size:
                raw = value[0]
        except (TypeError, ValueError, OSError):
            raw = None
    try:
        epsg = int(raw)
    except (TypeError, ValueError, OverflowError):
        raise ScreenError("GOFF_EPSG_CODE_MISSING:" + projection_path) from None
    if epsg <= 0:
        raise ScreenError("GOFF_EPSG_CODE_INVALID:" + str(epsg))
    return epsg


def _regular_axis(values: np.ndarray, label: str) -> tuple[np.ndarray, bool, float]:
    axis = np.asarray(values, dtype=np.float64).reshape(-1)
    if axis.size < 2 or not np.isfinite(axis).all():
        raise ScreenError("COORDINATE_AXIS_INVALID:" + label)
    diff = np.diff(axis)
    if np.all(diff > 0):
        ascending = True
    elif np.all(diff < 0):
        ascending = False
    else:
        raise ScreenError("COORDINATE_AXIS_NOT_MONOTONIC:" + label)
    step = float(np.median(np.abs(diff)))
    if not math.isfinite(step) or step <= 0:
        raise ScreenError("COORDINATE_AXIS_STEP_INVALID:" + label)
    tolerance = max(step * 1e-4, 1e-6)
    if float(np.max(np.abs(np.abs(diff) - step))) > tolerance:
        raise ScreenError("COORDINATE_AXIS_NOT_REGULAR:" + label)
    if not ascending:
        axis = axis[::-1]
    return axis, ascending, step


def _normalise_offset_grid(
    arrays: dict[str, np.ndarray], xs: np.ndarray, ys: np.ndarray
) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray, Affine]:
    xasc, x_was_ascending, dx = _regular_axis(xs, "GOFF_X")
    yasc, y_was_ascending, dy = _regular_axis(ys, "GOFF_Y")
    normalised: dict[str, np.ndarray] = {}
    for name, values in arrays.items():
        if values.ndim != 2 or values.shape != (len(ys), len(xs)):
            raise ScreenError(
                f"GOFF_ARRAY_SHAPE_MISMATCH:{name}:{values.shape}:expected={(len(ys), len(xs))}"
            )
        a = np.asarray(values)
        if not x_was_ascending:
            a = a[:, ::-1]
        # Raster row 0 is north. Normalize all arrays to descending northing.
        if y_was_ascending:
            a = a[::-1, :]
        normalised[name] = a
    xnorm = xasc
    ynorm = yasc[::-1]
    transform = Affine(dx, 0.0, float(xnorm[0] - dx / 2.0),
                       0.0, -dy, float(ynorm[0] + dy / 2.0))
    return normalised, xnorm, ynorm, transform


def _normalise_geometry_axes(
    hgt: np.ndarray,
    ys: np.ndarray,
    xs: np.ndarray,
    surfaces: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    h = np.asarray(hgt, dtype=np.float64).reshape(-1)
    y = np.asarray(ys, dtype=np.float64).reshape(-1)
    x = np.asarray(xs, dtype=np.float64).reshape(-1)
    if h.size < 2 or not np.isfinite(h).all():
        raise ScreenError("RADAR_GRID_HEIGHT_AXIS_INVALID")
    if np.all(np.diff(h) < 0):
        h = h[::-1]
        surfaces = {k: v[::-1, :, :] for k, v in surfaces.items()}
    elif not np.all(np.diff(h) > 0):
        raise ScreenError("RADAR_GRID_HEIGHT_AXIS_NOT_MONOTONIC")
    if np.all(np.diff(y) < 0):
        y = y[::-1]
        surfaces = {k: v[:, ::-1, :] for k, v in surfaces.items()}
    elif not np.all(np.diff(y) > 0):
        raise ScreenError("RADAR_GRID_Y_AXIS_NOT_MONOTONIC")
    if np.all(np.diff(x) < 0):
        x = x[::-1]
        surfaces = {k: v[:, :, ::-1] for k, v in surfaces.items()}
    elif not np.all(np.diff(x) > 0):
        raise ScreenError("RADAR_GRID_X_AXIS_NOT_MONOTONIC")
    expected = (h.size, y.size, x.size)
    for name, surface in surfaces.items():
        if surface.shape != expected:
            raise ScreenError(
                f"RADAR_GRID_SURFACE_SHAPE_MISMATCH:{name}:{surface.shape}:expected={expected}"
            )
    return h, y, x, surfaces


def interpolate_regular_volume(
    volume: np.ndarray,
    heights_axis: np.ndarray,
    y_axis_ascending: np.ndarray,
    x_axis_ascending: np.ndarray,
    elevation: np.ndarray,
    x_centres: np.ndarray,
    y_centres_descending: np.ndarray,
    row_start: int = 0,
) -> np.ndarray:
    """Vectorized trilinear interpolation on one row block, with no SciPy dependency."""
    h, y, x = heights_axis, y_axis_ascending, x_axis_ascending
    rows = elevation.shape[0]
    X = np.broadcast_to(x_centres[None, :], (rows, x_centres.size))
    Y = np.broadcast_to(y_centres_descending[row_start:row_start + rows, None],
                        (rows, x_centres.size))
    Z = np.asarray(elevation, dtype=np.float64)
    # Never extrapolate geometry at image edges: missing/out-of-range
    # radarGrid support becomes NoData, not copied edge coefficients.
    good = (
        np.isfinite(Z) & np.isfinite(X) & np.isfinite(Y)
        & (Z >= h[0]) & (Z <= h[-1])
        & (Y >= y[0]) & (Y <= y[-1])
        & (X >= x[0]) & (X <= x[-1])
    )
    out = np.full(Z.shape, np.nan, dtype=np.float32)
    if not good.any():
        return out

    xv, yv, zv = X[good], Y[good], Z[good]
    ix1 = np.searchsorted(x, xv, side="right").clip(1, x.size - 1)
    iy1 = np.searchsorted(y, yv, side="right").clip(1, y.size - 1)
    iz1 = np.searchsorted(h, zv, side="right").clip(1, h.size - 1)
    ix0, iy0, iz0 = ix1 - 1, iy1 - 1, iz1 - 1
    tx = np.clip((xv - x[ix0]) / (x[ix1] - x[ix0]), 0.0, 1.0)
    ty = np.clip((yv - y[iy0]) / (y[iy1] - y[iy0]), 0.0, 1.0)
    tz = np.clip((zv - h[iz0]) / (h[iz1] - h[iz0]), 0.0, 1.0)

    v000 = volume[iz0, iy0, ix0]
    v001 = volume[iz0, iy0, ix1]
    v010 = volume[iz0, iy1, ix0]
    v011 = volume[iz0, iy1, ix1]
    v100 = volume[iz1, iy0, ix0]
    v101 = volume[iz1, iy0, ix1]
    v110 = volume[iz1, iy1, ix0]
    v111 = volume[iz1, iy1, ix1]
    low = (v000 * (1 - tx) + v001 * tx) * (1 - ty) + (v010 * (1 - tx) + v011 * tx) * ty
    high = (v100 * (1 - tx) + v101 * tx) * (1 - ty) + (v110 * (1 - tx) + v111 * tx) * ty
    values = low * (1 - tz) + high * tz
    finite = np.isfinite(values)
    flat = np.flatnonzero(good)
    out.ravel()[flat[finite]] = values[finite].astype(np.float32)
    return out


def read_raw_pair(
    path: Path,
    dem_path: Path,
    ring: list[tuple[float, float]],
    polarization: str = "HH",
    frequency: str = "frequencyA",
    layer: str = "layer2",
    snr_min: float = 5.0,
    correlation_min: float = 0.30,
) -> RawPair:
    identity = parse_identity(path)
    if not path.is_file():
        raise ScreenError("INPUT_HDF5_NOT_FOUND:" + str(path))
    if not dem_path.is_file():
        raise ScreenError("DEM_NOT_FOUND:" + str(dem_path))
    expected_group = f"science/LSAR/GOFF/grids/{frequency}/pixelOffsets/{polarization}/{layer}"
    with h5py.File(path, "r") as h5:
        if expected_group not in h5:
            # Search exact layer2 groups in the product so group-layout changes
            # fail with a diagnostic instead of silently choosing an unrelated layer.
            alternatives = []
            prefix = f"science/LSAR/GOFF/grids/"
            def visitor(name: str, node: Any) -> None:
                if isinstance(node, h5py.Group) and name.endswith(f"/pixelOffsets/{polarization}/{layer}"):
                    alternatives.append(name)
            h5.visititems(visitor)
            if len(alternatives) != 1:
                raise ScreenError(
                    f"GOFF_LAYER_NOT_FOUND_OR_AMBIGUOUS:{expected_group}:alternatives={alternatives}"
                )
            expected_group = alternatives[0]

        def require(name: str) -> np.ndarray:
            key = expected_group + "/" + name
            if key not in h5:
                raise ScreenError("GOFF_REQUIRED_DATASET_MISSING:" + key)
            return np.asarray(h5[key][()], dtype=np.float32)

        at = require("alongTrackOffset")
        sr = require("slantRangeOffset")
        snr = require("snr")
        corr = require("correlationSurfacePeak")
        xs = np.asarray(h5[expected_group + "/xCoordinates"][()], dtype=np.float64)
        ys = np.asarray(h5[expected_group + "/yCoordinates"][()], dtype=np.float64)
        epsg = _epsg_from_group(h5, expected_group)
        geometry_group = "science/LSAR/GOFF/metadata/radarGrid"
        if geometry_group not in h5:
            raise ScreenError("RADAR_GRID_METADATA_MISSING")
        radar = h5[geometry_group]
        required_geometry = (
            "heightAboveEllipsoid", "xCoordinates", "yCoordinates", "incidenceAngle",
            "losUnitVectorX", "losUnitVectorY", "alongTrackUnitVectorX", "alongTrackUnitVectorY",
        )
        missing = [name for name in required_geometry if name not in radar]
        if missing:
            raise ScreenError("RADAR_GRID_FIELDS_MISSING:" + ",".join(missing))
        heights = np.asarray(radar["heightAboveEllipsoid"][()], dtype=np.float64)
        gx = np.asarray(radar["xCoordinates"][()], dtype=np.float64)
        gy = np.asarray(radar["yCoordinates"][()], dtype=np.float64)
        surfaces = {
            "incidence": np.asarray(radar["incidenceAngle"][()], dtype=np.float32),
            "los_x": np.asarray(radar["losUnitVectorX"][()], dtype=np.float32),
            "los_y": np.asarray(radar["losUnitVectorY"][()], dtype=np.float32),
            "along_x": np.asarray(radar["alongTrackUnitVectorX"][()], dtype=np.float32),
            "along_y": np.asarray(radar["alongTrackUnitVectorY"][()], dtype=np.float32),
        }

    if not (at.shape == sr.shape == snr.shape == corr.shape):
        raise ScreenError("GOFF_OFFSET_AND_QUALITY_SHAPES_DISAGREE")
    arrays, x, y, transform = _normalise_offset_grid(
        {"along": at, "range": sr, "snr": snr, "correlation": corr}, xs, ys
    )
    h_axis, gy_axis, gx_axis, surfaces = _normalise_geometry_axes(heights, gy, gx, surfaces)
    crs = CRS.from_epsg(epsg)
    if crs.is_geographic:
        raise ScreenError("GOFF_OFFSET_GRID_EXPECTS_PROJECTED_CRS_NOT_GEOGRAPHIC")
    polygon = {"type": "Polygon", "coordinates": [[list(point) for point in ring]]}
    projected_polygon = transform_geom("EPSG:4326", crs, polygon)
    aoi_mask = geometry_mask([projected_polygon], out_shape=arrays["range"].shape,
                             transform=transform, invert=True, all_touched=False)
    if not aoi_mask.any():
        raise ScreenError("AOI_HAS_NO_PIXEL_CENTRES_ON_GOFF_GRID")

    # Reproject the supplied DEM onto the native GOFF offset grid. Pixels with
    # missing DEM stay invalid; no mountain elevation is silently imputed.
    elevation = np.full(arrays["range"].shape, np.nan, dtype=np.float32)
    with rasterio.open(dem_path) as dem:
        if dem.crs is None:
            raise ScreenError("DEM_CRS_MISSING")
        if dem.count < 1:
            raise ScreenError("DEM_HAS_NO_BANDS")
        reproject(
            source=rasterio.band(dem, 1), destination=elevation,
            src_transform=dem.transform, src_crs=dem.crs, src_nodata=dem.nodata,
            dst_transform=transform, dst_crs=crs, dst_nodata=np.nan,
            resampling=Resampling.bilinear,
        )
    dem_valid = np.isfinite(elevation)
    if int((dem_valid & aoi_mask).sum()) < max(10, int(0.5 * aoi_mask.sum())):
        raise ScreenError("DEM_VALID_AOI_COVERAGE_BELOW_50_PERCENT")

    h, w = arrays["range"].shape
    east = np.full((h, w), np.nan, dtype=np.float32)
    north = np.full((h, w), np.nan, dtype=np.float32)
    block_rows = 96
    for r0 in range(0, h, block_rows):
        r1 = min(h, r0 + block_rows)
        elev_block = elevation[r0:r1]
        inc = interpolate_regular_volume(surfaces["incidence"], h_axis, gy_axis, gx_axis,
                                         elev_block, x, y, row_start=r0)
        lx = interpolate_regular_volume(surfaces["los_x"], h_axis, gy_axis, gx_axis,
                                        elev_block, x, y, row_start=r0)
        ly = interpolate_regular_volume(surfaces["los_y"], h_axis, gy_axis, gx_axis,
                                        elev_block, x, y, row_start=r0)
        ax = interpolate_regular_volume(surfaces["along_x"], h_axis, gy_axis, gx_axis,
                                        elev_block, x, y, row_start=r0)
        ay = interpolate_regular_volume(surfaces["along_y"], h_axis, gy_axis, gx_axis,
                                        elev_block, x, y, row_start=r0)
        inc_rad = np.deg2rad(inc)
        horizontal_los = np.hypot(lx, ly)
        good_geometry = (
            np.isfinite(inc_rad) & (inc_rad > 0) & (inc_rad < (np.pi / 2))
            & np.isfinite(horizontal_los) & (horizontal_los > 1e-6)
            & np.isfinite(ax) & np.isfinite(ay)
        )
        kx = np.full(inc.shape, np.nan, dtype=np.float32)
        ky = np.full(inc.shape, np.nan, dtype=np.float32)
        kx[good_geometry] = (-(lx[good_geometry] / horizontal_los[good_geometry])
                             / np.sin(inc_rad[good_geometry]))
        ky[good_geometry] = (-(ly[good_geometry] / horizontal_los[good_geometry])
                             / np.sin(inc_rad[good_geometry]))
        block_good = good_geometry & dem_valid[r0:r1]
        east[r0:r1] = np.where(
            block_good, arrays["along"][r0:r1] * ax + arrays["range"][r0:r1] * kx, np.nan
        )
        north[r0:r1] = np.where(
            block_good, arrays["along"][r0:r1] * ay + arrays["range"][r0:r1] * ky, np.nan
        )

    return RawPair(
        identity=identity, east_offsets_m=east, north_offsets_m=north,
        snr=arrays["snr"], correlation=arrays["correlation"],
        xs=x, ys=y, crs=crs, transform=transform, grid_mask=aoi_mask,
        source_sha256=sha256_file(path), source_bytes=path.stat().st_size,
        geometry_note=(
            "Per-pixel projection interpolated from this product's radarGrid at the supplied DEM elevations. "
            "DEM vertical datum compatibility with heightAboveEllipsoid must be verified by the operator."
        ),
    )


def box_sum(values: np.ndarray, window: int) -> np.ndarray:
    """NaN-free square box sum, implemented with float64 integral sums."""
    if window < 1 or window % 2 != 1:
        raise ValueError("WINDOW_MUST_BE_POSITIVE_ODD")
    radius = window // 2
    padded = np.pad(np.asarray(values, dtype=np.float64),
                    ((radius, radius), (radius, radius)), mode="constant")
    # First pass uses explicit leading-zero cumulative sums, so the output
    # shape is exactly HxW even for window=1 and at image boundaries.
    c0 = np.cumsum(padded, axis=0, dtype=np.float64)
    vertical = np.empty((values.shape[0], padded.shape[1]), dtype=np.float64)
    vertical[0] = c0[window - 1]
    if values.shape[0] > 1:
        vertical[1:] = c0[window:] - c0[:-window]
    del c0, padded
    c1 = np.cumsum(vertical, axis=1, dtype=np.float64)
    result = np.empty(values.shape, dtype=np.float32)
    result[:, 0] = c1[:, window - 1]
    if values.shape[1] > 1:
        result[:, 1:] = c1[:, window:] - c1[:, :-window]
    return result


def local_mean(
    values: np.ndarray, valid: np.ndarray, window: int, minimum_fraction: float = 0.15
) -> np.ndarray:
    m = np.asarray(valid, dtype=bool) & np.isfinite(values)
    weights = m.astype(np.float32)
    count = box_sum(weights, window)
    total = box_sum(np.where(m, values, 0.0).astype(np.float32), window)
    out = np.full(values.shape, np.nan, dtype=np.float32)
    enough = count >= (window * window * minimum_fraction)
    out[enough] = total[enough] / count[enough]
    return out


def directional_coherence(
    residual_east: np.ndarray,
    residual_north: np.ndarray,
    valid: np.ndarray,
    window: int,
) -> np.ndarray:
    magnitude = np.hypot(residual_east, residual_north)
    ok = valid & np.isfinite(magnitude) & (magnitude > 1e-6)
    ux = np.full(magnitude.shape, np.nan, dtype=np.float32)
    uy = np.full(magnitude.shape, np.nan, dtype=np.float32)
    ux[ok] = residual_east[ok] / magnitude[ok]
    uy[ok] = residual_north[ok] / magnitude[ok]
    mean_x = local_mean(ux, ok, window, minimum_fraction=0.51)
    mean_y = local_mean(uy, ok, window, minimum_fraction=0.51)
    return np.minimum(1.0, np.hypot(mean_x, mean_y)).astype(np.float32)


def score_pair_arrays(
    east: np.ndarray,
    north: np.ndarray,
    snr: np.ndarray,
    correlation: np.ndarray,
    aoi_mask: np.ndarray,
    snr_min: float = DEFAULTS["snr_min"],
    correlation_min: float = DEFAULTS["correlation_min"],
    magnitude_clip_m: float = DEFAULTS["magnitude_clip_m"],
    local_window_px: int = DEFAULTS["local_window_px"],
    direction_window_px: int = DEFAULTS["direction_window_px"],
    z_min: float = DEFAULTS["z_min"],
    anomaly_min_m: float = DEFAULTS["anomaly_min_m"],
    coherence_min: float = DEFAULTS["directional_coherence_min"],
) -> dict[str, Any]:
    """Score one pair on the whole valid scene, then restrict candidates to the AOI."""
    arrays = (east, north, snr, correlation, aoi_mask)
    if any(np.asarray(a).shape != np.asarray(east).shape for a in arrays):
        raise ScreenError("PAIR_ARRAY_SHAPES_DISAGREE")
    if local_window_px < 3 or local_window_px % 2 != 1:
        raise ScreenError("LOCAL_WINDOW_MUST_BE_ODD")
    if direction_window_px < 3 or direction_window_px % 2 != 1:
        raise ScreenError("DIRECTION_WINDOW_MUST_BE_ODD")
    magnitude = np.hypot(east, north).astype(np.float32)
    valid = (
        np.isfinite(east) & np.isfinite(north) & np.isfinite(snr) & (snr >= snr_min)
        & np.isfinite(correlation) & (correlation >= correlation_min)
        & np.isfinite(magnitude) & (magnitude <= magnitude_clip_m)
    )
    local_mag = local_mean(magnitude, valid, local_window_px)
    anomaly = magnitude - local_mag
    abs_dev = np.abs(anomaly)
    scale = local_mean(abs_dev, valid & np.isfinite(anomaly), local_window_px)
    z_score = anomaly / np.maximum(scale, 0.05)
    local_east = local_mean(east, valid, local_window_px)
    local_north = local_mean(north, valid, local_window_px)
    residual_east = east - local_east
    residual_north = north - local_north
    coherence = directional_coherence(residual_east, residual_north, valid, direction_window_px)
    hit = (
        valid & np.asarray(aoi_mask, dtype=bool)
        & np.isfinite(z_score) & (z_score >= z_min)
        & np.isfinite(anomaly) & (anomaly >= anomaly_min_m)
        & np.isfinite(coherence) & (coherence >= coherence_min)
    )
    return {
        "valid": valid,
        "candidate": hit,
        "magnitude_m": magnitude,
        "anomaly_m": anomaly,
        "z_score": z_score,
        "directional_coherence": coherence,
        "valid_aoi_pixels": int((valid & aoi_mask).sum()),
        "candidate_pixels": int(hit.sum()),
        "aoi_grid_pixels": int(np.asarray(aoi_mask, dtype=bool).sum()),
        "candidate_aoi_fraction": float(hit.sum() / max(1, np.asarray(aoi_mask, dtype=bool).sum())),
    }


def _grid_compatible(first: RawPair, other: RawPair) -> bool:
    return (
        first.crs == other.crs
        and first.east_offsets_m.shape == other.east_offsets_m.shape
        and np.allclose(first.xs, other.xs, rtol=0, atol=1e-3)
        and np.allclose(first.ys, other.ys, rtol=0, atol=1e-3)
        and np.allclose(tuple(first.transform), tuple(other.transform), rtol=0, atol=1e-6)
        and np.array_equal(first.grid_mask, other.grid_mask)
    )


def validate_series(pairs: list[RawPair]) -> None:
    if not pairs:
        raise ScreenError("NO_INPUT_PAIRS")
    key = pairs[0].identity.series_key
    for pair in pairs:
        if pair.identity.series_key != key:
            raise ScreenError(
                "MIXED_SERIES_IDENTITY: refusing to merge tracks, directions, frames, "
                "bandwidths, polarizations, processing types or CRIDs; got "
                + pair.identity.identity_string
            )
    if len({p.identity.pair_key for p in pairs}) != len(pairs):
        raise ScreenError("DUPLICATE_ACQUISITION_INTERVALS")
    reference = pairs[0]
    for pair in pairs[1:]:
        if not _grid_compatible(reference, pair):
            raise ScreenError("INCOMPATIBLE_PIXEL_GRID: no reprojection/resampling is done implicitly")
    for pair in pairs:
        if pair.identity.secondary <= pair.identity.reference:
            raise ScreenError("INVALID_ACQUISITION_INTERVAL:" + pair.identity.pair_key.__repr__())
    pairs.sort(key=lambda p: (p.identity.secondary, p.identity.reference))


def _component_rows(mask: np.ndarray, transform: Affine, crs: CRS, minimum_pixels: int) -> list[dict[str, Any]]:
    pixel_area = abs(transform.a * transform.e)
    rows: list[dict[str, Any]] = []
    to_wgs84 = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    next_label = 1
    # Use rasterio polygonization with 8-connectivity and count pixels from the
    # polygon area. Each output feature is a review polygon, not a landslide.
    for geometry, value in shapes(mask.astype(np.uint8), mask=mask, transform=transform, connectivity=8):
        if int(value) != 1:
            continue
        area_m2 = None
        try:
            # Pixel-aligned polygon area; all supported grids must be north-up,
            # projected and regular. Geometry coordinates are in CRS units.
            coords = geometry["coordinates"][0]
            twice = 0.0
            for i in range(len(coords) - 1):
                twice += coords[i][0] * coords[i + 1][1] - coords[i + 1][0] * coords[i][1]
            area_m2 = abs(twice) / 2.0
        except (KeyError, IndexError, TypeError):
            area_m2 = 0.0
        pixels = int(round((area_m2 or 0.0) / pixel_area)) if pixel_area > 0 else 0
        if pixels < minimum_pixels:
            continue
        bounds = _geometry_bounds(geometry)
        cx, cy = (bounds[0] + bounds[2]) / 2.0, (bounds[1] + bounds[3]) / 2.0
        lon, lat = to_wgs84.transform(cx, cy)
        rows.append({
            "candidate_id": f"C{next_label:04d}", "pixels_approx": pixels,
            "area_km2_approx": round((area_m2 or 0.0) / 1e6, 6),
            "bbox_center_lon": round(lon, 6), "bbox_center_lat": round(lat, 6),
            "geometry": geometry,
        })
        next_label += 1
    rows.sort(key=lambda r: (-r["pixels_approx"], r["bbox_center_lat"], r["bbox_center_lon"]))
    for i, row in enumerate(rows, 1):
        row["candidate_id"] = f"C{i:04d}"
    return rows


def _geometry_bounds(geometry: dict[str, Any]) -> tuple[float, float, float, float]:
    coordinates = geometry.get("coordinates", [])
    if geometry.get("type") == "Polygon":
        points = [point for ring in coordinates for point in ring]
    elif geometry.get("type") == "MultiPolygon":
        points = [point for polygon in coordinates for ring in polygon for point in ring]
    else:
        points = []
    if not points:
        return (0.0, 0.0, 0.0, 0.0)
    xs = [float(p[0]) for p in points]
    ys = [float(p[1]) for p in points]
    return (min(xs), min(ys), max(xs), max(ys))


def _write_raster(path: Path, data: np.ndarray, transform: Affine, crs: CRS,
                  nodata: float, description: str) -> None:
    arr = np.asarray(data)
    with rasterio.open(
        path, "w", driver="GTiff", width=arr.shape[1], height=arr.shape[0],
        count=1, dtype=arr.dtype, crs=crs, transform=transform, nodata=nodata,
        compress="deflate", tiled=True, blockxsize=256, blockysize=256,
        predictor=3 if arr.dtype.kind == "f" else 2,
    ) as dst:
        dst.write(arr, 1)
        dst.set_band_description(1, description)
        dst.update_tags(
            schema=SCHEMA, role="EXPLORATORY_RANKING_ONLY",
            no_hazard_alert="true",
            projection="Per-pixel horizontal East/North from GOFF range/azimuth offsets and radarGrid geometry",
        )


def run_screen(
    inputs: list[Path],
    dem_path: Path,
    aoi_wkt: str,
    out_dir: Path,
    *,
    polarization: str = "HH",
    frequency: str = "frequencyA",
    layer: str = "layer2",
    snr_min: float = DEFAULTS["snr_min"],
    correlation_min: float = DEFAULTS["correlation_min"],
    local_window_px: int = DEFAULTS["local_window_px"],
    direction_window_px: int = DEFAULTS["direction_window_px"],
    z_min: float = DEFAULTS["z_min"],
    anomaly_min_m: float = DEFAULTS["anomaly_min_m"],
    coherence_min: float = DEFAULTS["directional_coherence_min"],
    minimum_pairs: int = DEFAULTS["minimum_pairs"],
    minimum_component_pixels: int = DEFAULTS["minimum_component_pixels"],
) -> dict[str, Any]:
    if minimum_pairs < 2:
        raise ScreenError("MINIMUM_PAIRS_MUST_BE_AT_LEAST_TWO")
    if minimum_component_pixels < 1:
        raise ScreenError("MINIMUM_COMPONENT_PIXELS_MUST_BE_POSITIVE")
    if len(inputs) != len(set(p.resolve() for p in inputs)):
        raise ScreenError("DUPLICATE_INPUT_PATHS")
    ring = parse_polygon_wkt(aoi_wkt)
    pairs = [
        read_raw_pair(path, dem_path, ring, polarization=polarization,
                      frequency=frequency, layer=layer, snr_min=snr_min,
                      correlation_min=correlation_min)
        for path in inputs
    ]
    validate_series(pairs)
    out_dir.mkdir(parents=True, exist_ok=True)
    first = pairs[0]
    hit_count = np.zeros(first.east_offsets_m.shape, dtype=np.uint8)
    anomaly_sum = np.zeros(first.east_offsets_m.shape, dtype=np.float32)
    coherence_sum = np.zeros(first.east_offsets_m.shape, dtype=np.float32)
    reports = []

    for pair in pairs:
        result = score_pair_arrays(
            pair.east_offsets_m, pair.north_offsets_m, pair.snr, pair.correlation,
            pair.grid_mask, snr_min=snr_min, correlation_min=correlation_min,
            local_window_px=local_window_px, direction_window_px=direction_window_px,
            z_min=z_min, anomaly_min_m=anomaly_min_m, coherence_min=coherence_min,
        )
        hit = result["candidate"]
        hit_count += hit.astype(np.uint8)
        anomaly_sum += np.where(hit, result["anomaly_m"], 0.0).astype(np.float32)
        coherence_sum += np.where(hit, result["directional_coherence"], 0.0).astype(np.float32)
        reports.append({
            "granule": pair.identity.path.name,
            "granule_identity": pair.identity.identity_string,
            "subswath": pair.identity.subswath,
            "reference": pair.identity.reference,
            "secondary": pair.identity.secondary,
            "interval_days": (
                datetime.strptime(pair.identity.secondary, "%Y%m%dT%H%M%S")
                - datetime.strptime(pair.identity.reference, "%Y%m%dT%H%M%S")
            ).total_seconds() / 86400.0,
            "source_sha256": pair.source_sha256,
            "source_bytes": pair.source_bytes,
            "aoi_grid_pixels": result["aoi_grid_pixels"],
            "valid_aoi_pixels": result["valid_aoi_pixels"],
            "candidate_pixels": result["candidate_pixels"],
            "candidate_aoi_fraction": result["candidate_aoi_fraction"],
            "status": "SCREENED_NOT_DETECTION",
        })

    persistent = hit_count >= minimum_pairs
    anomaly_mean = np.full(hit_count.shape, np.nan, dtype=np.float32)
    coherence_mean = np.full(hit_count.shape, np.nan, dtype=np.float32)
    flagged = hit_count > 0
    anomaly_mean[flagged] = anomaly_sum[flagged] / hit_count[flagged]
    coherence_mean[flagged] = coherence_sum[flagged] / hit_count[flagged]
    # Retain every pair's count for audit, but the candidate mask is gated by
    # persistence and minimum connected support. It is not a hazard label.
    components = _component_rows(persistent, first.transform, first.crs, minimum_component_pixels)
    qualified_mask = np.zeros(persistent.shape, dtype=bool)
    if components:
        polygon_mask = np.zeros(persistent.shape, dtype=np.uint8)
        # Reconstruct the retained component mask from pixel counts / polygon
        # shape without using bbox: rasterize actual polygons back to the grid.
        from rasterio.features import rasterize
        burn = [(row["geometry"], 1) for row in components]
        rasterize(burn, out=polygon_mask, transform=first.transform, fill=0, dtype=np.uint8)
        qualified_mask = polygon_mask.astype(bool)

    status = (
        "EXPLORATORY_RANKING_ONLY"
        if len(pairs) >= minimum_pairs
        else "NO_TEMPORAL_PERSISTENCE_TEST"
    )
    reasons = [
        "Candidate thresholds are adapted from a single published NISAR GOFF case study; not calibrated for arbitrary AOIs.",
        "Correlation>=0.3 is an additional SARGuardian gate; results are not an exact reproduction of the published screen.",
        "A candidate is a repeatable coherent offset anomaly for manual review, not proof of landslide or glacier-collapse motion.",
        "No slope-direction validation is performed without a separate terrain-derived downslope/aspect assessment.",
        "DEM vertical datum must be checked against radarGrid heightAboveEllipsoid before scientific magnitude interpretation.",
        "No spatial-null/multiple-testing calibration is applied in this prototype; component ranks are exploratory only.",
    ]
    if len(pairs) < minimum_pairs:
        reasons.append(
            f"Only {len(pairs)} compatible pair(s) were supplied; at least {minimum_pairs} are required "
            "for even an exploratory persistence screen."
        )

    hit_count_path = out_dir / "pair_persistence_count.tif"
    anomaly_path = out_dir / "mean_anomaly_m.tif"
    coherence_path = out_dir / "mean_directional_coherence.tif"
    candidates_path = out_dir / "persistent_candidate_mask.tif"
    _write_raster(hit_count_path, hit_count, first.transform, first.crs, 255,
                  "Number of GOFF pairs with a coherent local anomaly")
    _write_raster(anomaly_path, anomaly_mean, first.transform, first.crs, np.nan,
                  "Mean positive local displacement-magnitude anomaly among flagged pairs (m)")
    _write_raster(coherence_path, coherence_mean, first.transform, first.crs, np.nan,
                  "Mean residual directional coherence among flagged pairs")
    _write_raster(candidates_path, qualified_mask.astype(np.uint8), first.transform, first.crs, 255,
                  "Exploratory persistent-component mask; 1=review candidate, not hazard")

    candidate_features = []
    for row in components:
        # RFC 7946 GeoJSON coordinates are WGS84 longitude/latitude. Rasters
        # stay in the native projected CRS; polygon coordinates are transformed.
        geometry_wgs84 = transform_geom(first.crs, "EPSG:4326", row["geometry"])
        feature = {
            "type": "Feature",
            "properties": {k: v for k, v in row.items() if k != "geometry"},
            "geometry": geometry_wgs84,
        }
        candidate_features.append(feature)
    geojson_path = out_dir / "candidates.geojson"
    geojson_path.write_text(json.dumps({
        "type": "FeatureCollection", "features": candidate_features,
        "properties": {"status": status, "interpretation": "Exploratory review candidates; not hazard detections."},
    }, indent=2, ensure_ascii=False) + "\\n", encoding="utf-8")
    csv_path = out_dir / "candidates.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=["candidate_id", "pixels_approx", "area_km2_approx", "bbox_center_lon", "bbox_center_lat"],
        )
        writer.writeheader()
        for row in components:
            writer.writerow({k: row[k] for k in writer.fieldnames})

    summary = {
        "schema_version": SCHEMA,
        "status": status,
        "interpretation": "Exploratory target-free ranking only. No hazard alert, failure prediction or validated detection is issued.",
        "series_identity": first.identity.identity_string,
        "processing": first.identity.processing,
        "crid": first.identity.crid,
        "pair_count": len(pairs),
        "minimum_pairs_for_persistence": minimum_pairs,
        "aoi_wkt": aoi_wkt,
        "aoi_grid_pixels": int(first.grid_mask.sum()),
        "grid_shape": list(first.grid_mask.shape),
        "grid_crs": first.crs.to_string(),
        "grid_transform": list(first.transform),
        "dem": {"path": str(dem_path), "sha256": sha256_file(dem_path),
                "vertical_datum_status": "MUST_BE_VERIFIED_COMPATIBLE_WITH_HEIGHT_ABOVE_ELLIPSOID"},
        "parameters": {
            "snr_min": snr_min, "correlation_min": correlation_min,
            "magnitude_clip_m": DEFAULTS["magnitude_clip_m"],
            "local_window_px": local_window_px, "direction_window_px": direction_window_px,
            "z_min": z_min, "anomaly_min_m": anomaly_min_m,
            "directional_coherence_min": coherence_min,
            "minimum_pairs": minimum_pairs, "minimum_component_pixels": minimum_component_pixels,
        },
        "parameter_source": {
            "adapted_from": "Meena (2026), NISAR GOFF source screening, Zenodo 10.5281/zenodo.22730419",
            "license": "CC BY 4.0 for code/data outside separately licensed extent_CC-BY-NC-4.0 folder",
            "custom_aoi_calibration": "NOT_PERFORMED",
            "exact_paper_reproduction": False,
        },
        "source_pairs": reports,
        "persistent_pixels_before_component_size_filter": int(persistent.sum()),
        "review_candidate_components": len(components),
        "qualified_candidate_pixels": int(qualified_mask.sum()),
        "reasons_and_limitations": reasons,
        "outputs": {
            "pair_persistence_count": hit_count_path.name,
            "mean_anomaly": anomaly_path.name,
            "mean_directional_coherence": coherence_path.name,
            "persistent_candidate_mask": candidates_path.name,
            "candidates_geojson": geojson_path.name,
            "candidates_csv": csv_path.name,
        },
    }
    summary["output_sha256"] = {
        name: sha256_file(out_dir / filename) for name, filename in summary["outputs"].items()
    }
    (out_dir / "screen-summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + "\\n", encoding="utf-8"
    )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--read", nargs="+", required=True, type=Path,
                        help="Raw NISAR L2 GOFF HDF5 products from one compatible series")
    parser.add_argument("--dem", required=True, type=Path,
                        help="Georeferenced DEM raster with a verified vertical datum")
    parser.add_argument("--aoi-wkt", required=True,
                        help="One WGS84 POLYGON WKT. Never defaults to the AOI centre.")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--polarization", default="HH")
    parser.add_argument("--frequency", default="frequencyA")
    parser.add_argument("--layer", default="layer2")
    parser.add_argument("--snr-min", type=float, default=DEFAULTS["snr_min"])
    parser.add_argument("--correlation-min", type=float, default=DEFAULTS["correlation_min"])
    parser.add_argument("--local-window-px", type=int, default=DEFAULTS["local_window_px"])
    parser.add_argument("--direction-window-px", type=int, default=DEFAULTS["direction_window_px"])
    parser.add_argument("--z-min", type=float, default=DEFAULTS["z_min"])
    parser.add_argument("--anomaly-min-m", type=float, default=DEFAULTS["anomaly_min_m"])
    parser.add_argument("--coherence-min", type=float, default=DEFAULTS["directional_coherence_min"])
    parser.add_argument("--minimum-pairs", type=int, default=DEFAULTS["minimum_pairs"])
    parser.add_argument("--minimum-component-pixels", type=int, default=DEFAULTS["minimum_component_pixels"])
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = run_screen(
            args.read, args.dem, args.aoi_wkt, args.out,
            polarization=args.polarization, frequency=args.frequency, layer=args.layer,
            snr_min=args.snr_min, correlation_min=args.correlation_min,
            local_window_px=args.local_window_px, direction_window_px=args.direction_window_px,
            z_min=args.z_min, anomaly_min_m=args.anomaly_min_m,
            coherence_min=args.coherence_min, minimum_pairs=args.minimum_pairs,
            minimum_component_pixels=args.minimum_component_pixels,
        )
    except (ScreenError, OSError, ValueError, KeyError) as exc:
        print(f"SCREEN_BLOCKED:{type(exc).__name__}:{exc}", file=sys.stderr)
        return 2
    print(f"GOFF_SCREEN_STATUS={result['status']}")
    print(f"GOFF_SCREEN_PAIRS={result['pair_count']}")
    print(f"GOFF_SCREEN_REVIEW_CANDIDATES={result['review_candidate_components']}")
    print("No hazard alert or validated detection is issued.")
    print(f"GOFF_SCREEN_SUMMARY={args.out / 'screen-summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
