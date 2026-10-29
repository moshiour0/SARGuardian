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
from pyproj.transformer import TransformerGroup

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
        # Acquisitions across adjacent subswaths can differ by seconds while
        # representing the same repeat pair; do not count them twice.
        return self.reference[:8], self.secondary[:8]

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
    elevation_ellipsoid_m: np.ndarray | None = None
    range_azimuth_deg: np.ndarray | None = None


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
    dem_vertical_crs: str,
    polarization: str = "HH",
    frequency: str = "frequencyA",
    layer: str = "layer2",
    snr_min: float = 5.0,
    correlation_min: float | None = 0.30,
) -> RawPair:
    identity = parse_identity(path)
    if dem_vertical_crs not in {"EPSG:3855", "ELLIPSOIDAL_WGS84"}:
        raise ScreenError("DEM_VERTICAL_DATUM_REQUIRED: choose EPSG:3855 or ELLIPSOIDAL_WGS84")
    if not path.is_file():
        raise ScreenError("INPUT_HDF5_NOT_FOUND:" + str(path))
    if not dem_path.is_file():
        raise ScreenError("DEM_NOT_FOUND:" + str(dem_path))
    expected_group = f"science/LSAR/GOFF/grids/{frequency}/pixelOffsets/{polarization}/{layer}"
    with h5py.File(path, "r") as h5:
        if expected_group not in h5:
            raise ScreenError("GOFF_LAYER_NOT_FOUND:" + expected_group)

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
    # GOFF radarGrid heightAboveEllipsoid expects ellipsoidal elevation.
    # Use the standard vertical CRS operation and record the actual selected
    # operation. If its required geoid grid is unavailable, fail closed.
    vertical_conversion = "NONE_ALREADY_ELLIPSOIDAL_WGS84"
    if dem_vertical_crs == "EPSG:3855":
        group = TransformerGroup("EPSG:4326+3855", "EPSG:4979", always_xy=True)
        transformer = next(
            (candidate for candidate in group.transformers
             if "EGM2008" in candidate.description.upper()
             or "EGM08" in candidate.description.upper()),
            None,
        )
        if transformer is None or not group.best_available:
            raise ScreenError(
                "EGM2008_VERTICAL_GRID_UNAVAILABLE: install a current PROJ grid "
                "covering the AOI; orthometric DEM cannot be used as ellipsoidal height"
            )
        to_lonlat = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
        xcentres = transform.c + (np.arange(elevation.shape[1], dtype=np.float64) + 0.5) * transform.a
        rows_per_block = 96
        for r0 in range(0, elevation.shape[0], rows_per_block):
            r1 = min(elevation.shape[0], r0 + rows_per_block)
            block = elevation[r0:r1]
            good_block = np.isfinite(block)
            if not good_block.any():
                continue
            ycentres = transform.f + (np.arange(r0, r1, dtype=np.float64) + 0.5) * transform.e
            xx, yy = np.meshgrid(xcentres, ycentres)
            lon, lat = to_lonlat.transform(xx[good_block], yy[good_block])
            _, _, h_ellipsoid = transformer.transform(
                lon, lat, block[good_block].astype(np.float64)
            )
            converted = np.asarray(h_ellipsoid, dtype=np.float64)
            if not np.isfinite(converted).all():
                raise ScreenError("DEM_VERTICAL_CONVERSION_RETURNED_NONFINITE_VALUES")
            elevation[r0:r1][good_block] = converted.astype(np.float32)
        vertical_conversion = "EGM2008_ORTHOMETRIC_TO_WGS84_ELLIPSOIDAL:" + transformer.description
    dem_valid = np.isfinite(elevation)
    if int((dem_valid & aoi_mask).sum()) < max(10, int(0.5 * aoi_mask.sum())):
        raise ScreenError("DEM_VALID_AOI_COVERAGE_BELOW_50_PERCENT")

    h, w = arrays["range"].shape
    east = np.full((h, w), np.nan, dtype=np.float32)
    north = np.full((h, w), np.nan, dtype=np.float32)
    range_azimuth = np.full((h, w), np.nan, dtype=np.float32)
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
        range_azimuth[r0:r1] = np.where(
            block_good, np.degrees(np.arctan2(kx, ky)) % 360.0, np.nan
        ).astype(np.float32)

    return RawPair(
        identity=identity, east_offsets_m=east, north_offsets_m=north,
        snr=arrays["snr"], correlation=arrays["correlation"],
        xs=x, ys=y, crs=crs, transform=transform, grid_mask=aoi_mask,
        source_sha256=sha256_file(path), source_bytes=path.stat().st_size,
        elevation_ellipsoid_m=elevation.astype(np.float32, copy=False),
        range_azimuth_deg=range_azimuth,
        geometry_note=(
            "Per-pixel projection interpolated from this product's radarGrid at DEM elevations converted "
            "to the WGS84 ellipsoid; vertical_conversion=" + vertical_conversion
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
    correlation_min: float | None = DEFAULTS["correlation_min"],
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
        & np.isfinite(magnitude) & (magnitude <= magnitude_clip_m)
    )
    # Exact author-like screening does not gate on correlationSurfacePeak.
    # The stricter SARGuardian mode can opt in to that additional filter.
    if correlation_min is not None:
        valid &= np.isfinite(correlation) & (correlation >= correlation_min)
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
        # Retain only the fields consumed by the streaming accumulator. The
        # intermediate magnitude/z/residual arrays are intentionally released
        # when this function returns; returning them would multiply peak RAM.
        "candidate": hit,
        "anomaly_m": anomaly,
        "directional_coherence": coherence,
        "residual_east_m": residual_east.astype(np.float32, copy=False),
        "residual_north_m": residual_north.astype(np.float32, copy=False),
        "valid": valid & np.asarray(aoi_mask, dtype=bool),
        "valid_aoi_pixels": int((valid & aoi_mask).sum()),
        "candidate_pixels": int(hit.sum()),
        "aoi_grid_pixels": int(np.asarray(aoi_mask, dtype=bool).sum()),
        "candidate_aoi_fraction": float(hit.sum() / max(1, np.asarray(aoi_mask, dtype=bool).sum())),
    }


def validate_identity_series(identities: list[GranuleIdentity]) -> list[GranuleIdentity]:
    """Validate series metadata before reading large source rasters into memory."""
    if not identities:
        raise ScreenError("NO_INPUT_PAIRS")
    key = identities[0].series_key
    for identity in identities:
        if identity.series_key != key:
            raise ScreenError(
                "MIXED_SERIES_IDENTITY: refusing to merge tracks, directions, frames, "
                "bandwidths, polarizations, processing types or CRIDs; got "
                + identity.identity_string
            )
        if identity.secondary <= identity.reference:
            raise ScreenError("INVALID_ACQUISITION_INTERVAL:" + repr(identity.pair_key))
    if len({identity.pair_key for identity in identities}) != len(identities):
        raise ScreenError("DUPLICATE_ACQUISITION_INTERVALS")
    return sorted(identities, key=lambda item: (item.secondary, item.reference))


def _grid_signature(pair: RawPair) -> dict[str, Any]:
    """Keep only small grid metadata plus AOI mask; never retain pair displacement arrays."""
    return {
        "crs": pair.crs,
        "shape": pair.east_offsets_m.shape,
        "xs": pair.xs,
        "ys": pair.ys,
        "transform": pair.transform,
        "grid_mask": pair.grid_mask,
    }


def _grid_compatible_signature(reference: dict[str, Any], other: RawPair) -> bool:
    return (
        reference["crs"] == other.crs
        and reference["shape"] == other.east_offsets_m.shape
        and np.allclose(reference["xs"], other.xs, rtol=0, atol=1e-3)
        and np.allclose(reference["ys"], other.ys, rtol=0, atol=1e-3)
        and np.allclose(tuple(reference["transform"]), tuple(other.transform), rtol=0, atol=1e-6)
        and np.array_equal(reference["grid_mask"], other.grid_mask)
    )


def _grid_compatible(first: RawPair, other: RawPair) -> bool:
    return _grid_compatible_signature(_grid_signature(first), other)


def validate_series(pairs: list[RawPair]) -> None:
    """Compatibility helper for tests; production screening uses streaming validation."""
    if not pairs:
        raise ScreenError("NO_INPUT_PAIRS")
    sorted_ids = validate_identity_series([pair.identity for pair in pairs])
    ordered = {identity.path: pair for identity, pair in zip([p.identity for p in pairs], pairs)}
    pairs[:] = [ordered[identity.path] for identity in sorted_ids]
    reference = pairs[0]
    signature = _grid_signature(reference)
    for pair in pairs[1:]:
        if not _grid_compatible_signature(signature, pair):
            raise ScreenError("INCOMPATIBLE_PIXEL_GRID: no reprojection/resampling is done implicitly")


def terrain_slope_and_downslope(
    elevation_m: np.ndarray, pixel_size_m: float = 80.0
) -> tuple[np.ndarray, np.ndarray]:
    """Horn-like 3x3 slope and downslope azimuth on a north-up projected DEM grid."""
    z = np.asarray(elevation_m, dtype=np.float32)
    if z.ndim != 2 or not np.isfinite(pixel_size_m) or pixel_size_m <= 0:
        raise ScreenError("TERRAIN_GRID_OR_RESOLUTION_INVALID")
    # A central gradient is stable on the regular grid and preserves NoData
    # near missing terrain instead of filling it with an arbitrary elevation.
    dz_north = np.gradient(z, -pixel_size_m, axis=0)
    dz_east = np.gradient(z, pixel_size_m, axis=1)
    slope = np.degrees(np.arctan(np.hypot(dz_east, dz_north))).astype(np.float32)
    # Downslope vector points along negative gradient. Azimuth is clockwise
    # from north; x=east, y=north.
    downslope = (np.degrees(np.arctan2(-dz_east, -dz_north)) % 360.0).astype(np.float32)
    valid = np.isfinite(z) & np.isfinite(dz_east) & np.isfinite(dz_north)
    slope[~valid] = np.nan
    downslope[(~valid) | (slope < 1.0)] = np.nan
    return slope, downslope


def circular_mean_degrees(values: np.ndarray) -> float:
    a = np.radians(np.asarray(values, dtype=np.float64))
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan")
    return float(np.degrees(np.arctan2(np.sin(a).mean(), np.cos(a).mean())) % 360.0)


def angular_difference_degrees(first: float, second: float) -> float:
    if not (math.isfinite(first) and math.isfinite(second)):
        return float("nan")
    return float(abs((first - second + 180.0) % 360.0 - 180.0))


def _component_rows(mask: np.ndarray, transform: Affine, crs: CRS, minimum_pixels: int) -> list[dict[str, Any]]:
    pixel_area = abs(transform.a * transform.e)
    rows: list[dict[str, Any]] = []
    to_wgs84 = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    next_label = 1
    # Use rasterio polygonization with 8-connectivity and count pixels from the
    # polygon area. Each output feature is a review polygon, not a landslide.
    for geometry, value in shapes(mask.astype(np.uint8), mask=mask, transform=transform, connectivity=4):
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


def measure_candidate_pair_series(
    identities: list[GranuleIdentity],
    dem_path: Path,
    ring: list[tuple[float, float]],
    dem_vertical_crs: str,
    components: list[dict[str, Any]],
    grid: dict[str, Any],
    *,
    polarization: str,
    frequency: str,
    layer: str,
    snr_min: float,
    correlation_min: float | None,
    local_window_px: int,
    direction_window_px: int,
    z_min: float,
    anomaly_min_m: float,
    coherence_min: float,
) -> list[dict[str, Any]]:
    """Post-hoc pairwise measurements inside discovered components, never cumulative LOS/3-D velocity."""
    if not components:
        return []
    from rasterio.features import rasterize

    zone_masks = [
        rasterize(
            [(component["geometry"], 1)],
            out_shape=grid["shape"],
            transform=grid["transform"],
            fill=0,
            dtype=np.uint8,
        ).astype(bool)
        for component in components
    ]
    all_zone_mask = np.zeros(grid["shape"], dtype=bool)
    for mask in zone_masks:
        all_zone_mask |= mask
    rows: list[dict[str, Any]] = []

    for identity in identities:
        pair = read_raw_pair(
            identity.path, dem_path, ring, dem_vertical_crs=dem_vertical_crs,
            polarization=polarization, frequency=frequency, layer=layer,
            snr_min=snr_min, correlation_min=correlation_min,
        )
        if pair.identity.series_key != identities[0].series_key:
            raise ScreenError("PAIR_SERIES_IDENTITY_CHANGED_DURING_MEASUREMENT")
        if not _grid_compatible_signature(grid, pair):
            raise ScreenError("INCOMPATIBLE_PIXEL_GRID_DURING_MEASUREMENT:" + identity.path.name)

        result = score_pair_arrays(
            pair.east_offsets_m, pair.north_offsets_m, pair.snr, pair.correlation,
            pair.grid_mask, snr_min=snr_min, correlation_min=correlation_min,
            local_window_px=local_window_px, direction_window_px=direction_window_px,
            z_min=z_min, anomaly_min_m=anomaly_min_m, coherence_min=coherence_min,
        )
        residual_e = result["residual_east_m"]
        residual_n = result["residual_north_m"]
        residual_mag = np.hypot(residual_e, residual_n)
        control = (
            result["valid"] & ~result["candidate"] & ~all_zone_mask
            & np.isfinite(residual_mag)
        )
        control_values = residual_mag[control]
        control_median = float(np.median(control_values)) if control_values.size else None
        control_p95 = float(np.percentile(control_values, 95)) if control_values.size else None
        interval_days = (
            datetime.strptime(identity.secondary, "%Y%m%dT%H%M%S")
            - datetime.strptime(identity.reference, "%Y%m%dT%H%M%S")
        ).total_seconds() / 86400.0

        for component, zone_mask in zip(components, zone_masks):
            zone_pixels = int(zone_mask.sum())
            zone_valid = zone_mask & result["valid"]
            zone_hits = zone_mask & result["candidate"]
            n_valid = int(zone_valid.sum())
            n_hits = int(zone_hits.sum())
            zone_residual_values = residual_mag[zone_valid]
            zone_anomaly_values = result["anomaly_m"][zone_valid]
            hit_e = residual_e[zone_hits]
            hit_n = residual_n[zone_hits]
            hit_coherence = result["directional_coherence"][zone_hits]
            hit_anomaly = result["anomaly_m"][zone_hits]

            mean_e = float(np.mean(hit_e)) if hit_e.size else None
            mean_n = float(np.mean(hit_n)) if hit_n.size else None
            mean_vec = float(math.hypot(mean_e, mean_n)) if mean_e is not None and mean_n is not None else None
            azimuth = (
                float(np.degrees(np.arctan2(mean_e, mean_n)) % 360.0)
                if mean_e is not None and mean_n is not None else None
            )
            rows.append({
                "candidate_id": component["candidate_id"],
                "reference_date": identity.reference[:8],
                "secondary_date": identity.secondary[:8],
                "interval_days": round(interval_days, 6),
                "series_identity": identity.identity_string,
                "crid": identity.crid,
                "zone_pixels": zone_pixels,
                "valid_zone_pixels": n_valid,
                "valid_zone_fraction": round(n_valid / max(1, zone_pixels), 6),
                "screen_hit_pixels_in_zone": n_hits,
                "screen_hit_fraction_of_zone": round(n_hits / max(1, zone_pixels), 6),
                "zone_median_local_residual_magnitude_m": (
                    round(float(np.median(zone_residual_values)), 6) if zone_residual_values.size else None
                ),
                "zone_median_local_anomaly_m": (
                    round(float(np.median(zone_anomaly_values)), 6) if zone_anomaly_values.size else None
                ),
                "median_anomaly_on_screen_hits_m": (
                    round(float(np.median(hit_anomaly)), 6) if hit_anomaly.size else None
                ),
                "screen_hit_mean_residual_east_m": round(mean_e, 6) if mean_e is not None else None,
                "screen_hit_mean_residual_north_m": round(mean_n, 6) if mean_n is not None else None,
                "screen_hit_mean_residual_vector_m": round(mean_vec, 6) if mean_vec is not None else None,
                "screen_hit_relative_rate_m_per_day": (
                    round(mean_vec / interval_days, 8) if mean_vec is not None and interval_days > 0 else None
                ),
                "screen_hit_mean_vector_azimuth_deg": round(azimuth, 3) if azimuth is not None else None,
                "screen_hit_median_directional_coherence": (
                    round(float(np.median(hit_coherence)), 6) if hit_coherence.size else None
                ),
                "control_pixels": int(control_values.size),
                "control_median_local_residual_magnitude_m": round(control_median, 6) if control_median is not None else None,
                "control_p95_local_residual_magnitude_m": round(control_p95, 6) if control_p95 is not None else None,
                "zone_minus_control_median_residual_m": (
                    round(float(np.median(zone_residual_values)) - control_median, 6)
                    if zone_residual_values.size and control_median is not None else None
                ),
                "interpretation": "PAIRWISE_RELATIVE_DIAGNOSTIC_NOT_CUMULATIVE_VELOCITY",
            })
        del pair, result, residual_e, residual_n, residual_mag

    return rows


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
    correlation_min: float | None = DEFAULTS["correlation_min"],
    local_window_px: int = DEFAULTS["local_window_px"],
    direction_window_px: int = DEFAULTS["direction_window_px"],
    z_min: float = DEFAULTS["z_min"],
    anomaly_min_m: float = DEFAULTS["anomaly_min_m"],
    coherence_min: float = DEFAULTS["directional_coherence_min"],
    minimum_pairs: int = DEFAULTS["minimum_pairs"],
    minimum_component_pixels: int = DEFAULTS["minimum_component_pixels"],
    dem_vertical_crs: str | None = None,
) -> dict[str, Any]:
    if dem_vertical_crs not in {"EPSG:3855", "ELLIPSOIDAL_WGS84"}:
        raise ScreenError("DEM_VERTICAL_DATUM_REQUIRED: choose EPSG:3855 or ELLIPSOIDAL_WGS84")
    if minimum_pairs < 2:
        raise ScreenError("MINIMUM_PAIRS_MUST_BE_AT_LEAST_TWO")
    if minimum_component_pixels < 1:
        raise ScreenError("MINIMUM_COMPONENT_PIXELS_MUST_BE_POSITIVE")
    if len(inputs) != len(set(p.resolve() for p in inputs)):
        raise ScreenError("DUPLICATE_INPUT_PATHS")
    ring = parse_polygon_wkt(aoi_wkt)
    identities = validate_identity_series([parse_identity(path) for path in inputs])
    out_dir.mkdir(parents=True, exist_ok=True)

    # Stream one raw HDF5 pair at a time. GOFF granules and derived 80 m arrays
    # are large; retaining every pair's E/N/SNR/correlation arrays at once is
    # needlessly memory-hungry and can make a real multi-pair run fail.
    grid = None
    hit_count = anomaly_sum = coherence_sum = None
    hit_east_sum = hit_north_sum = None
    terrain_elevation = terrain_range_azimuth = None
    reports: list[dict[str, Any]] = []

    for identity in identities:
        pair = read_raw_pair(
            identity.path, dem_path, ring, dem_vertical_crs=dem_vertical_crs,
            polarization=polarization, frequency=frequency, layer=layer, snr_min=snr_min,
            correlation_min=correlation_min,
        )
        if pair.identity.series_key != identities[0].series_key or pair.identity.pair_key != identity.pair_key:
            raise ScreenError("INPUT_IDENTITY_CHANGED_DURING_READ:" + identity.path.name)

        if grid is None:
            grid = _grid_signature(pair)
            hit_count = np.zeros(pair.east_offsets_m.shape, dtype=np.uint16)
            anomaly_sum = np.zeros(pair.east_offsets_m.shape, dtype=np.float32)
            coherence_sum = np.zeros(pair.east_offsets_m.shape, dtype=np.float32)
            hit_east_sum = np.zeros(pair.east_offsets_m.shape, dtype=np.float32)
            hit_north_sum = np.zeros(pair.east_offsets_m.shape, dtype=np.float32)
            terrain_elevation = pair.elevation_ellipsoid_m
            terrain_range_azimuth = pair.range_azimuth_deg
        elif not _grid_compatible_signature(grid, pair):
            raise ScreenError(
                "INCOMPATIBLE_PIXEL_GRID: no reprojection/resampling is done implicitly: "
                + identity.path.name
            )

        result = score_pair_arrays(
            pair.east_offsets_m, pair.north_offsets_m, pair.snr, pair.correlation,
            pair.grid_mask, snr_min=snr_min, correlation_min=correlation_min,
            local_window_px=local_window_px, direction_window_px=direction_window_px,
            z_min=z_min, anomaly_min_m=anomaly_min_m, coherence_min=coherence_min,
        )
        hit = result["candidate"]
        hit_count += hit.astype(np.uint16)
        anomaly_sum += np.where(hit, result["anomaly_m"], 0.0).astype(np.float32)
        coherence_sum += np.where(hit, result["directional_coherence"], 0.0).astype(np.float32)
        hit_east_sum += np.where(hit, result["residual_east_m"], 0.0).astype(np.float32)
        hit_north_sum += np.where(hit, result["residual_north_m"], 0.0).astype(np.float32)
        reports.append({
            "granule": identity.path.name,
            "granule_identity": identity.identity_string,
            "subswath": identity.subswath,
            "reference": identity.reference,
            "secondary": identity.secondary,
            "interval_days": (
                datetime.strptime(identity.secondary, "%Y%m%dT%H%M%S")
                - datetime.strptime(identity.reference, "%Y%m%dT%H%M%S")
            ).total_seconds() / 86400.0,
            "source_sha256": pair.source_sha256,
            "source_bytes": pair.source_bytes,
            "aoi_grid_pixels": result["aoi_grid_pixels"],
            "valid_aoi_pixels": result["valid_aoi_pixels"],
            "candidate_pixels": result["candidate_pixels"],
            "candidate_aoi_fraction": result["candidate_aoi_fraction"],
            "status": "SCREENED_NOT_DETECTION",
        })
        # Release the per-pair rasters before opening the next raw granule.
        del result, hit, pair

    if (grid is None or hit_count is None or anomaly_sum is None or coherence_sum is None
            or hit_east_sum is None or hit_north_sum is None
            or terrain_elevation is None or terrain_range_azimuth is None):
        raise ScreenError("NO_PAIRS_PROCESSED")
    pair_count = len(reports)
    persistent = hit_count >= minimum_pairs
    anomaly_mean = np.full(hit_count.shape, np.nan, dtype=np.float32)
    coherence_mean = np.full(hit_count.shape, np.nan, dtype=np.float32)
    flagged = hit_count > 0
    anomaly_mean[flagged] = anomaly_sum[flagged] / hit_count[flagged]
    coherence_mean[flagged] = coherence_sum[flagged] / hit_count[flagged]
    components = _component_rows(persistent, grid["transform"], grid["crs"], minimum_component_pixels)
    pair_series_rows = measure_candidate_pair_series(
        identities, dem_path, ring, dem_vertical_crs, components, grid,
        polarization=polarization, frequency=frequency, layer=layer,
        snr_min=snr_min, correlation_min=correlation_min,
        local_window_px=local_window_px, direction_window_px=direction_window_px,
        z_min=z_min, anomaly_min_m=anomaly_min_m, coherence_min=coherence_min,
    )
    terrain_slope, terrain_downslope = terrain_slope_and_downslope(
        terrain_elevation, pixel_size_m=abs(float(grid["transform"].a))
    )
    # Independent post-hoc directional attribution. It does NOT gate the
    # candidate mask: terrain-following motion is evidence for review, not proof.
    candidate_hit = (hit_count > 0) & persistent
    pixel_disp_az = (np.degrees(np.arctan2(hit_east_sum, hit_north_sum)) % 360.0).astype(np.float32)
    valid_direction = candidate_hit & np.isfinite(pixel_disp_az) & np.isfinite(terrain_downslope)
    pixel_misfit = np.full(hit_count.shape, np.nan, dtype=np.float32)
    pixel_misfit[valid_direction] = np.abs(
        (pixel_disp_az[valid_direction] - terrain_downslope[valid_direction] + 180.0) % 360.0 - 180.0
    )
    for component in components:
        from rasterio.features import rasterize
        component_mask = rasterize(
            [(component["geometry"], 1)], out_shape=persistent.shape,
            transform=grid["transform"], fill=0, dtype=np.uint8,
        ).astype(bool)
        pixels = component_mask & candidate_hit
        component["mean_residual_azimuth_deg"] = round(
            float(np.degrees(np.arctan2(hit_east_sum[pixels].sum(), hit_north_sum[pixels].sum())) % 360.0), 2
        ) if pixels.any() else None
        component["mean_downslope_azimuth_deg"] = round(circular_mean_degrees(terrain_downslope[component_mask]), 2)
        component["mean_slope_deg"] = round(float(np.nanmean(terrain_slope[component_mask])), 2) if np.isfinite(terrain_slope[component_mask]).any() else None
        component["mean_elevation_ellipsoidal_m"] = round(float(np.nanmean(terrain_elevation[component_mask])), 1) if np.isfinite(terrain_elevation[component_mask]).any() else None
        component["downsloping_direction_misfit_deg"] = round(
            angular_difference_degrees(
                component["mean_residual_azimuth_deg"]
                if component["mean_residual_azimuth_deg"] is not None else float("nan"),
                component["mean_downslope_azimuth_deg"],
            ), 2
        )
        component["pixel_fraction_within_45deg_downslope"] = round(
            float(np.mean(pixel_misfit[pixels & valid_direction] <= 45.0)), 3
        ) if (pixels & valid_direction).any() else None
        range_mean = circular_mean_degrees(terrain_range_azimuth[component_mask])
        component["mean_ground_range_azimuth_deg"] = round(range_mean, 2) if math.isfinite(range_mean) else None
        component["range_direction_misfit_deg"] = round(
            angular_difference_degrees(
                component["mean_residual_azimuth_deg"]
                if component["mean_residual_azimuth_deg"] is not None else float("nan"),
                component["mean_ground_range_azimuth_deg"],
            ), 2
        )
    qualified_mask = np.zeros(persistent.shape, dtype=bool)
    if components:
        polygon_mask = np.zeros(persistent.shape, dtype=np.uint8)
        from rasterio.features import rasterize
        burn = [(row["geometry"], 1) for row in components]
        rasterize(burn, out=polygon_mask, transform=grid["transform"], fill=0, dtype=np.uint8)
        qualified_mask = polygon_mask.astype(bool)

    status = "EXPLORATORY_RANKING_ONLY" if pair_count >= minimum_pairs else "NO_TEMPORAL_PERSISTENCE_TEST"
    reasons = [
        "Candidate thresholds are adapted from one published NISAR GOFF case study and are not calibrated for arbitrary AOIs.",
        (
            "Correlation quality gate is disabled to match the author's SNR-only screen; this remains an implementation comparison, not proof of physical movement."
            if correlation_min is None else
            f"Correlation >= {correlation_min:g} is an added SARGuardian gate, so this run differs from the author's SNR-only screen."
        ),
        "Candidate-zone pair metrics are local, pairwise relative diagnostics; they are not cumulative displacement or validated velocity.",
        "Terrain downslope/range alignment is reported post hoc, but no matched terrain/surface-type spatial-null significance or glacier attribution is calculated here.",
        "The declared DEM vertical datum is converted to ellipsoidal height using the recorded EGM2008 grid; DEM/geoid uncertainty is not quantified in this prototype.",
        "No spatial-null/multiple-testing correction or held-out positive/negative-control benchmark is applied; candidate ranks are exploratory only.",
    ]
    if pair_count < minimum_pairs:
        reasons.append(
            f"Only {pair_count} compatible pair(s) were supplied; at least {minimum_pairs} are required "
            "for even an exploratory persistence screen."
        )

    hit_count_path = out_dir / "pair_persistence_count.tif"
    anomaly_path = out_dir / "mean_anomaly_m.tif"
    coherence_path = out_dir / "mean_directional_coherence.tif"
    candidates_path = out_dir / "persistent_candidate_mask.tif"
    _write_raster(hit_count_path, hit_count, grid["transform"], grid["crs"], 65535,
                  "Number of GOFF pairs with a coherent local anomaly")
    _write_raster(anomaly_path, anomaly_mean, grid["transform"], grid["crs"], np.nan,
                  "Mean positive local displacement-magnitude anomaly among flagged pairs (m)")
    _write_raster(coherence_path, coherence_mean, grid["transform"], grid["crs"], np.nan,
                  "Mean residual directional coherence among flagged pairs")
    _write_raster(candidates_path, qualified_mask.astype(np.uint8), grid["transform"], grid["crs"], 255,
                  "Exploratory persistent-component mask; 1=review candidate, not hazard")

    candidate_features = []
    for row in components:
        geometry_wgs84 = transform_geom(grid["crs"], "EPSG:4326", row["geometry"])
        candidate_features.append({
            "type": "Feature",
            "properties": {k: v for k, v in row.items() if k != "geometry"},
            "geometry": geometry_wgs84,
        })
    geojson_path = out_dir / "candidates.geojson"
    geojson_path.write_text(json.dumps({
        "type": "FeatureCollection", "features": candidate_features,
        "properties": {"status": status, "interpretation": "Exploratory review candidates; not hazard detections."},
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    csv_path = out_dir / "candidates.csv"
    pair_series_path = out_dir / "candidate_pair_series.csv"
    pair_series_fields = [
        "candidate_id", "reference_date", "secondary_date", "interval_days", "series_identity", "crid",
        "zone_pixels", "valid_zone_pixels", "valid_zone_fraction", "screen_hit_pixels_in_zone",
        "screen_hit_fraction_of_zone", "zone_median_local_residual_magnitude_m",
        "zone_median_local_anomaly_m", "median_anomaly_on_screen_hits_m",
        "screen_hit_mean_residual_east_m", "screen_hit_mean_residual_north_m",
        "screen_hit_mean_residual_vector_m", "screen_hit_relative_rate_m_per_day",
        "screen_hit_mean_vector_azimuth_deg", "screen_hit_median_directional_coherence",
        "control_pixels", "control_median_local_residual_magnitude_m",
        "control_p95_local_residual_magnitude_m", "zone_minus_control_median_residual_m", "interpretation",
    ]
    with pair_series_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=pair_series_fields)
        writer.writeheader()
        writer.writerows(pair_series_rows)
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "candidate_id", "pixels_approx", "area_km2_approx",
                "bbox_center_lon", "bbox_center_lat",
                "mean_residual_azimuth_deg", "mean_downslope_azimuth_deg",
                "mean_slope_deg", "mean_elevation_ellipsoidal_m",
                "downsloping_direction_misfit_deg",
                "pixel_fraction_within_45deg_downslope",
                "mean_ground_range_azimuth_deg", "range_direction_misfit_deg",
            ],
        )
        writer.writeheader()
        for row in components:
            writer.writerow({k: row[k] for k in writer.fieldnames})

    summary = {
        "schema_version": SCHEMA,
        "status": status,
        "interpretation": "Exploratory target-free ranking only. No hazard alert, failure prediction or validated detection is issued.",
        "series_identity": identities[0].identity_string,
        "processing": identities[0].processing,
        "crid": identities[0].crid,
        "pair_count": pair_count,
        "minimum_pairs_for_persistence": minimum_pairs,
        "aoi_wkt": aoi_wkt,
        "aoi_grid_pixels": int(grid["grid_mask"].sum()),
        "grid_shape": list(grid["shape"]),
        "grid_crs": grid["crs"].to_string(),
        "grid_transform": list(grid["transform"]),
        "dem": {
            "path": str(dem_path), "sha256": sha256_file(dem_path),
            "declared_vertical_crs": dem_vertical_crs,
            "vertical_conversion": (
                "EGM2008_ORTHOMETRIC_TO_WGS84_ELLIPSOIDAL_US_NGA_EGM08_25"
                if dem_vertical_crs == "EPSG:3855" else "NONE_ALREADY_ELLIPSOIDAL_WGS84"
            ),
        },
        "parameters": {
            "snr_min": snr_min, "correlation_min": correlation_min,
            "magnitude_clip_m": DEFAULTS["magnitude_clip_m"],
            "local_window_px": local_window_px, "direction_window_px": direction_window_px,
            "z_min": z_min, "anomaly_min_m": anomaly_min_m,
            "directional_coherence_min": coherence_min,
            "minimum_pairs": minimum_pairs, "minimum_component_pixels": minimum_component_pixels,
        },
        "per_pair_statistics": [
            {
                "granule": report["granule"],
                "valid_aoi_pixels": report["valid_aoi_pixels"],
                "candidate_pixels": report["candidate_pixels"],
                "candidate_aoi_fraction": report["candidate_aoi_fraction"],
                "status": report["status"],
            }
            for report in reports
        ],
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
        "candidate_pair_series": {
            "path": pair_series_path.name,
            "rows": len(pair_series_rows),
            "status": "PAIRWISE_RELATIVE_DIAGNOSTICS_NOT_CUMULATIVE_VELOCITY",
            "control_reference": "AOI-valid pixels outside all candidate zones and non-hit in that pair; descriptive comparator, NOT a calibrated detection floor.",
        },
        "terrain_directional_review": {
            "status": "POST_HOC_ONLY_NO_ALERT",
            "definition": "Angle between mean residual displacement azimuth and DEM-derived downslope azimuth; ground-range alignment is reported as an artefact comparator.",
            "pixels_with_terrain_and_motion_direction": int(valid_direction.sum()),
            "fraction_of_review_pixels_within_45deg_downslope": (
                round(float(np.mean(pixel_misfit[valid_direction] <= 45.0)), 4)
                if valid_direction.any() else None
            ),
            "method": "Central DEM gradient on the native projected grid; flat (<1 degree) or invalid terrain has no aspect.",
            "limitations": "No matched spatial-null p-value, Rayleigh/V-test, glacier inventory join or false-alarm calibration is implemented.",
        },
        "reasons_and_limitations": reasons,
        "outputs": {
            "pair_persistence_count": hit_count_path.name,
            "mean_anomaly": anomaly_path.name,
            "mean_directional_coherence": coherence_path.name,
            "persistent_candidate_mask": candidates_path.name,
            "candidates_geojson": geojson_path.name,
            "candidates_csv": csv_path.name,
            "candidate_pair_series_csv": pair_series_path.name,
        },
    }
    summary["output_sha256"] = {
        name: sha256_file(out_dir / filename) for name, filename in summary["outputs"].items()
    }
    (out_dir / "screen-summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8"
    )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--read", nargs="+", required=True, type=Path,
                        help="Raw NISAR L2 GOFF HDF5 products from one compatible series")
    parser.add_argument("--dem", required=True, type=Path,
                        help="Georeferenced DEM raster")
    parser.add_argument("--dem-vertical-crs", required=True,
                        choices=["EPSG:3855", "ELLIPSOIDAL_WGS84"],
                        help="Explicit DEM height datum; EPSG:3855 is converted with the installed EGM2008 PROJ grid")
    parser.add_argument("--aoi-wkt", required=True,
                        help="One WGS84 POLYGON WKT. Never defaults to the AOI centre.")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--polarization", default="HH")
    parser.add_argument("--frequency", default="frequencyA")
    parser.add_argument("--layer", default="layer2")
    parser.add_argument("--snr-min", type=float, default=DEFAULTS["snr_min"])
    def parse_optional_correlation(value: str) -> float | None:
        if value.strip().lower() in {"none", "off", "disabled"}:
            return None
        threshold = float(value)
        if not math.isfinite(threshold):
            raise argparse.ArgumentTypeError("correlation threshold must be finite or 'none'")
        return threshold

    parser.add_argument(
        "--correlation-min", type=parse_optional_correlation,
        default=DEFAULTS["correlation_min"],
        help="Correlation quality threshold; pass 'none' to match the author screen's SNR-only gate",
    )
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
            dem_vertical_crs=args.dem_vertical_crs,
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
