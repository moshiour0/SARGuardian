#!/usr/bin/env python3
"""Refresh the checked-in browser bundle after the exact source AOI is enforced.

This is intentionally an all-or-nothing transform: if the current grid cannot
be polygon-masked, candidate IDs no longer correspond to the saved bound table,
scene cloud statistics cannot be recalculated, or essential validation fails,
the script exits before writing either site_data.json or index.html.
"""
from __future__ import annotations

import argparse
import base64
import io
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(ROOT / "docs" / "site") not in sys.path:
    sys.path.insert(0, str(ROOT / "docs" / "site"))

import build_data  # noqa: E402
from gunw_reader import AOIS  # noqa: E402
from scar_map import (GRID, enumerate_candidates, classify, grid_transform,
                      source_aoi_mask)  # noqa: E402
from local_floor import point_in_aoi  # noqa: E402


SCHEMA = 2
AOI_VERSION = "source-ring-v1"
EVENT_GRID_BBOX = [85.4645, 28.2453, 85.5562, 28.3529]


def polygon_and_affine() -> tuple[list[list[float]], list[float]]:
    from pyproj import Transformer
    from gunw_reader import SOURCE_RING

    aff = grid_transform()
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:32645", always_xy=True)
    projected = [list(transformer.transform(lon, lat)) for lon, lat in SOURCE_RING]
    affine = [float(aff.c), float(aff.a), float(aff.b),
              float(aff.f), float(aff.d), float(aff.e)]
    return projected, affine


def refresh(data: dict, template: str, *, scenes_cache: Path) -> tuple[dict, str]:
    if not isinstance(data.get("meta"), dict) or not isinstance(data.get("clusters"), list):
        raise ValueError("site_data.json schema is incomplete.")
    grid = data["meta"].get("grid", {})
    width, height = int(grid.get("w", 0)), int(grid.get("h", 0))
    if (height, width) != GRID:
        raise ValueError(f"Grid mismatch: JSON says {width}x{height}, expected {GRID[1]}x{GRID[0]}.")
    image_bytes = base64.b64decode(data["change_png"], validate=True)
    with Image.open(io.BytesIO(image_bytes)) as image:
        image = image.convert("RGB")
        if image.size != (width, height):
            raise ValueError(f"PNG grid mismatch: {image.size}, expected {(width, height)}.")
        rgb = np.asarray(image, dtype=np.uint8).copy()

    exact_aoi = np.asarray(source_aoi_mask(), dtype=bool)
    if exact_aoi.shape != (height, width):
        raise ValueError(f"Polygon mask shape mismatch: {exact_aoi.shape}.")
    original_both = rgb[:, :, 2] > 128
    both = original_both & exact_aoi
    rgb[~exact_aoi, 0] = 127
    rgb[~exact_aoi, 1] = 127
    rgb[:, :, 2] = np.where(both, 255, 0).astype(np.uint8)

    # Cluster thresholds are run on the same quantized channels used by the
    # browser. After applying the exact polygon mask, reclassify those clusters
    # rather than merely hiding an old out-of-AOI centre.
    d = rgb[:, :, 0].astype(np.float32) / 127.5 - 1.0
    pre = rgb[:, :, 1].astype(np.float32) / 127.5 - 1.0
    d[~both] = np.nan
    pre[~exact_aoi] = np.nan
    candidate_cells = enumerate_candidates(both, pre, d, -0.25, 0.35, 40)
    classified = classify(candidate_cells, d, pre, grid_transform(), limit=12)
    for row in classified:
        for key in ("area_km2", "lat", "lon", "d_ndsi_median", "elev_m",
                    "slope_deg", "aspect_deg"):
            row[key] = round(float(row[key]), 5)
        if not point_in_aoi(float(row["lon"]), float(row["lat"]), "source"):
            raise ValueError(f"Candidate {row.get('id')} is outside SOURCE_RING after polygon masking.")
    if len(classified) < 2 or sum(r["reading"] == "DETACHMENT-like" for r in classified) < 2:
        raise ValueError("Reclassified polygon-masked grid no longer supports the site's multiple-candidate conclusion; refusing to publish.")

    # Candidate ids used by the bound table refer to the earlier component
    # ranking. Do not publish a mixed set of IDs if exact clipping reorders one.
    old_candidates_path = ROOT / "outputs" / "scar_candidates_source.csv"
    with old_candidates_path.open(newline="", encoding="utf-8") as stream:
        old_rows = {str(r["id"]): r for r in __import__("csv").DictReader(stream)}
    with (ROOT / "outputs" / "bound_source.csv").open(newline="", encoding="utf-8") as stream:
        bound_rows = list(__import__("csv").DictReader(stream))
    bound_rows = [
        r for r in bound_rows
        if point_in_aoi(float(r["lon"]), float(r["lat"]), "source")
    ]
    for bound_row in bound_rows:
        cid = str(bound_row["id"])
        if cid not in old_rows:
            raise ValueError(f"Bound candidate {cid} has no source candidate row.")
        prior = old_rows[cid]
        matches = [r for r in classified if str(r["id"]) == cid]
        if len(matches) != 1:
            raise ValueError(f"Candidate ID {cid} changed ranking under exact masking; do not publish old bound values against it.")
        old_lat, old_lon = float(prior["lat"]), float(prior["lon"])
        new_lat, new_lon = float(matches[0]["lat"]), float(matches[0]["lon"])
        if max(abs(new_lat-old_lat), abs(new_lon-old_lon)) > 0.003:
            raise ValueError(f"Candidate ID {cid} moved materially after exact masking; regenerate its floor/bound inputs before publishing.")

    # Refresh AOI cloud statistics with the current exact-polygon implementation,
    # ignoring any legacy scene cache that counted the bounding rectangle.
    scenes = build_data.scenes_block(scenes_cache)
    if not scenes:
        raise ValueError("Scene refresh returned no Sentinel-2 scenes.")

    png_io = io.BytesIO()
    Image.fromarray(rgb, "RGB").save(png_io, format="PNG", optimize=True)
    data["change_png"] = base64.b64encode(png_io.getvalue()).decode("ascii")

    projected_ring, affine = polygon_and_affine()
    aoi_px = int(exact_aoi.sum())
    usable_px = int(both.sum())
    if aoi_px <= 0 or usable_px > aoi_px:
        raise ValueError(f"Invalid polygon-grid counts: usable={usable_px}, AOI={aoi_px}.")
    data["meta"]["built"] = __import__("datetime").date.today().isoformat()
    grid.update({
        "w": width, "h": height,
        "aoi_px": aoi_px, "usable_px": usable_px,
        "usable_pct": round(100.0 * usable_px / aoi_px, 1),
        "aoi_geometry": AOI_VERSION,
        "affine_utm": affine,
        "aoi_polygon_projected": projected_ring,
        "aoi_polygon_wgs84": [[float(lon), float(lat)] for lon, lat in AOIS["source"]],
        "aoi": EVENT_GRID_BBOX,
    })
    data["scenes"] = scenes
    data["clusters"] = classified
    data["bound"]["candidates"] = bound_rows

    output = json.dumps(data, separators=(",", ":"), ensure_ascii=False)
    if "</script" in output.lower():
        raise ValueError("Updated site payload contains a closing script tag.")
    if "{{DATA}}" not in template:
        raise ValueError("Template no longer contains the {{DATA}} placeholder.")
    html = template.replace("{{DATA}}", output)
    return data, html


def validate_current(site_path: Path, template_path: Path, html_path: Path) -> dict:
    data = json.loads(site_path.read_text(encoding="utf-8"))
    grid = data.get("meta", {}).get("grid", {})
    if grid.get("aoi_geometry") != AOI_VERSION:
        raise ValueError("site_data.json is not refreshed to the exact source polygon.")
    projected, affine = polygon_and_affine()
    if grid.get("affine_utm") != affine or grid.get("aoi_polygon_projected") != projected:
        raise ValueError("site_data.json projected polygon/affine metadata does not match the source geometry.")
    image_bytes = base64.b64decode(data["change_png"], validate=True)
    with Image.open(io.BytesIO(image_bytes)) as image:
        rgb = np.asarray(image.convert("RGB"), dtype=np.uint8)
    exact_aoi = np.asarray(source_aoi_mask(), dtype=bool)
    if rgb.shape[:2] != exact_aoi.shape:
        raise ValueError("Site PNG and source polygon mask dimensions differ.")
    if np.any(rgb[:, :, 2][~exact_aoi] > 128):
        raise ValueError("Site change PNG still marks pixels outside SOURCE_RING as usable.")
    if int(grid.get("aoi_px", -1)) != int(exact_aoi.sum()):
        raise ValueError("Site AOI pixel denominator does not match SOURCE_RING.")
    usable = int(np.count_nonzero((rgb[:, :, 2] > 128) & exact_aoi))
    if int(grid.get("usable_px", -1)) != usable:
        raise ValueError("Site usable pixel count does not match the polygon-masked PNG.")
    expect_pct = round(100.0 * usable / int(exact_aoi.sum()), 1)
    if float(grid.get("usable_pct", -1)) != expect_pct:
        raise ValueError("Site usable percentage does not use the exact polygon denominator.")
    for row in data.get("clusters", []):
        if not point_in_aoi(float(row["lon"]), float(row["lat"]), "source"):
            raise ValueError(f"Out-of-AOI candidate {row.get('id')} remains in site cluster table.")
    for row in data.get("bound", {}).get("candidates", []):
        if not point_in_aoi(float(row["lon"]), float(row["lat"]), "source"):
            raise ValueError(f"Out-of-AOI candidate {row.get('id')} remains in site bound table.")
    if sum(row.get("reading") == "DETACHMENT-like" for row in data.get("clusters", [])) < 2:
        raise ValueError("Site bundle no longer supports multiple detachment-like candidates.")
    payload = json.dumps(data, separators=(",", ":"), ensure_ascii=False)
    if "</script" in payload.lower():
        raise ValueError("Site payload contains a closing script tag.")
    expected_html = template_path.read_text(encoding="utf-8").replace("{{DATA}}", payload)
    if html_path.read_text(encoding="utf-8") != expected_html:
        raise ValueError("index.html is not generated from current site_data.json and template.")
    return data


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-data", type=Path, default=ROOT / "docs/site/site_data.json")
    parser.add_argument("--template", type=Path, default=ROOT / "docs/site/index.tpl.html")
    parser.add_argument("--html", type=Path, default=ROOT / "docs/site/index.html")
    parser.add_argument("--scenes-cache", type=Path, default=ROOT / "docs/site/_scenes_polygon_v2.json")
    parser.add_argument("--check", action="store_true", help="validate existing polygon-aware bundle without network calls")
    args = parser.parse_args()

    template = args.template.read_text(encoding="utf-8")
    if args.check:
        data = validate_current(args.site_data, args.template, args.html)
        print(
            "SITE_POLYGON_BUNDLE=PASS "
            f"clusters={len(data['clusters'])} bounds={len(data['bound']['candidates'])} "
            f"grid_usable={data['meta']['grid']['usable_px']}/{data['meta']['grid']['aoi_px']}"
        )
        return 0

    old_data = json.loads(args.site_data.read_text(encoding="utf-8"))
    refreshed, html = refresh(old_data, template, scenes_cache=args.scenes_cache)
    # Validate the in-memory output before writing either file.
    temp_payload = json.dumps(refreshed, separators=(",", ":"), ensure_ascii=False)
    if "</script" in temp_payload.lower():
        raise ValueError("Updated site payload contains a closing script tag.")
    if html != args.template.read_text(encoding="utf-8").replace("{{DATA}}", temp_payload):
        raise ValueError("Generated HTML does not embed the exact generated JSON.")

    # Write atomically so an interruption cannot leave mismatched HTML/JSON.
    for path, value in ((args.site_data, temp_payload), (args.html, html)):
        path.parent.mkdir(parents=True, exist_ok=True)
        temp = path.with_suffix(path.suffix + ".tmp")
        temp.write_text(value, encoding="utf-8")
        temp.replace(path)
    validate_current(args.site_data, args.template, args.html)
    print(
        "SITE_POLYGON_BUNDLE=REFRESHED "
        f"clusters={len(refreshed['clusters'])} "
        f"detachment_like={sum(r['reading']=='DETACHMENT-like' for r in refreshed['clusters'])} "
        f"grid_usable={refreshed['meta']['grid']['usable_px']}/{refreshed['meta']['grid']['aoi_px']} "
        f"bounds={len(refreshed['bound']['candidates'])} scenes={len(refreshed['scenes'])}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
