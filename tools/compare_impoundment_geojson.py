#!/usr/bin/env python3
"""Compare baseline and regenerated dam GeoJSONs without unsafe nearest-site pairing.

Coordinates define site identity for this audit. A nearest neighbour is never
silently treated as the same site, because a hydrology model may select a
different channel cell and thereby change the inferred sill and pool.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_geojson(path: Path) -> list[dict[str, Any]]:
    obj = json.loads(path.read_text(encoding="utf-8"))
    if obj.get("type") != "FeatureCollection" or not isinstance(obj.get("features"), list):
        raise ValueError(f"{path}: expected a GeoJSON FeatureCollection.")
    rows = []
    for index, feature in enumerate(obj["features"]):
        if feature.get("type") != "Feature":
            raise ValueError(f"{path}: feature {index} has invalid type.")
        geometry = feature.get("geometry") or {}
        if geometry.get("type") != "Point":
            raise ValueError(f"{path}: feature {index} is not a Point.")
        coords = geometry.get("coordinates")
        if not isinstance(coords, list) or len(coords) < 2:
            raise ValueError(f"{path}: feature {index} has invalid coordinates.")
        lon, lat = float(coords[0]), float(coords[1])
        if not (math.isfinite(lon) and math.isfinite(lat) and -180 <= lon <= 180 and -90 <= lat <= 90):
            raise ValueError(f"{path}: feature {index} has non-finite/out-of-range WGS84 coordinates.")
        props = feature.get("properties") or {}
        if not isinstance(props, dict):
            raise ValueError(f"{path}: feature {index} properties must be an object.")
        key = (lon, lat)
        if any(row["coordinate"] == key for row in rows):
            raise ValueError(f"{path}: duplicate exact site coordinate {key}.")
        rows.append({"coordinate": key, "lon": lon, "lat": lat, "properties": props})
    return rows


def _flatten(value: Any, prefix: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {}
    if isinstance(value, dict):
        for key in sorted(value):
            name = f"{prefix}.{key}" if prefix else str(key)
            out.update(_flatten(value[key], name))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            name = f"{prefix}[{index}]"
            out.update(_flatten(item, name))
    else:
        out[prefix or "$"] = value
    return out


def _values_equal(a: Any, b: Any, tolerance: float) -> bool:
    if isinstance(a, (int, float)) and not isinstance(a, bool) and isinstance(b, (int, float)) and not isinstance(b, bool):
        if not (math.isfinite(float(a)) and math.isfinite(float(b))):
            return (math.isnan(float(a)) and math.isnan(float(b))) or a == b
        return math.isclose(float(a), float(b), rel_tol=tolerance, abs_tol=tolerance)
    return a == b


def response_counts(rows: list[dict[str, Any]], heights: list[float]) -> dict[str, int]:
    """Count sites with a complete by_height result at each requested height.

    The top-level threshold_height_m field is a current-output eligibility or
    ranking property and is absent from historical GeoJSONs. It is not a stable
    cross-version schema. Count a site only when its by_height map contains the
    requested numeric key and finite volume, area and backwater values.
    """
    counts: dict[str, int] = {}
    required_fields = ("volume_Mm3", "area_km2", "backwater_km")
    for height in heights:
        key = f"{height:g}"
        count = 0
        for row in rows:
            by_height = row["properties"].get("by_height")
            if not isinstance(by_height, dict):
                continue
            result = None
            for raw_key, candidate in by_height.items():
                try:
                    parsed_height = float(raw_key)
                except (TypeError, ValueError, OverflowError):
                    continue
                if math.isfinite(parsed_height) and math.isclose(
                    parsed_height, height, rel_tol=0.0, abs_tol=1e-9
                ):
                    result = candidate
                    break
            if not isinstance(result, dict):
                continue
            complete = True
            for field in required_fields:
                try:
                    value = float(result[field])
                except (KeyError, TypeError, ValueError, OverflowError):
                    complete = False
                    break
                if not math.isfinite(value):
                    complete = False
                    break
            if complete:
                count += 1
        counts[key] = count
    return counts


def distance_m(a: tuple[float, float], b: tuple[float, float]) -> float:
    """Great-circle distance between (lon, lat) WGS84 points, in metres."""
    lon1, lat1 = map(math.radians, a)
    lon2, lat2 = map(math.radians, b)
    dlat, dlon = lat2 - lat1, lon2 - lon1
    hav = (math.sin(dlat / 2.0) ** 2
           + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2.0) ** 2)
    hav = min(1.0, max(0.0, hav))
    return 2.0 * 6_371_008.8 * math.asin(math.sqrt(hav))


def nearest_unmatched(source: dict[tuple[float, float], dict],
                      target: dict[tuple[float, float], dict],
                      unmatched: set[tuple[float, float]],
                      *, target_name: str) -> tuple[list[dict], list[float]]:
    """Provide nearest-site suggestions for review without treating them as matches."""
    suggestions, distances = [], []
    if not target:
        return suggestions, distances
    for coordinate in sorted(unmatched):
        source_row = source[coordinate]
        distance, nearest_coordinate = min(
            (distance_m(coordinate, candidate_coordinate), candidate_coordinate)
            for candidate_coordinate in target
        )
        distances.append(distance)
        suggestion = {
            f"nearest_{target_name}_lon": nearest_coordinate[0],
            f"nearest_{target_name}_lat": nearest_coordinate[1],
            f"nearest_{target_name}_distance_m": round(distance, 3),
            "nearest_neighbour_is_identity_match": False,
        }
        source_row = source_row.copy()
        source_row.update(suggestion)
        suggestions.append(source_row)
    return suggestions, distances


def distance_summary(distances: list[float]) -> dict:
    """Compact distance distribution; it is descriptive, never an identity rule."""
    if not distances:
        return {"count": 0, "min_m": None, "median_m": None, "max_m": None,
                "within_100m": 0, "within_500m": 0, "within_1000m": 0}
    ordered = sorted(distances)
    n = len(ordered)
    median = (ordered[n // 2] if n % 2 else
              (ordered[n // 2 - 1] + ordered[n // 2]) / 2.0)
    return {
        "count": n,
        "min_m": round(ordered[0], 3),
        "median_m": round(median, 3),
        "max_m": round(ordered[-1], 3),
        "within_100m": sum(x <= 100 for x in ordered),
        "within_500m": sum(x <= 500 for x in ordered),
        "within_1000m": sum(x <= 1000 for x in ordered),
    }


def compare_geojson(baseline_path: Path, candidate_path: Path, *,
                    tolerance: float = 1e-9, heights: list[float] | None = None) -> tuple[dict, list[dict]]:
    baseline = _read_geojson(baseline_path)
    candidate = _read_geojson(candidate_path)
    heights = heights or [10.0, 25.0, 50.0, 100.0, 150.0]
    bmap = {row["coordinate"]: row for row in baseline}
    cmap = {row["coordinate"]: row for row in candidate}
    bkeys, ckeys = set(bmap), set(cmap)
    common = sorted(bkeys & ckeys)
    changes = []
    for key in common:
        bf, cf = _flatten(bmap[key]["properties"]), _flatten(cmap[key]["properties"])
        for field in sorted(set(bf) | set(cf)):
            a, b = bf.get(field), cf.get(field)
            if field not in bf or field not in cf or not _values_equal(a, b, tolerance):
                changes.append({
                    "status": "PROPERTY_CHANGED",
                    "lon": key[0], "lat": key[1], "field": field,
                    "baseline": json.dumps(a, sort_keys=True, allow_nan=False) if a is not None else "",
                    "candidate": json.dumps(b, sort_keys=True, allow_nan=False) if b is not None else "",
                })
    baseline_unmatched = bkeys - ckeys
    candidate_unmatched = ckeys - bkeys
    baseline_nearest_rows, baseline_nearest_distances = nearest_unmatched(
        bmap, cmap, baseline_unmatched, target_name="candidate"
    )
    candidate_nearest_rows, candidate_nearest_distances = nearest_unmatched(
        cmap, bmap, candidate_unmatched, target_name="baseline"
    )
    nearest_baseline_by_coord = {row["coordinate"]: row for row in baseline_nearest_rows}
    nearest_candidate_by_coord = {row["coordinate"]: row for row in candidate_nearest_rows}

    for lon, lat in sorted(baseline_unmatched):
        row = {"status": "BASELINE_SITE_MISSING", "lon": lon, "lat": lat,
               "field": "", "baseline": "present", "candidate": "missing"}
        row.update({k: v for k, v in nearest_baseline_by_coord[(lon, lat)].items()
                    if k.startswith("nearest_") or k == "nearest_neighbour_is_identity_match"})
        changes.append(row)
    for lon, lat in sorted(candidate_unmatched):
        row = {"status": "NEW_CANDIDATE_SITE", "lon": lon, "lat": lat,
               "field": "", "baseline": "missing", "candidate": "present"}
        row.update({k: v for k, v in nearest_candidate_by_coord[(lon, lat)].items()
                    if k.startswith("nearest_") or k == "nearest_neighbour_is_identity_match"})
        changes.append(row)

    baseline_count = len(baseline)
    candidate_count = len(candidate)
    summary = {
        "schema_version": 1,
        "status": "EXACT_SEMANTIC_MATCH" if not changes and baseline_count == candidate_count else "REVIEW_REQUIRED",
        "baseline_path": str(baseline_path),
        "candidate_path": str(candidate_path),
        "baseline_sha256": sha256_file(baseline_path),
        "candidate_sha256": sha256_file(candidate_path),
        "baseline_site_count": baseline_count,
        "candidate_site_count": candidate_count,
        "exact_coordinate_matches": len(common),
        "baseline_sites_without_exact_match": len(bkeys - ckeys),
        "candidate_sites_without_exact_match": len(ckeys - bkeys),
        "coordinate_overlap_fraction_of_baseline": len(common) / baseline_count if baseline_count else None,
        "nearest_neighbour_review_only": {
            "baseline_unmatched_to_candidate": distance_summary(baseline_nearest_distances),
            "candidate_unmatched_to_baseline": distance_summary(candidate_nearest_distances),
            "warning": "Nearest neighbours are suggestions for human review only. They are not paired as identical sites, do not change exact_coordinate_matches, and can never produce an EXACT_SEMANTIC_MATCH status."
        },
        "property_change_count_at_exactly_matched_sites": sum(r["status"] == "PROPERTY_CHANGED" for r in changes),
        "review_row_count": len(changes),
        "response_counts_by_threshold_height_m": {
            "baseline": response_counts(baseline, heights),
            "candidate": response_counts(candidate, heights),
        },
        "response_count_definition": (
            "Count a site only when by_height contains the requested numeric height "
            "and a complete finite volume_Mm3/area_km2/backwater_km result. The "
            "top-level threshold_height_m property is not used because it is not "
            "present in historical GeoJSONs."
        ),
        "comparison_tolerance": tolerance,
        "identity_rule": "Only exact WGS84 coordinate pairs are treated as the same site; nearest neighbours are never paired as if identical.",
        "interpretation": "A semantic match is a regression check, not geophysical validation. A difference is flagged for review and does not by itself identify which result is scientifically correct.",
    }
    return summary, changes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True, dest="json_out")
    parser.add_argument("--csv", type=Path, required=True, dest="csv_out")
    parser.add_argument("--tolerance", type=float, default=1e-9)
    parser.add_argument("--heights", type=float, nargs="+", default=[10, 25, 50, 100, 150])
    args = parser.parse_args()
    if not math.isfinite(args.tolerance) or args.tolerance < 0:
        parser.error("--tolerance must be finite and non-negative")
    if not args.heights or any(not math.isfinite(h) or h <= 0 for h in args.heights):
        parser.error("--heights must contain finite positive blockage heights")
    summary, changes = compare_geojson(args.baseline, args.candidate,
                                       tolerance=args.tolerance, heights=args.heights)
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.csv_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    fields = ["status", "lon", "lat", "field", "baseline", "candidate"]
    with args.csv_out.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(changes)
    print(f"IMPOUNDMENT_COMPARE status={summary['status']} baseline_sites={summary['baseline_site_count']} candidate_sites={summary['candidate_site_count']} exact_matches={summary['exact_coordinate_matches']} review_rows={summary['review_row_count']}")
    print(f"BASELINE_H10={summary['response_counts_by_threshold_height_m']['baseline'].get('10',0)} CANDIDATE_H10={summary['response_counts_by_threshold_height_m']['candidate'].get('10',0)}")
    print(f"WROTE_JSON={args.json_out} WROTE_CSV={args.csv_out}")
    # Differences are scientifically important but must be reviewed, not hidden.
    # Workflow remains successful so the diagnostic artifact is retained.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
