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


def complete_height_result(row: dict[str, Any], height: float) -> dict[str, Any] | None:
    """Return a complete finite by-height result, or None if it is not measured."""
    by_height = row.get("properties", {}).get("by_height")
    if not isinstance(by_height, dict):
        return None
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
        return None
    for field in ("volume_Mm3", "area_km2", "backwater_km"):
        try:
            value = float(result[field])
        except (KeyError, TypeError, ValueError, OverflowError):
            return None
        if not math.isfinite(value) or value < 0:
            return None
    return result


def response_counts(rows: list[dict[str, Any]], heights: list[float]) -> dict[str, int]:
    """Count sites with a complete by_height result at each requested height.

    The top-level threshold_height_m field is a current-output ranking property
    and is absent from historical GeoJSONs. It is not a stable cross-version
    schema. A height counts only when its by_height entry contains finite,
    non-negative volume, area and backwater values. This is a result-coverage
    count, not proof that the DEM or the physical dam scenario is accurate.
    """
    return {
        f"{height:g}": sum(complete_height_result(row, height) is not None for row in rows)
        for height in heights
    }


def truncated_response_counts(rows: list[dict[str, Any]], heights: list[float]) -> dict[str, int]:
    """Count complete height results whose modeled water surface hits the grid edge.

    Such volumes are lower bounds. Keep their count separate from the number
    of computed responses so a complete-looking value is not mistaken for a
    closed-basin estimate.
    """
    return {
        f"{height:g}": sum(
            (result := complete_height_result(row, height)) is not None
            and result.get("truncated") is True
            for row in rows
        )
        for height in heights
    }


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


def _normalised_properties(row: dict[str, Any]) -> dict[str, Any]:
    """Normalize only a known schema addition derivable from by_height.

    Old tracked GeoJSONs predate top-level threshold_height_m. Derive it from
    the smallest computed height rather than treating a missing redundant field
    as a science mismatch. Existing values are never overwritten.
    """
    props = dict(row["properties"])
    if "threshold_height_m" not in props and isinstance(props.get("by_height"), dict):
        eligible = []
        for raw_height in props["by_height"]:
            try:
                height = float(raw_height)
            except (TypeError, ValueError, OverflowError):
                continue
            if math.isfinite(height) and complete_height_result(row, height) is not None:
                eligible.append(height)
        if eligible:
            props["threshold_height_m"] = min(eligible)
    return props


def compare_geojson(
    baseline_path: Path,
    candidate_path: Path,
    *,
    tolerance: float = 1e-9,
    coordinate_tolerance_m: float = 2.0,
    heights: list[float] | None = None,
) -> tuple[dict, list[dict]]:
    baseline = _read_geojson(baseline_path)
    candidate = _read_geojson(candidate_path)
    heights = heights or [10.0, 25.0, 50.0, 100.0, 150.0]
    if not math.isfinite(coordinate_tolerance_m) or coordinate_tolerance_m < 0:
        raise ValueError("coordinate_tolerance_m must be finite and non-negative")

    bmap = {row["coordinate"]: row for row in baseline}
    cmap = {row["coordinate"]: row for row in candidate}
    bkeys, ckeys = set(bmap), set(cmap)
    exact_keys = sorted(bkeys & ckeys)
    changes: list[dict] = []
    property_changes_exact = 0
    property_changes_rounded = 0

    # Exact coordinate pairs take precedence.
    matched_baseline = set(exact_keys)
    matched_candidate = set(exact_keys)
    rounded_pairs: list[tuple[tuple[float, float], tuple[float, float], float]] = []

    # Historical GeoJSONs keep full-precision coordinates; the current writer
    # serializes lon/lat to five decimal places. Use a small, explicit 2 m
    # tolerance as a secondary review category, never as an exact identity rule.
    candidate_pairs = []
    for bcoord in sorted(bkeys - matched_baseline):
        for ccoord in sorted(ckeys - matched_candidate):
            distance = distance_m(bcoord, ccoord)
            if distance <= coordinate_tolerance_m:
                candidate_pairs.append((distance, bcoord, ccoord))
    for distance, bcoord, ccoord in sorted(candidate_pairs):
        if bcoord in matched_baseline or ccoord in matched_candidate:
            continue
        matched_baseline.add(bcoord)
        matched_candidate.add(ccoord)
        rounded_pairs.append((bcoord, ccoord, distance))
        changes.append({
            "status": "COORDINATE_ROUNDING_MATCH",
            "lon": bcoord[0], "lat": bcoord[1],
            "candidate_lon": ccoord[0], "candidate_lat": ccoord[1],
            "coordinate_distance_m": round(distance, 3),
            "field": "geometry.coordinates",
            "baseline": json.dumps(list(bcoord)),
            "candidate": json.dumps(list(ccoord)),
            "nearest_neighbour_is_identity_match": False,
        })

    # Compare exact pairs and precision-tolerant pairs with separate labels.
    matched_pairs = [(key, key, "EXACT") for key in exact_keys]
    matched_pairs += [(bcoord, ccoord, "ROUNDED")
                      for bcoord, ccoord, _ in rounded_pairs]
    for bcoord, ccoord, match_kind in matched_pairs:
        bf = _flatten(_normalised_properties(bmap[bcoord]))
        cf = _flatten(_normalised_properties(cmap[ccoord]))
        for field in sorted(set(bf) | set(cf)):
            a, b = bf.get(field), cf.get(field)
            if field not in bf or field not in cf or not _values_equal(a, b, tolerance):
                status = ("PROPERTY_CHANGED_AT_ROUNDED_COORDINATE_MATCH"
                          if match_kind == "ROUNDED" else "PROPERTY_CHANGED")
                changes.append({
                    "status": status,
                    "lon": bcoord[0], "lat": bcoord[1],
                    "candidate_lon": ccoord[0], "candidate_lat": ccoord[1],
                    "coordinate_distance_m": round(distance_m(bcoord, ccoord), 3),
                    "field": field,
                    "baseline": json.dumps(a, sort_keys=True, allow_nan=False) if a is not None else "",
                    "candidate": json.dumps(b, sort_keys=True, allow_nan=False) if b is not None else "",
                    "nearest_neighbour_is_identity_match": match_kind == "EXACT",
                })
                if match_kind == "ROUNDED":
                    property_changes_rounded += 1
                else:
                    property_changes_exact += 1

    baseline_unmatched = bkeys - matched_baseline
    candidate_unmatched = ckeys - matched_candidate
    baseline_target = {key: cmap[key] for key in candidate_unmatched}
    candidate_target = {key: bmap[key] for key in baseline_unmatched}
    baseline_nearest_rows, baseline_nearest_distances = nearest_unmatched(
        bmap, baseline_target, baseline_unmatched, target_name="candidate"
    )
    candidate_nearest_rows, candidate_nearest_distances = nearest_unmatched(
        cmap, candidate_target, candidate_unmatched, target_name="baseline"
    )
    nearest_baseline_by_coord = {row["coordinate"]: row for row in baseline_nearest_rows}
    nearest_candidate_by_coord = {row["coordinate"]: row for row in candidate_nearest_rows}

    for lon, lat in sorted(baseline_unmatched):
        row = {"status": "BASELINE_SITE_MISSING", "lon": lon, "lat": lat,
               "field": "", "baseline": "present", "candidate": "missing",
               "nearest_neighbour_is_identity_match": False}
        suggestion = nearest_baseline_by_coord.get((lon, lat))
        if suggestion:
            row.update({k: v for k, v in suggestion.items()
                        if k.startswith("nearest_") or k == "nearest_neighbour_is_identity_match"})
        changes.append(row)
    for lon, lat in sorted(candidate_unmatched):
        row = {"status": "NEW_CANDIDATE_SITE", "lon": lon, "lat": lat,
               "field": "", "baseline": "missing", "candidate": "present",
               "nearest_neighbour_is_identity_match": False}
        suggestion = nearest_candidate_by_coord.get((lon, lat))
        if suggestion:
            row.update({k: v for k, v in suggestion.items()
                        if k.startswith("nearest_") or k == "nearest_neighbour_is_identity_match"})
        changes.append(row)

    baseline_count = len(baseline)
    candidate_count = len(candidate)
    substantive_changes = [row for row in changes if row["status"] != "COORDINATE_ROUNDING_MATCH"]
    if substantive_changes or baseline_unmatched or candidate_unmatched:
        status = "REVIEW_REQUIRED"
    elif rounded_pairs:
        status = "SEMANTIC_MATCH_WITH_COORDINATE_ROUNDING"
    else:
        status = "EXACT_SEMANTIC_MATCH"

    summary = {
        "schema_version": 2,
        "status": status,
        "baseline_path": str(baseline_path),
        "candidate_path": str(candidate_path),
        "baseline_sha256": sha256_file(baseline_path),
        "candidate_sha256": sha256_file(candidate_path),
        "baseline_site_count": baseline_count,
        "candidate_site_count": candidate_count,
        "exact_coordinate_matches": len(exact_keys),
        "coordinate_rounding_matches": len(rounded_pairs),
        "coordinate_tolerance_m": coordinate_tolerance_m,
        "baseline_sites_without_exact_match": len(bkeys - ckeys),
        "candidate_sites_without_exact_match": len(ckeys - bkeys),
        "baseline_sites_unmatched_after_coordinate_tolerance": len(baseline_unmatched),
        "candidate_sites_unmatched_after_coordinate_tolerance": len(candidate_unmatched),
        "coordinate_overlap_fraction_of_baseline": len(exact_keys) / baseline_count if baseline_count else None,
        "matched_fraction_with_rounding_tolerance": len(matched_baseline) / baseline_count if baseline_count else None,
        "property_change_count_at_exactly_matched_sites": property_changes_exact,
        "property_change_count_at_rounded_coordinate_matches": property_changes_rounded,
        "review_row_count": len(substantive_changes),
        "comparison_row_count": len(changes),
        "response_counts_by_threshold_height_m": {
            "baseline": response_counts(baseline, heights),
            "candidate": response_counts(candidate, heights),
        },
        "truncated_response_counts_by_threshold_height_m": {
            "baseline": truncated_response_counts(baseline, heights),
            "candidate": truncated_response_counts(candidate, heights),
        },
        "response_count_definition": (
            "Count a site only when by_height contains the requested numeric height "
            "and a complete finite non-negative volume_Mm3/area_km2/backwater_km "
            "result. This is computed-result coverage, not independent hydrologic "
            "validation. Truncated counts are reported separately because their "
            "volumes are lower bounds. Missing historical threshold_height_m is "
            "derived from the lowest complete by_height key for comparison."
        ),
        "coordinate_match_policy": (
            "Exact coordinate pairs are exact matches. A one-to-one pair within "
            f"{coordinate_tolerance_m:g} m is labelled COORDINATE_ROUNDING_MATCH, "
            "not exact identity; the threshold accommodates five-decimal coordinate "
            "serialization only. All other nearest neighbours are review hints only."
        ),
        "nearest_neighbour_review_only": {
            "baseline_unmatched_to_candidate": distance_summary(baseline_nearest_distances),
            "candidate_unmatched_to_baseline": distance_summary(candidate_nearest_distances),
            "warning": "Nearest neighbours are suggestions for human review only. They are not paired as identical sites and cannot change semantic status."
        },
        "comparison_tolerance": tolerance,
        "interpretation": (
            "A comparison pass is not geophysical validation. Differences in unmatched "
            "locations, per-height water surfaces, sill, upstream area, volume, truncation "
            "or ranks require review; no result is auto-accepted."
        ),
    }
    return summary, changes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True, dest="json_out")
    parser.add_argument("--csv", type=Path, required=True, dest="csv_out")
    parser.add_argument("--tolerance", type=float, default=1e-9,
                        help="relative/absolute numeric-property tolerance")
    parser.add_argument("--coordinate-tolerance-m", type=float, default=2.0,
                        help="secondary coordinate-rounding review tolerance; nearest neighbours are never exact matches")
    parser.add_argument("--heights", type=float, nargs="+", default=[10, 25, 50, 100, 150])
    args = parser.parse_args()
    if not math.isfinite(args.tolerance) or args.tolerance < 0:
        parser.error("--tolerance must be finite and non-negative")
    if not math.isfinite(args.coordinate_tolerance_m) or args.coordinate_tolerance_m < 0:
        parser.error("--coordinate-tolerance-m must be finite and non-negative")
    if not args.heights or any(not math.isfinite(h) or h <= 0 for h in args.heights):
        parser.error("--heights must contain finite positive blockage heights")
    summary, changes = compare_geojson(args.baseline, args.candidate,
                                       tolerance=args.tolerance,
                                       coordinate_tolerance_m=args.coordinate_tolerance_m,
                                       heights=args.heights)
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.csv_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    fields = sorted({key for row in changes for key in row}) or [
        "status", "lon", "lat", "field", "baseline", "candidate"
    ]
    with args.csv_out.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(changes)
    h10 = summary["response_counts_by_threshold_height_m"]
    t10 = summary["truncated_response_counts_by_threshold_height_m"]
    print(
        f"IMPOUNDMENT_COMPARE status={summary['status']} "
        f"baseline_sites={summary['baseline_site_count']} "
        f"candidate_sites={summary['candidate_site_count']} "
        f"exact_matches={summary['exact_coordinate_matches']} "
        f"rounding_matches={summary['coordinate_rounding_matches']} "
        f"unmatched_baseline={summary['baseline_sites_unmatched_after_coordinate_tolerance']} "
        f"unmatched_candidate={summary['candidate_sites_unmatched_after_coordinate_tolerance']} "
        f"review_rows={summary['review_row_count']}"
    )
    print(
        f"BASELINE_H10={h10['baseline'].get('10',0)} CANDIDATE_H10={h10['candidate'].get('10',0)} "
        f"BASELINE_H10_TRUNCATED={t10['baseline'].get('10',0)} "
        f"CANDIDATE_H10_TRUNCATED={t10['candidate'].get('10',0)}"
    )
    print(f"WROTE_JSON={args.json_out} WROTE_CSV={args.csv_out}")
    # Differences are scientifically important but must be reviewed, not hidden.
    # Workflow remains successful so the diagnostic artifact is retained.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
