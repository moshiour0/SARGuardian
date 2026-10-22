#!/usr/bin/env python3
"""Controlled attribution of the impoundment GeoJSON changes.

Runs the pre-baseline impoundment algorithm and current algorithm against the
same immutable DEM samples from a prior audited run, then compares each result
with the committed GeoJSON baseline. A green execution means the experiment
completed; semantic differences are reported as REVIEW_REQUIRED, not passed.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402
from tools.compare_impoundment_geojson import compare_geojson  # noqa: E402

AOIS = ("source", "langtang", "lhende", "blatten")
HEIGHTS = (10.0, 25.0, 50.0, 100.0, 150.0)
MIN_ACCUM_CELLS = 40
TOP_N = 12
MIN_SEPARATION_KM = 1.5
RANK_BY = "efficiency"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load Python module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def finite_count_by_height(rows: list[dict], heights=HEIGHTS) -> dict[str, int]:
    result = {}
    for height in heights:
        count = 0
        for site in rows:
            by_height = site.get("by_height", {})
            value = next((v for k, v in by_height.items()
                          if math.isclose(float(k), height, rel_tol=0, abs_tol=1e-9)), None)
            if not isinstance(value, dict):
                continue
            try:
                nums = [float(value[key]) for key in ("volume_Mm3", "area_km2", "backwater_km")]
            except (KeyError, TypeError, ValueError, OverflowError):
                continue
            if all(math.isfinite(v) and v >= 0 for v in nums):
                count += 1
        result[f"{height:g}"] = count
    return result


def legacy_efficiency_rank(sites: list[dict], rank_by: str, top_n: int,
                          min_separation_km: float) -> list[dict]:
    """Reproduce the rank/NMS policy from the baseline generator commit.

    This is audit-only. Historical code sorted all sites by efficiency (or
    volume), then applied 1.5-km non-maximum suppression and top-N truncation.
    The function intentionally does not adopt the current threshold-first rank.
    """
    if rank_by not in {"volume", "efficiency"}:
        raise ValueError(f"Unsupported ranking metric: {rank_by}")
    metric = "efficiency_Mm3_per_m" if rank_by == "efficiency" else "max_volume_Mm3"
    ranked = sorted(sites, key=lambda row: -float(row[metric]))
    keep = []
    for row in ranked:
        if any(math.hypot(
            (float(row["lat"]) - float(other["lat"])) * 110.57,
            (float(row["lon"]) - float(other["lon"])) * 111.32
            * math.cos(math.radians(float(row["lat"])))
        ) < min_separation_km for other in keep):
            continue
        keep.append(row)
        if len(keep) >= top_n:
            break
    return keep


def one_model_summary(module, z, lats, lons, label: str, out_dir: Path, aoi: str):
    rows = module.analyse(
        z, lats, lons, list(HEIGHTS), MIN_ACCUM_CELLS, TOP_N,
        min_separation_km=MIN_SEPARATION_KM, rank_by=RANK_BY,
    )
    out = out_dir / f"{aoi}__{label}.geojson"
    module.write_geojson(rows, out)
    return rows, out, {
        "model": label,
        "sites": len(rows),
        "scenario_site_counts_by_height_m": finite_count_by_height(rows),
        "parameters": {
            "heights_m": list(HEIGHTS), "min_accum_cells": MIN_ACCUM_CELLS,
            "top_n": TOP_N, "min_separation_km": MIN_SEPARATION_KM,
            "rank_by": RANK_BY,
        },
        "output_geojson": str(out),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--historical-module", required=True, type=Path)
    parser.add_argument("--dem-dir", required=True, type=Path)
    parser.add_argument("--baseline-dir", type=Path, default=ROOT / "outputs")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--source-run-id", default="39777265576")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    historical = load_module(args.historical_module.resolve(), "impoundment_historical_baseline")
    current = load_module(ROOT / "src" / "impoundment.py", "impoundment_current_model")

    report = {
        "schema_version": 1,
        "status": "REVIEW_REQUIRED",
        "experiment": "same DEM samples + same parameters; compare historical and current hydrology algorithms",
        "historical_source_commit": "298b6b18c53511821850517822a4a3b7b045a1e2",
        "current_source_commit": "filled_by_workflow",
        "dem_source_run_id": str(args.source_run_id),
        "parameters": {
            "heights_m": list(HEIGHTS), "min_accum_cells": MIN_ACCUM_CELLS,
            "top_n": TOP_N, "min_separation_km": MIN_SEPARATION_KM,
            "rank_by": RANK_BY,
        },
        "aois": {},
        "interpretation_guard": (
            "The old code is an attribution control, not a preferred scientific model. "
            "Similarity to historical output does not establish physical accuracy; a "
            "model-change explanation is not a scientific endorsement."
        ),
    }

    for aoi in AOIS:
        dem_path = args.dem_dir / f"{aoi}_srtm30m.npz"
        baseline_path = args.baseline_dir / f"dam_sites_{aoi}.geojson"
        if not dem_path.is_file() or not baseline_path.is_file():
            raise FileNotFoundError(f"Required controlled input is missing: DEM={dem_path}; baseline={baseline_path}")
        with np.load(dem_path, allow_pickle=False) as npz:
            z = np.asarray(npz["z"], dtype=float)
            lats = np.asarray(npz["lats"], dtype=float)
            lons = np.asarray(npz["lons"], dtype=float)
            metadata = json.loads(str(npz["metadata_json"].item()))
        expected_bbox = list(historical.AOIS[aoi])
        if metadata.get("bbox") != expected_bbox and metadata.get("bbox_wgs84") != expected_bbox:
            # Metadata schema has changed over time; fail if the numeric bbox differs.
            observed = metadata.get("bbox", metadata.get("bbox_wgs84"))
            if observed is None or not np.allclose(observed, expected_bbox, rtol=0, atol=1e-8):
                raise ValueError(f"{aoi}: cached DEM bbox {observed} != historical AOI bbox {expected_bbox}")
        if z.shape != (len(lats), len(lons)):
            raise ValueError(f"{aoi}: DEM array shape does not match coordinate axes.")
        logging.info("%s DEM grid shape=%s hash=%s", aoi, z.shape, metadata.get("grid_sha256"))

        hist_rows, hist_path, hist_summary = one_model_summary(
            historical, z.copy(), lats.copy(), lons.copy(), "historical_same_dem", args.out_dir, aoi
        )
        current_rows, current_path, current_summary = one_model_summary(
            current, z.copy(), lats.copy(), lons.copy(), "current_threshold_first_ranking", args.out_dir, aoi
        )
        # Isolate the current version's ranking-policy change without changing
        # current hydrology. The historical rank_sites sorts by volume-per-metre
        # alone; current rank_sites prioritizes minimum tested height first.
        original_rank_sites = current.rank_sites
        try:
            current.rank_sites = legacy_efficiency_rank
            legacy_rank_rows, legacy_rank_path, legacy_rank_summary = one_model_summary(
                current, z.copy(), lats.copy(), lons.copy(),
                "current_hydrology_historical_efficiency_ranking", args.out_dir, aoi
            )
        finally:
            current.rank_sites = original_rank_sites

        hist_vs_baseline, hist_changes = compare_geojson(baseline_path, hist_path)
        curr_legacy_vs_baseline, curr_legacy_changes = compare_geojson(baseline_path, legacy_rank_path)
        curr_vs_baseline, curr_changes = compare_geojson(baseline_path, current_path)
        hist_vs_current, cross_changes = compare_geojson(hist_path, current_path)
        legacy_vs_current, legacy_cross_changes = compare_geojson(legacy_rank_path, current_path)

        report["aois"][aoi] = {
            "dem": {
                "path": str(dem_path),
                "sha256": metadata.get("grid_sha256"),
                "shape_rows_cols": list(z.shape),
                "bbox_wgs84": metadata.get("bbox", metadata.get("bbox_wgs84")),
                "source": metadata.get("source"),
            },
            "baseline_path": str(baseline_path),
            "historical_model": hist_summary,
            "current_hydrology_historical_ranking": legacy_rank_summary,
            "current_model_threshold_first_ranking": current_summary,
            "comparisons": {
                "historical_vs_committed_baseline": hist_vs_baseline,
                "current_hydrology_historical_ranking_vs_committed_baseline": curr_legacy_vs_baseline,
                "current_threshold_first_ranking_vs_committed_baseline": curr_vs_baseline,
                "historical_vs_current_threshold_first_same_dem": hist_vs_current,
                "current_hydrology_historical_ranking_vs_current_threshold_first_same_dem": legacy_vs_current,
            },
            "change_csvs": {
                "historical_vs_baseline": f"{aoi}__historical_vs_baseline.csv",
                "current_historical_ranking_vs_baseline": f"{aoi}__current_historical_ranking_vs_baseline.csv",
                "current_threshold_first_vs_baseline": f"{aoi}__current_threshold_first_vs_baseline.csv",
                "historical_vs_current": f"{aoi}__historical_vs_current.csv",
                "current_historical_ranking_vs_threshold_first": f"{aoi}__current_historical_ranking_vs_threshold_first.csv",
            },
        }
        for name, rows in (
            (f"{aoi}__historical_vs_baseline.csv", hist_changes),
            (f"{aoi}__current_historical_ranking_vs_baseline.csv", curr_legacy_changes),
            (f"{aoi}__current_threshold_first_vs_baseline.csv", curr_changes),
            (f"{aoi}__historical_vs_current.csv", cross_changes),
            (f"{aoi}__current_historical_ranking_vs_threshold_first.csv", legacy_cross_changes),
        ):
            import csv
            keys = sorted({key for row in rows for key in row})
            with (args.out_dir / name).open("w", newline="", encoding="utf-8") as fh:
                writer = csv.DictWriter(fh, fieldnames=keys or ["status"])
                writer.writeheader()
                writer.writerows(rows)

    # Attribution is decisive only if the historical model is materially closer
    # to baseline on the same grid in every AOI; otherwise DEM/grid provenance
    # remains a confounder and this experiment remains REVIEW_REQUIRED.
    closer_count = 0
    ranking_effect_count = 0
    per_aoi_assessment = {}
    for aoi, data in report["aois"].items():
        hist = data["comparisons"]["historical_vs_committed_baseline"]
        curr_legacy = data["comparisons"]["current_hydrology_historical_ranking_vs_committed_baseline"]
        curr = data["comparisons"]["current_threshold_first_ranking_vs_committed_baseline"]
        hist_unmatched = hist["baseline_sites_unmatched_after_coordinate_tolerance"] + hist["candidate_sites_unmatched_after_coordinate_tolerance"]
        legacy_unmatched = curr_legacy["baseline_sites_unmatched_after_coordinate_tolerance"] + curr_legacy["candidate_sites_unmatched_after_coordinate_tolerance"]
        curr_unmatched = curr["baseline_sites_unmatched_after_coordinate_tolerance"] + curr["candidate_sites_unmatched_after_coordinate_tolerance"]
        hist_round = hist["coordinate_rounding_matches"]
        legacy_round = curr_legacy["coordinate_rounding_matches"]
        curr_round = curr["coordinate_rounding_matches"]
        historical_closer = hist_unmatched < curr_unmatched or (
            hist_unmatched == curr_unmatched and hist_round > curr_round
        )
        legacy_rank_closer = legacy_unmatched < curr_unmatched or (
            legacy_unmatched == curr_unmatched and legacy_round > curr_round
        )
        if historical_closer:
            closer_count += 1
        if legacy_rank_closer:
            ranking_effect_count += 1
        per_aoi_assessment[aoi] = {
            "historical_unmatched_site_total": hist_unmatched,
            "current_legacy_ranking_unmatched_site_total": legacy_unmatched,
            "current_threshold_first_ranking_unmatched_site_total": curr_unmatched,
            "historical_rounding_matches": hist_round,
            "current_legacy_ranking_rounding_matches": legacy_round,
            "current_threshold_first_ranking_rounding_matches": curr_round,
            "historical_model_closer_to_baseline": historical_closer,
            "current_hydrology_with_historical_ranking_closer_to_baseline": legacy_rank_closer,
            "ranking_policy_effect_is_material_for_site_identity": legacy_rank_closer,
            "comparison_basis": "unmatched site identities and coordinate-rounding matches only; no nearest-neighbour identity assumptions",
        }
    report["historical_model_closer_to_baseline_aois"] = closer_count
    report["current_hydrology_historical_ranking_closer_to_baseline_aois"] = ranking_effect_count
    report["historical_model_closeness_by_aoi"] = per_aoi_assessment
    if closer_count == len(AOIS) and ranking_effect_count == len(AOIS):
        report["attribution_result"] = "BOTH_HYDROLOGY_AND_RANKING_POLICY_CHANGES_AFFECT_SITE_CATALOGUE"
    elif closer_count == len(AOIS):
        report["attribution_result"] = "MODEL_CHANGE_IS_PRIMARY_EXPLANATION_ON_THE_TESTED_GRIDS; RANKING_EFFECT_QUANTIFIED_SEPARATELY"
    elif closer_count == 0:
        report["attribution_result"] = "DEM_GRID_OR_OTHER_FACTOR_REMAINS_A_MAJOR_CONFOUNDER"
    else:
        report["attribution_result"] = "MIXED_EVIDENCE_REVIEW_REQUIRED"

    report_path = args.out_dir / "historical_model_attribution.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"HISTORICAL_MODEL_ATTRIBUTION={report['attribution_result']}")
    print(f"HISTORICAL_CLOSER_AOIS={closer_count}/{len(AOIS)}")
    for aoi, item in report["historical_model_closeness_by_aoi"].items():
        print(f"AOI={aoi} historical_unmatched={item['historical_unmatched_site_total']} current_unmatched={item['current_unmatched_site_total']} historical_rounding={item['historical_rounding_matches']} current_rounding={item['current_rounding_matches']} historical_closer={item['historical_model_closer_to_baseline']}")
    print(f"WROTE_REPORT={report_path}")
    # A completed experiment may legitimately conclude REVIEW_REQUIRED.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
