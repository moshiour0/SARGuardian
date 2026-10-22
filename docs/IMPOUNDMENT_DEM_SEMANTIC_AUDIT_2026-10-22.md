# Impoundment DEM semantic audit — 22 October 2026

**Disposition: REVIEW REQUIRED; do not replace committed GeoJSON outputs.**  
The new repeatability audit completed successfully, but all four regenerated GeoJSON datasets differ materially from the committed references. The workflow was changed so that reproducible execution and scientific semantic parity are no longer represented by the same green status: unresolved semantic differences now fail a dedicated review gate while preserving the artifact for inspection.

Evidence:
- [Audit run #28](https://github.com/moshiour0/SARGuardian/actions/runs/39783389849)
- Artifact: `impoundment-semantic-audit-39783389849`, ID `12427699540`
- Current isolated branch: `science/valid-pixel-target4-9-audit`
- The audit generated each AOI twice from the same saved elevation grid, identical parameters and pinned implementation; the two GeoJSON outputs and their provenance sidecars were byte-identical for each AOI.

## Results against committed reference

| AOI | Baseline / regenerated sites | Exact coordinate matches | Rounding-only coordinate pairs | Unmatched baseline / regenerated sites | Complete H=10 scenarios, baseline → regenerated | Nearest unmatched baseline-to-current median / max distance | Semantic status |
|---|---:|---:|---:|---:|---:|---:|---|
| source | 12 / 12 | 0 | 1 | 11 / 11 | 4 → 12 | 769 m / 3,248 m | REVIEW_REQUIRED |
| langtang | 12 / 12 | 0 | 0 | 12 / 12 | 2 → 12 | 658 m / 3,841 m | REVIEW_REQUIRED |
| lhende | 12 / 12 | 0 | 0 | 12 / 12 | 1 → 12 | 1,533 m / 4,038 m | REVIEW_REQUIRED |
| blatten | 12 / 12 | 0 | 2 | 10 / 10 | 6 → 12 | 504 m / 2,509 m | REVIEW_REQUIRED |

Nearest-neighbour distances are **review hints only**; they do not establish site identity. Site identity in the comparison is exact coordinate equality, with a narrowly defined coordinate-rounding exception. Different channels, saddles or grid cells can be geographically near but hydrologically distinct.

All four regenerated collections contain 12 returned features and are produced with `top=12`, `rank_by=efficiency`, and `min_separation_km=1.5`. These are ranked/spatially-suppressed candidates, not a census of all channel sites. Counts at imposed water heights describe complete modelled scenarios among those returned features only. They do not establish that a real landslide can construct the blockage or that a real impoundment would behave like this model.

## DEM sample grid is not a 30 m analysis grid

The elevation source is OpenTopoData's SRTM 30 m service, but this experiment samples point elevations on a much sparser regular analysis grid:

| AOI | Stored grid shape | Approximate east-west posting | Approximate north-south posting | Nodata samples |
|---|---:|---:|---:|---:|
| source | 126 × 95 | 96 m | 95 m | 0 |
| langtang | 104 × 115 | 91 m | 91 m | 0 |
| lhende | 99 × 121 | 147 m | 147 m | 0 |
| blatten | 96 × 125 | 93 m | 93 m | 0 |

The `*.npz` cache records the exact sampled elevations, coordinate axes, requested bounding box, endpoint, grid shape and a SHA-256 fingerprint. The reproducibility claim applies to **these saved point elevations and this analysis grid**, not to a native 30 m raster. A nominal 30 m input model interpolated/sampled onto roughly 90–147 m cells cannot support a blanket claim of 30 m hydrologic resolution.

## What the result does—and does not—prove

**Proved:** the current algorithm is deterministic for a fixed saved DEM sample grid and fixed parameters. The four regenerated files were byte-identical to their repeated runs, including provenance. The exact sampled grid hashes and implementation hashes are in the artifact.

**Not proved:** why the regenerated sites differ from the committed reference. The committed reference GeoJSONs do not carry the old elevation-grid samples or a comparable DEM-grid hash, so this experiment cannot distinguish among a change in the original elevation sample locations, changes in hydrology implementation/ranking, and differences in earlier execution inputs. The current provenance captures the *new* sample; it cannot retroactively recover the old one. That missing baseline is an important reproducibility gap.

The current recorded algorithm is `priority-flood-d8-flat-routing-tree-connected-ponding-v2`; it uses depression filling for routing, acyclic equal-height routing, Kahn-ordered flow accumulation, and an upstream receiver-tree-connected ponding model. The provenance records no lateral-spill solver and describes backwater extent using axis-aligned grid bounds. Thus the output is a **screening/scenario model**, not a validated hydraulic inundation simulation. The algorithm and assumptions should be compared against hydrology literature and independent terrain/control cases before interpreting volumes as realistic.

## Required disposition

1. Preserve the committed GeoJSON references until the old DEM inputs or an independently approved expected-site dataset are available; do not overwrite them just to remove a diff.
2. Retain the current cached DEM grids and hashes as one named, reproducible input version. Keep source-product nominal resolution separate from analysis-grid posting in all provenance and displays.
3. Run exact per-site/parameter comparison against an approved reference, including sill elevation, catchment/pond volume by each height, truncation flags, connected-component membership and backwater metrics—not just feature counts or coordinates.
4. Validate routing and ponding on synthetic basins with known answers and on an independently checked DEM; include flats, pits, boundary outlets, narrow channels, nodata gaps, saddles and disconnected depressions.
5. Explicitly report imposed dam height, DEM and grid resolution, truncation/lower-bound status, ranking/suppression settings, no lateral spill solver, and the distinction between a modelled scenario and observed physical blockage.
6. Keep the science PR draft and unmerged. The blocking audit gate now ensures future workflow runs cannot signal semantic parity while any of the four AOIs still require human review.

**Authoritative implementation/data references:** [OpenTopoData API](https://www.opentopodata.org/api/), [SRTM mission overview (NASA/JPL)](https://www.earthdata.nasa.gov/data/instruments/srtm), and the project-local `src/impoundment.py` provenance fields and comparator contract in `tools/compare_impoundment_geojson.py`.