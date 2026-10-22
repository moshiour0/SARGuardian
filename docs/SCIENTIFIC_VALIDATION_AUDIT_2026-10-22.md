# SARGuardian scientific validation audit — 2026-10-22

**Isolation:** branch science/valid-pixel-target4-9-audit in moshiour0/SARGuardian. Base was main at 9fcb543ef28a7f6fe70a379aed9cfe3e43e277ec. This work does not modify SARGuardian-backend, its feature branch, PR #35, or either repository's main branch.

**Processing evidence:** SARGuardian-backend Run #55, 39657764262; backend workflow revision e30ce4c6af9a62c194873b658bfbdb5738993408; pinned scientific code revision ab547ecc82bcc4d6ccef545a59852536c2beacc6.

## Executive result

- The 15 baseline CSVs previously compared between the uploaded source ZIP and Run #55 preserve expected logical row keys and scientific numeric/text values within floating-point tolerance; largest numeric difference was about 3.41e-13. Some export_mb file-size metadata changed by 0.001 MB.
- Four newly generated CSVs have no committed reference baseline: ts_gunw_source.csv, geometry_sensitivity_source.csv, detectability_blatten_phase.csv and detectability_blatten_offsets.csv.
- Four impoundment GeoJSON files changed materially. Each contains 12 sites, but exact-coordinate overlap with committed values is only source 4/12, Langtang 2/12, Lhende 1/12 and Blatten 6/12. The number reported as having a 10 m water response is now 12/12 in each AOI, compared with 4/12, 2/12, 1/12 and 6/12 in the reference. Treat as an unresolved scientific change, not formatting.
- Target 9 is outside the exact source AOI polygon by 670 m but inside langtang by 183 m. The 0–1 valid-cell counts in source-clipped outputs cannot be interpreted as native HDF5 lack of signal. Target 4 is inside source AOI but only two of 15 source-clipped GUNW windows have 3 finite cells each; all are below MIN_PX=6. Target 4 has 23–119 finite GOFF cells in each of its 18 current derived-TIFF windows.
- The current Run #55 raw HDF5 is not retained in its derived-only artifact. In this execution environment, Earthdata Login was not configured, the ASF package was not installed, and the Earthdata host could not be resolved. Therefore **raw HDF5 files were not downloaded here**; the branch adds a credential-safe downloader and gate diagnostic so the exact investigation can be run in an authenticated environment. Do not claim root cause for Target 4 or the P05023 mask attrition until that run is done.

## Target 4 and 9: the exact current evidence

A WGS84-to-EPSG:32645 polygon check gives:

| Candidate | In source? | Distance from source boundary | In langtang? | Distance from langtang boundary |
|---|---|---:|---|---:|
| 4 (85.52983E, 28.27799N) | Yes | 2,237.9 m inside | Yes | 3,307.3 m inside |
| 9 (85.47591E, 28.25460N) | No | 670.0 m outside | Yes | 183.1 m inside |

A 13×13 window has radius six 80 m cells (~480 m along each grid axis). Because target 9 lies outside the source polygon, occasional edge cells may cross the polygon, but they are not valid point support. The candidate generator in src/scar_map.py uses a rectangular STAC bbox, not a check against the exact SOURCE_RING. This allows a candidate to be listed for source-AOI analysis even though it is outside the source polygon. The science branch now masks optical NDSI/coverage/clusters with the exact SOURCE_RING, and local_floor.py requires explicit or unambiguously inferred AOI eligibility before it can report a noise floor.

The relevant Run #55 GUNW gates are:

| Product | CRID | Finite nonzero phase | Coherence >= 0.3 | Connected components > 0 | Current water/subswath mask | Source AOI |
|---|---|---:|---:|---:|---:|---:|
| GUNW-006 | X05010 | 9,421,891 | 8,284,468 | 8,157,627 | 8,137,813 | 6,314 |
| GUNW-007 | P05023 | 9,061,090 | 7,460,176 | 7,372,633 | 77,074 | 174 |
| GUNW-009 | P05023 | 8,896,083 | 7,356,248 | 7,258,483 | 79,801 | 363 |

The frame IDs differ (016_010, 016_024, 016_025) so whole-product counts are not a controlled same-footprint comparison and do not prove a maturity-specific mask defect. The sudden gate reduction is still a warning that warrants raw-data diagnosis. Do not remove/relax masks to make more pixels survive.

## Script-by-script validation register

Status labels distinguish software tests from real scientific validation. Tests do not establish geophysical truth; no script has been certified against in-situ deformation observations in this audit.

| Script / component | Scientific role | Important assumption or risk | Status / next validation |
|---|---|---|---|
| src/gunw_reader.py | Unwrapped phase to LOS displacement using -lambda/(4*pi), coherence, connected components, water/subswath mask, AOI clipping and reference subtraction. Optional ionosphere phase-screen subtraction is on by default. | Product-version dataset paths/grid alignment; mask values; screen uncertainty; disconnected phase offsets; reference stability. GUNW documentation says the ionosphere screen is present but not applied to data layers by default. | Unresolved P0 for target support. Added CRID-first maturity logic and regression tests. Raw target-specific gate audit, mask histogram, dataset path/shape, screen-off/on sensitivity and same-component reference checks remain required. |
| src/goff_reader.py | Range/along-track offsets, correlation and SNR gates, reference/deramp and MAD outlier removal. | ASF says GOFF is distributed unculled/not low-pass filtered and may contain outliers. Current MAD keep-mask is computed from range offset and also applied to azimuth. | Unresolved P1. Add separate per-axis diagnostics; test independent range and azimuth outliers; compare correlation/SNR thresholds and layer selection on raw data before changing estimator. |
| src/local_floor.py | MAD-derived per-pixel/window/block noise floor. | MIN_PX=6 is only eligibility, not proof of statistical stability; neighbouring pixels are spatially correlated; clipped target windows can be outside AOI. | Exact AOI guard implemented: point must be inside named polygon; AOI inferred only from unambiguous export-folder suffix or supplied as `--aoi`; CSV writes `TARGET_OUTSIDE_AOI` and `INSUFFICIENT_LOCAL_PIXELS` rows rather than silently omitting them. Spatial-correlation/effective-sample-size validation remains. |
| src/bound.py | Derives an observability bound from measured floor and assumed precursor rate/time window. | Conditional on valid local floor, correct event window and LOS-to-LOS comparisons; sensitivity bound is not proof of no precursor. | Regression-tested construction. Do not report bound when inputs have no valid local estimate. |
| src/candidate_check.py | Candidate geometry/sensitivity gates. | Thresholds, DEM slope/aspect and nominal look-vector assumptions. | Partial. Sensitivity tests exist; check native per-pixel LOS vectors, DEM perturbations, and explicit target-AOI eligibility. |
| src/geometry_merge.py | Projects LOS measurements and merges viewing geometry. | Low sensitivity amplifies errors; nominal tracks may not represent actual product vectors; DEM/aspect and common atmospheric errors propagate. | Partial. Unit tests cover sign/geometry. Add conditioning numbers, per-pixel product look vectors, uncertainty propagation, and null-signal simulations. |
| src/pixel_stack.py | Grid alignment, temporal stack, scene common-mode removal and velocity fit. | Median/common-mode removal can erase broad real motion; missing epochs, autocorrelation and non-uniform temporal spans affect uncertainty. | Partial. Grid regressions exist. Run robustness analysis with/without common-mode removal and autocorrelated residual tests. |
| src/common_ref.py | Establishes consistent relative reference and stable-pixel buffer. | Reference must be in compatible connected component and genuinely stable; an unrelated component adds whole-fringe bias. | Partial. Regression tests exist. In every product verify component ID, reference sample count and post-reference residual. |
| src/troposphere.py | Fits height-associated atmospheric trend and compares before/after noise. | Real deformation can covary with elevation; correction may remove signal; low sample count causes unstable regression. | Unresolved scientifically. Compare raw/corrected series with weather/reanalysis, sample leverage, confidence interval and independent control. |
| src/scar_map.py | Sentinel-2 pre/post NDSI/reflectance change and clusters candidate areas. | Snow/cloud/shadow, seasonal melt, scene registration and bounding-box reads may generate false candidates; clusters are not ground-truth scar polygons. | Exact SOURCE_RING mask implemented for NDSI validity, scene coverage denominators and change clusters. Regression test asserts Target 4 is inside and Target 9 is outside. Optical thresholds, registration and independent event-footprint validation remain. |
| src/impoundment.py | Elevation sampling / local geometry and 10 m water-response flags. | Live OpenTopoData sampling/interpolation and missing elevations change output unless exact grid/coordinates are pinned. | Unresolved P1. Four GeoJSONs changed coordinates/response count. Pin query locations, DEM pixels, CRS/transform, response status and SHA; add semantic baseline tests. |
| src/nisar_acquisition.py | ASF search/download and useful/sanitized error diagnostics. | New releases/reprocessed exact products can change search resolution; search success is not file identity. | Current main API regressions passed in CI. Keep exact basename, CRID, size and SHA-256 in every data manifest; downloader added for Run #55's exact granules. |
| src/paths.py | Resolve product paths and prevent cache/environment traversal. | Multiple files for the same pair can be chosen ambiguously. | Partial. Keep unique basename, CRID, maturity and hash in input manifests. |
| src/organise.py | Groups files by pair date/season/event period. | Organization date is not authoritative product maturity; backdated supplemental Provisional granules exist. | Partial. Reader now prioritizes CRID; manifest/presentation should expose actual CRID and maturity separately. |
| src/detectability.py | Monte Carlo detectability for assumed precursor shapes/noise/cadence. | Detector rates are conditional on chosen simulator, thresholds, seed and noise model, not event probability. | Code/test checked, assumptions not externally calibrated. Publish seeds/parameters; calibrate from independent controls and known positive/null periods. |
| src/impoundment.py / src/geometry_merge.py interaction | Uses DEM-derived terrain and projected signal to infer response geometry. | DEM vertical datum, interpolation and coordinate/CRS choices affect derived values. | Requires local-DEM repeat; live API mismatch still unresolved. |
| demo.py | Orchestrates demo and outputs. | Demo can imply evidence is empirical even when an input/output is synthetic or stale. | Software role. Print exact input hashes, CRIDs and synthetic/real status. |
| docs/figures/make_figures.py | Creates charts/illustrations. | Figures can lag source CSV and imply unvalidated values. | Presentation only. Tie figure generation to output manifests and hashes. |
| docs/site/build_data.py and docs/site/render.py | Populate and render static site. | Site may misrepresent stale data, maturity, uncertainty or missing observations. | Presentation only. Surface no-data vs zero, CRID and last verified output hash. |
| tests/*.py | Synthetic/software regression tests. | Synthetic tests do not prove physical interpretation of real NISAR products. | Latest CI [run 39706771204](https://github.com/moshiour0/SARGuardian/actions/runs/39706771204) passed 196 tests with 17 warnings at code commit `a2e1e009`; it includes the AOI/floor guard regression tests. CI tests software behavior, not physical truth against in-situ displacement. |
| tools/target49_gate_diagnostic.py (new) | Chunked raw-HDF5 per-stage counts, dataset paths/shapes, CRS/grid extents, mask histogram, local target windows and SHA-256. | Reports current-reader mask plus a diagnostic alternative; does not assert either mask is scientifically correct. | Completed raw HDF5 run [39693169804](https://github.com/moshiour0/SARGuardian-backend/actions/runs/39693169804): 12/12 selected exact files, SHA/size verified, source and Langtang AOIs plus GUNW no-ionosphere sensitivity; raw files not retained. A broader reproducibility workflow is running on the isolated backend branch. |
| tools/download_target49_products.py (new) | Downloads exact Run #55 basenames via ASF search, records bytes and SHA-256. | Requires authenticated network, disk and Earthdata login. | Dry-run checked only. No HDF5 was downloaded in this environment. |

## Reproducible Target 4/9 experiment

Seven exact basenames from Run #55 are pinned in tools/download_target49_products.py: one X05010 GUNW control (GUNW-006), two P05023 GUNW granules (GUNW-007/009), one X05010 GOFF control (GOFF-023) and three P05023 GOFF granules (GOFF-025/027/029). GOFF-026 is a descending track and must not be substituted into the ascending track series.

On an internet-enabled host with Earthdata Login configured by the ASF-supported .netrc or ASFSession mechanism, run:
- python tools/download_target49_products.py --profile target4-9 --outdir data/target49_raw
- python tools/target49_gate_diagnostic.py data/target49_raw/*.h5 --aoi source --out outputs/target49_source_aoi.json
- python tools/target49_gate_diagnostic.py data/target49_raw/*.h5 --aoi langtang --out outputs/target49_langtang_aoi.json

Inspect each raw mask's dataset path, dtype, shape, attributes and value histogram. Compare gates in order: native-grid support, finite nonzero phase, coherence, connected components, mask, AOI. For GOFF, report range/azimuth finite values separately plus correlation and SNR gates. Do not remove a gate or lower thresholds merely to get a floor.

## Official references checked 2026-10-22

- ASF/NASA availability and CRID maturity, including supplemental Provisional dates and caution mixing maturity: https://nisar-docs.asf.alaska.edu/availability-overview/
- ASF/NASA Provisional product limitations (updated 2026-07-30): https://nisar-docs.asf.alaska.edu/provisional-known-issues/
- ASF/NASA BETA known issues, including valid-data mask alignment: https://nisar-docs.asf.alaska.edu/product-known-issues/
- GUNW overview (grid and ionosphere screen): https://nisar-docs.asf.alaska.edu/gunw/
- GOFF overview (raw/unculled offsets and possible outliers): https://nisar-docs.asf.alaska.edu/goff/
- ISRO/NISAR L1/L2 product format document, GUNW three-digit mask definition: https://bhoonidhi.nrsc.gov.in/NISAR/NISAR_Data_Product_Format_Document_V1.2.1_digisigned.pdf
- ASF authenticated downloader documentation: https://docs.asf.alaska.edu/asf_search/downloading/

This audit separates software behavior, product documentation and physical validation. Raw HDF5 inspection is complete for 12 selected products. The main unresolved physical question is whether interpolated ionosphere-screen pixels should be excluded or represented with a defensible uncertainty model; GOFF axis/outlier calibration and changed impoundment GeoJSONs remain scientifically unresolved.
