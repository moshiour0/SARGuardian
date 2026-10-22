# SARGuardian scientific validation audit — 2026-10-22

**Isolation:** branch science/valid-pixel-target4-9-audit in moshiour0/SARGuardian. Base was main at 9fcb543ef28a7f6fe70a379aed9cfe3e43e277ec. This work does not modify SARGuardian-backend, its feature branch, PR #35, or either repository's main branch.

**Baseline processing evidence:** SARGuardian-backend Run #55, 39657764262; backend workflow revision e30ce4c6af9a62c194873b658bfbdb5738993408.  
**Raw-data follow-up:** [authenticated Target 4/9 audit run #5](https://github.com/moshiour0/SARGuardian-backend/actions/runs/39712719257), workflow revision 9aa829c2c32b706320d4e5bc24ed710654598142; science diagnostics pinned to 0f9b9328db4ddd14c38180d05e9dad9ca10f64cb; artifact ID 12398535967.

## Executive result

- The 15 baseline CSVs previously compared between the uploaded source ZIP and Run #55 preserve expected logical row keys and scientific numeric/text values within floating-point tolerance; largest numeric difference was about 3.41e-13. Some export_mb file-size metadata changed by 0.001 MB.
- Four newly generated CSVs have no committed reference baseline: ts_gunw_source.csv, geometry_sensitivity_source.csv, detectability_blatten_phase.csv and detectability_blatten_offsets.csv.
- Four impoundment GeoJSON files changed materially. Each contains 12 sites, but exact-coordinate overlap with committed values is only source 4/12, Langtang 2/12, Lhende 1/12 and Blatten 6/12. The number reported as having a 10 m water response is now 12/12 in each AOI, compared with 4/12, 2/12, 1/12 and 6/12 in the reference. Treat as an unresolved scientific change, not formatting.
- Target 9 is outside the exact source AOI by about 670 m but inside Langtang by about 183 m. In the source AOI, clipped 0/1-cell windows are polygon ineligibility, not evidence that the native HDF5 lacks signal. Target 4 is inside the source AOI; across seven exact GUNW files, five local windows have zero cells after connected-component filtering (one has no finite nonzero phase), one has four cells, and one has ten. In that one 10-cell P05023 pair, only three survive the default ionosphere-screen gate. All are below `MIN_PX=6`.
- **Raw HDF5 inspection is now complete for 13 exact products** (7 GUNW, 6 GOFF). The authenticated follow-up downloaded 23,142,072,320 bytes (21.553 GiB), validated 40/40 JSON reports and 128/128 target/layer rows against per-file identity, byte-size and SHA-256, with zero validation errors. HDF5 files were deleted from the ephemeral runner after each product; artifact 12398535967 contains diagnostic JSON/CSV/manifests, not raw HDF5. [Run #5](https://github.com/moshiour0/SARGuardian-backend/actions/runs/39712719257).

## Target 4 and 9: the exact current evidence

A WGS84-to-EPSG:32645 polygon check gives:

| Candidate | In source? | Distance from source boundary | In langtang? | Distance from langtang boundary |
|---|---|---:|---|---:|
| 4 (85.52983E, 28.27799N) | Yes | 2,237.9 m inside | Yes | 3,307.3 m inside |
| 9 (85.47591E, 28.25460N) | No | 670.0 m outside | Yes | 183.1 m inside |

A 13×13 window has radius six 80 m cells (~480 m along each grid axis). Because target 9 lies outside the source polygon, occasional edge cells may cross the polygon, but they are not valid point support. The candidate generator in src/scar_map.py uses a rectangular STAC bbox, not a check against the exact SOURCE_RING. This allows a candidate to be listed for source-AOI analysis even though it is outside the source polygon. The science branch now masks optical NDSI/coverage/clusters with the exact SOURCE_RING, and local_floor.py requires explicit or unambiguously inferred AOI eligibility before it can report a noise floor.

The raw HDF5 audit confirms that the selected mask dataset is `science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/mask` for these products, but the representation differs by release. X05010 BETA is `uint8`; P05023 PROVISIONAL is `uint32`. P05023 uses low-byte water/subswath digits, bits 8–23 for anomaly flags, bit 24 for ionosphere-screen interpolation, and bits 25–31 reserved. The current reader excludes bit 24 when screen correction is applied.

| Product | CRID / dtype | Connected components > 0 | After low-byte land/subswath + anomaly | After bit-24 ionosphere-screen gate | After source AOI |
|---|---|---:|---:|---:|---:|
| GUNW-006 | X05010 / uint8 | 8,157,627 | 8,137,813 | 8,137,813 | 6,314 |
| GUNW-007 | P05023 / uint32 | 7,372,633 | 7,350,651 | 77,074 | 174 |
| GUNW-009 | P05023 / uint32 | 7,258,483 | 7,242,135 | 79,801 | 363 |
| GUNW-010 (descending) | P05023 / uint32 | 6,487,623 | 3,739,028 | 146,004 | 1,444 |
| GUNW-011 | P05023 / uint32 | 7,153,109 | 7,141,049 | 80,156 | 154 |
| GUNW-012 (descending) | P05023 / uint32 | 7,455,180 | 7,445,169 | 71,379 | 387 |
| GUNW-013 | P05023 / uint32 | 7,183,273 | 7,173,834 | 112,518 | 389 |

The source AOI denominator is small compared with the native GUNW grid; its final counts are AOI coverage, not global product validity. P05023's bit-24 interpolation gate removes about 96.1–99.0% of pixels remaining after the prior gates. This is confirmed as the stage of greatest attrition, but **not** proof that the bit should be ignored: the flag marks pixels for which the correction screen was interpolated, and the conservative current policy excludes them. Changing this policy requires an uncertainty/residual study.

The seven products are not a matched same-footprint maturity experiment, so global counts do not prove a CRID-specific defect. At Target 4, five products have zero local support after connected components (one has no finite nonzero phase), one has four cells, and one has ten cells before the ionosphere gate (three after it). At Target 9, the source polygon masks out the candidate; in Langtang, three ascending P05023 GUNW pairs leave 27, 29 and 33 cells after default gates. Six GOFF products have layer-2 support after finite/correlation/SNR/AOI gates of 100–160 cells for Target 4 and 52–99 for Target 9 in Langtang. These counts show measurement availability, not geophysical accuracy.


## Script-by-script validation register

Status labels distinguish software tests from real scientific validation. Tests do not establish geophysical truth; no script has been certified against in-situ deformation observations in this audit.

| Script / component | Scientific role | Important assumption or risk | Status / next validation |
|---|---|---|---|
| src/gunw_reader.py | Unwrapped phase to LOS displacement using -lambda/(4*pi), coherence, connected components, water/subswath mask, AOI clipping and reference subtraction. Optional ionosphere phase-screen subtraction is on by default. | Product-version dataset paths/grid alignment; mask values; screen uncertainty; disconnected phase offsets; reference stability. GUNW documentation says the ionosphere screen is present but not applied to data layers by default. | Unresolved P0 for target support. Added CRID-first maturity logic and regression tests. Raw target-specific gate audit, mask histogram, dataset path/shape, screen-off/on sensitivity and same-component reference checks remain required. |
| src/goff_reader.py | Range/along-track offsets, correlation and SNR gates, reference/deramp and MAD outlier removal. | ASF says GOFF is distributed unculled/not low-pass filtered and may contain outliers. Range and azimuth MAD masks are now computed independently and intersected; the shared deramp fit uses common inliers. | Regression tests pass, but no real matched-scene estimator validation yet. Compare residuals and axis-specific rejection against independent offsets before accepting changed output. |
| src/local_floor.py | MAD-derived per-pixel/window/block noise floor. | MIN_PX=6 is only eligibility, not proof of statistical stability; neighbouring pixels are spatially correlated; clipped target windows can be outside AOI. | Exact AOI guard implemented: point must be inside named polygon; AOI inferred only from unambiguous export-folder suffix or supplied as `--aoi`; CSV writes `TARGET_OUTSIDE_AOI` and `INSUFFICIENT_LOCAL_PIXELS` rows rather than silently omitting them. Spatial-correlation/effective-sample-size validation remains. |
| src/bound.py | Derives an observability bound from measured floor and assumed precursor rate/time window. | Conditional on valid local floor, correct event window and LOS-to-LOS comparisons; sensitivity bound is not proof of no precursor. | Regression-tested construction. Do not report bound when inputs have no valid local estimate. |
| src/candidate_check.py | Candidate geometry/sensitivity gates. | Thresholds, DEM slope/aspect and nominal look-vector assumptions. | Partial. Sensitivity tests exist; check native per-pixel LOS vectors, DEM perturbations, and explicit target-AOI eligibility. |
| src/geometry_merge.py | Projects LOS measurements and merges viewing geometry. | Low sensitivity amplifies errors; nominal tracks may not represent actual product vectors; DEM/aspect and common atmospheric errors propagate. | Partial. Unit tests cover sign/geometry. Add conditioning numbers, per-pixel product look vectors, uncertainty propagation, and null-signal simulations. |
| src/pixel_stack.py | Grid alignment, temporal stack, scene common-mode removal and velocity fit. | Median/common-mode removal can erase broad real motion; missing epochs, autocorrelation and non-uniform temporal spans affect uncertainty. | Partial. Grid regressions exist. Run robustness analysis with/without common-mode removal and autocorrelated residual tests. |
| src/common_ref.py | Establishes consistent relative reference and stable-pixel buffer. | Reference must be in compatible connected component and genuinely stable; an unrelated component adds whole-fringe bias. | Partial. Regression tests exist. In every product verify component ID, reference sample count and post-reference residual. |
| src/troposphere.py | Fits height-associated atmospheric trend and compares before/after noise. | Real deformation can covary with elevation; correction may remove signal; low sample count causes unstable regression. | Unresolved scientifically. Compare raw/corrected series with weather/reanalysis, sample leverage, confidence interval and independent control. |
| src/scar_map.py | Sentinel-2 pre/post NDSI/reflectance change and clusters candidate areas. | Snow/cloud/shadow, seasonal melt, scene registration and bounding-box reads may generate false candidates; clusters are not ground-truth scar polygons. | Exact SOURCE_RING mask implemented for NDSI validity, scene coverage denominators and change clusters. Regression test asserts Target 4 is inside and Target 9 is outside. Optical thresholds, registration and independent event-footprint validation remain. |
| src/impoundment.py | Elevation sampling / local geometry and 10 m water-response flags. | Live OpenTopoData sampling/interpolation and missing elevations change output unless exact grid/coordinates are pinned. | Reproducibility safeguards added; historical outputs remain scientifically unresolved. Four GeoJSONs still differ from Run #55 references. Explicit --api-cache now pins exact grid bytes and records SHA/provenance; non-EPSG:4326 local DEMs fail closed; no existing GeoJSON outputs were regenerated. |
| src/nisar_acquisition.py | ASF search/download and useful/sanitized error diagnostics. | New releases/reprocessed exact products can change search resolution; search success is not file identity. | Current main API regressions passed in CI. Keep exact basename, CRID, size and SHA-256 in every data manifest; downloader added for Run #55's exact granules. |
| src/paths.py | Resolve product paths and prevent cache/environment traversal. | Multiple files for the same pair can be chosen ambiguously. | Partial. Keep unique basename, CRID, maturity and hash in input manifests. |
| src/organise.py | Groups files by pair date/season/event period. | Organization date is not authoritative product maturity; backdated supplemental Provisional granules exist. | Partial. Reader now prioritizes CRID; manifest/presentation should expose actual CRID and maturity separately. |
| src/detectability.py | Monte Carlo detectability for assumed precursor shapes/noise/cadence. | Detector rates are conditional on chosen simulator, thresholds, seed and noise model, not event probability. | Code/test checked, assumptions not externally calibrated. Publish seeds/parameters; calibrate from independent controls and known positive/null periods. |
| src/impoundment.py / src/geometry_merge.py interaction | Uses DEM-derived terrain and projected signal to infer response geometry. | DEM vertical datum, interpolation and coordinate/CRS choices affect derived values. | Local DEM handling now has CRS and bbox-crop guards; live API input-grid mismatch remains unresolved until cached-grid reruns and per-site review. |
| demo.py | Orchestrates demo and outputs. | Demo can imply evidence is empirical even when an input/output is synthetic or stale. | Software role. Print exact input hashes, CRIDs and synthetic/real status. |
| docs/figures/make_figures.py | Creates charts/illustrations. | Figures can lag source CSV and imply unvalidated values. | Presentation only. Tie figure generation to output manifests and hashes. |
| docs/site/build_data.py and docs/site/render.py | Populate and render static site. | Site may misrepresent stale data, maturity, uncertainty or missing observations. | Presentation only. Surface no-data vs zero, CRID and last verified output hash. |
| tests/*.py | Synthetic/software regression tests. | Synthetic tests do not prove physical interpretation of real NISAR products. | Latest completed CI [run 39714565560](https://github.com/moshiour0/SARGuardian/actions/runs/39714565560) passed 203 tests with 21 warnings at commit `f30a0bff`, including exact AOI, common-reference grid and site-bundle checks. CI tests software behavior, not geophysical truth. |
| tools/target49_gate_diagnostic.py | Chunked raw-HDF5 per-stage counts, dataset paths/shapes, CRS/grid extent, mask histogram, local target windows and SHA-256. | Diagnostic no-ionosphere mode is sensitivity analysis, not an alternative production measurement. | Completed [follow-up run #5](https://github.com/moshiour0/SARGuardian-backend/actions/runs/39712719257): 13/13 files, 40/40 reports, 128/128 rows and 0 validation errors. |
| tools/download_target49_products.py | Downloads exact Run #55 basenames through a canonical HTTPS URL with timeout, content-length and SHA-256 checks. | Requires authenticated network, disk and Earthdata login; no signed URL or credentials may be retained. | Used successfully by backend workflow; raw HDF5 was removed after each diagnostic and never entered the artifact. |
| tools/refresh_site_polygon_guard.py and site-polygon-refresh.yml | Exact-mask the browser PNG; regenerate scene cloud statistics and cluster classification; filter source-bound candidates; rewrite and validate static HTML/data atomically. | Scene asset and SRTM service availability is external; a failed refresh must not publish a partially updated bundle. | [Site refresh run](https://github.com/moshiour0/SARGuardian/actions/runs/39714388544) passed; committed static bundle has polygon metadata/mask, no out-of-AOI source candidates, and HTML matches JSON. |

## Remaining scientific validation

1. **Ionosphere-screen uncertainty:** quantify interpolation uncertainty/residuals and compare the corrected/current gate against the no-ionosphere diagnostic without treating them as equivalent measurements. Do not disable bit 24 solely to increase Target 4 pixel count.
2. **Target 4:** maintain `MIN_PX=6`; do not publish a local GUNW floor when 0–4 cells survive. Check reference pixels against the same connected component and run a known-displacement LOS sign/scale test.
3. **Target 9:** do not compute a source-zone floor/bound for this point. If analyzing the Langtang corridor, state that the AOI is event/runout affected and is not automatically a stable control.
4. **GOFF estimator:** independent per-axis outlier regression tests pass; next validate offset sign/units, correlations/SNR, residuals, and bias against matched scenes or independent references.
5. **Impoundment outputs:** four GeoJSONs still differ materially from committed baselines. New explicit DEM-grid caching and provenance do not retroactively reproduce the missing historic API samples; run one pinned cache per AOI and review all site-level shifts before regenerating outputs.
6. **Unmatched controls:** compare BETA and P05023 only on overlapping/matched frames; current product-wide counts cannot isolate maturity effects.
7. **Static site:** bundle [refresh run](https://github.com/moshiour0/SARGuardian/actions/runs/39714388544) regenerated from the exact polygon; re-run the refresh workflow whenever AOI geometry, scene data or classifier thresholds change.

The evidence establishes software behavior and exact NISAR pixel-gate counts. It does **not** settle whether interpolated ionosphere-screen samples should be hard-excluded, nor certify physical deformation against in-situ measurements.

## Official references checked 2026-10-22

- ASF/NASA availability and CRID maturity, including supplemental Provisional dates and caution mixing maturity: https://nisar-docs.asf.alaska.edu/availability-overview/
- ASF/NASA Provisional product limitations (updated 2026-07-30): https://nisar-docs.asf.alaska.edu/provisional-known-issues/
- ASF/NASA BETA known issues, including valid-data mask alignment: https://nisar-docs.asf.alaska.edu/product-known-issues/
- GUNW overview (grid and ionosphere screen): https://nisar-docs.asf.alaska.edu/gunw/
- GOFF overview (raw/unculled offsets and possible outliers): https://nisar-docs.asf.alaska.edu/goff/
- ISRO/NISAR L1/L2 product format document, GUNW three-digit mask definition: https://bhoonidhi.nrsc.gov.in/NISAR/NISAR_Data_Product_Format_Document_V1.2.1_digisigned.pdf
- ASF authenticated downloader documentation: https://docs.asf.alaska.edu/asf_search/downloading/

This audit separates software behavior, product documentation and physical validation. Raw HDF5 inspection is complete for 13 exact products. The main unresolved physical question is whether interpolated ionosphere-screen pixels should be excluded or represented with a defensible uncertainty model; GOFF axis/outlier calibration and changed impoundment GeoJSONs remain scientifically unresolved.
