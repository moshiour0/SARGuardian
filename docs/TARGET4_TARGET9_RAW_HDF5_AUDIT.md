# Target 4 / Target 9 — raw NISAR HDF5 support audit

**Scope:** read-only investigation of the exact 13 NISAR granules listed in the Target 4/9 investigation manifest, using the authenticated raw-HDF5 diagnostic workflow in SARGuardian-backend Run #5.

- Workflow: https://github.com/moshiour0/SARGuardian-backend/actions/runs/39712719257
- Source processing run: SARGuardian-backend Run #55, ID 39657764262
- Downloaded and identity-checked: **13/13 granules**, **40/40 diagnostic reports**, **128/128 target/layer rows**, zero diagnostic execution errors.
- Raw input volume streamed: **23,142,072,320 bytes (21.553 GiB)**. SHA-256 and byte size were checked for each input. Raw HDF5 files were deleted after diagnostics; only manifests and small summaries were retained in the artifact.
- This is a pixel-support and mask-gate audit, **not** proof of geophysical accuracy or a field-validated deformation estimate.

## Findings

### Target eligibility

| Target | Coordinate | SOURCE AOI polygon | LANGTANG AOI polygon | Interpretation |
|---|---|---:|---:|---|
| 4 | 28.27799 N, 85.52983 E | Inside | Inside | Eligible for a source-zone point analysis, subject to actual data support |
| 9 | 28.25460 N, 85.47591 E | **Outside** | Inside | Not eligible for a SOURCE-AOI point result. It can be investigated in the broader LANGTANG AOI, but that is a different spatial question. |

The source AOI and Langtang AOI are not interchangeable. In some GOFF source-AOI windows around Target 9, one neighbouring raster cell intersects the polygon even though the target coordinate itself is outside it. A nearby valid cell must not be presented as a measurement at Target 9.

The diagnostic now exposes separate window/point flags: **has_nearby_aoi_gated_support** reports any eligible pixel in the local window; **has_aoi_gated_support** and **has_target_point_support** require both polygon eligibility and the nearest grid cell representing the requested coordinate to pass the product quality gates and AOI mask. A nearby in-polygon cell cannot satisfy the point-level flag. The raw cumulative counts and nearest-cell gate results remain available for transparent review.

### GUNW support at Target 4

The table below uses counts in the 13×13 local window (169 cells) after each quality gate. These are not interchangeable with a local noise-floor estimate; **MIN_PX=6** is a minimum support rule, not a guarantee of independence or physical validity.

| Granule label | Finite nonzero phase | Coherence ≥ 0.3 | Connected components > 0 | Current mask gate | SOURCE-AOI gated |
|---|---:|---:|---:|---:|---:|
| GUNW-006-BETA-control | 169 | 40 | 0 | 0 | 0 |
| GUNW-007-PROVISIONAL | 169 | 17 | 0 | 0 | 0 |
| GUNW-009-PROVISIONAL | 169 | 22 | 0 | 0 | 0 |
| GUNW-010-PROVISIONAL | 0 | 0 | 0 | 0 | 0 |
| GUNW-011-TARGET4-LOW-SUPPORT | 169 | 31 | 10 | 3 | 3 |
| GUNW-012-TARGET4-LOW-SUPPORT | 169 | 66 | 4 | 3 | 3 |
| GUNW-013-PROVISIONAL-FOLLOWUP | 169 | 35 | 0 | 0 | 0 |

**Result:** none of these GUNW windows provides six cells surviving all current gates at Target 4. The strongest current support is only **3 cells** in each of GUNW-011 and GUNW-012. For GUNW-011, the no-ionosphere-gate sensitivity retained 10 cells before the screen-interpolation exclusion and 3 after it; that is a sensitivity comparison, not permission to publish the no-ionosphere result as the corrected product. Do not calculate or export a GUNW local floor/bound at Target 4 from these windows.

### GUNW support at Target 9

Target 9 is outside SOURCE, so **all SOURCE-AOI Target 9 rows are ineligible regardless of neighbouring-window counts**. In the LANGTANG AOI, the final after-AOI counts were:

| Granule label | Final LANGTANG-AOI gated cells |
|---|---:|
| GUNW-006-BETA-control | 0 |
| GUNW-007-PROVISIONAL | 27 |
| GUNW-009-PROVISIONAL | 29 |
| GUNW-010-PROVISIONAL | 0 |
| GUNW-011-TARGET4-LOW-SUPPORT | 33 |
| GUNW-012-TARGET4-LOW-SUPPORT | 0 |
| GUNW-013-PROVISIONAL-FOLLOWUP | 0 |

These counts establish that support is product- and AOI-dependent. They do not establish that the pixels are independent, correctly referenced, or a physically accurate measurement.

### GOFF offset support

GOFF is distributed without the same GUNW unwrapping mask. The diagnostic counts are after the finite-offset, correlation, SNR, and AOI gates. The three numbers per row below are layer 1 / layer 2 / layer 3.

| Granule label | Target 4, SOURCE AOI | Target 9, LANGTANG AOI |
|---|---:|---:|
| GOFF-023-BETA-control | 119 / 138 / 135 | 81 / 59 / 24 |
| GOFF-025-PROVISIONAL | 93 / 111 / 142 | 71 / 70 / 69 |
| GOFF-026-PROVISIONAL | 130 / 136 / 41 | 51 / 58 / 54 |
| GOFF-027-PROVISIONAL | 90 / 100 / 151 | 106 / 99 / 102 |
| GOFF-028-PROVISIONAL-MATCHED-DESCENDING | 158 / 160 / 121 | 47 / 52 / 51 |
| GOFF-029-PROVISIONAL | 93 / 107 / 150 | 93 / 91 / 93 |

These are **usable sample counts, not accuracy scores**. Different frames, look directions, CRIDs, and layers must not be pooled as if they were one matched time series. Target 9 must be analysed against the LANGTANG polygon, not reported as a SOURCE-AOI point result.

## Corrective actions and safeguards

1. **Out-of-AOI point guard:** src/local_floor.py rejects targets outside the selected polygon, emits explicit non-measurement rows, and does not export those rows as inverse-velocity floors.
2. **Nearby-pixel/point distinction:** tools/target49_gate_diagnostic.py distinguishes nearby in-AOI window pixels from valid target-point support. Regression coverage reproduces the edge case where an out-of-AOI point's square window intersects the polygon.
3. **Exact optical polygon:** src/scar_map.py applies the exact source polygon to coverage percentages, NDSI/change masks, and candidate clusters; the bounding box is only used for data discovery and raster I/O.
4. **Fail closed on unsupported floors:** no zero floor is substituted for missing/degenerate samples; GUNW Target 4's three-cell support remains insufficient.
5. **Independent physical validation still required:** matched geometry/CRID selection, spatial correlation/effective sample size, block-bootstrap floor uncertainty, phase-reference stability, and external deformation validation remain open. Do not lower thresholds to manufacture support.

## Validation state

- The raw diagnostic workflow passed and its artifact confirms the counts above.
- The scientific regression workflow passed on the immediately preceding commits. The latest commits that add the target-point/window distinction and repair the planned-path classifier must be checked for their final CI result before this branch is considered validated.
- No changes were made to main, the earlier feature branch, or PR #35. No scientific output has been promoted or overwritten.
