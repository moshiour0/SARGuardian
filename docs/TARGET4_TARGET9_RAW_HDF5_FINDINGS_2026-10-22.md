# Raw NISAR HDF5 investigation: Target 4 and Target 9

**Date:** 2026-10-22  
**Status:** 13 exact Run #55 raw products (7 GUNW, 6 GOFF) were downloaded, SHA-256 and byte sizes were checked, both AOIs were diagnosed, and each raw HDF5 was deleted from the isolated runner immediately afterward. The retained artifact contains diagnostic JSON/CSV and download manifests only, never HDF5.  
**Evidence run:** [authenticated raw audit run #5](https://github.com/moshiour0/SARGuardian-backend/actions/runs/39712719257)  
**Artifact:** `target49-raw-science-diagnostic` (artifact ID 12398535967; 126,850 bytes; 13/13 products, 40/40 diagnostic reports and 128/128 target/layer rows verified; 0 validation errors).  
**Downloaded bytes:** 23,142,072,320 bytes (21.553 GiB) in total across the sequential one-file-at-a-time downloads; this is not retained storage.

## Executive findings

1. **Target 9 is not eligible for analysis against the `source` polygon.** Its centre (85.47591°E, 28.25460°N) is about 670 m outside the exact source AOI, but 183 m inside the `langtang` AOI. In the raw files, Target 9 lies inside the coordinate-grid bounds, but source AOI masking removes its signal samples: GUNW `after_aoi=0`; GOFF may leave 0–1 edge cells, which are not support at the target point. Do not interpret those clipped zeros as absence of raw measurements. Any real Target 9 analysis must explicitly use the `langtang` or a justified target-specific polygon. Langtang may contain event-affected ground and is not automatically a stable control.
2. **Target 4's GUNW support shortage is real in the tested local windows, not merely an AOI-centre error.** The point (85.52983°E, 28.27799°N) is inside both polygons and within the GUNW coordinate-grid bounds for all six sampled files. Across all seven GUNW HDF5s, local windows have 0–10 cells after coherence and connected-component filtering. Five produce zero cells at the connected-component gate (one of these has no finite nonzero phase), one produces only four, and one produces ten before the ionosphere-screen interpolation gate.
3. **The default ionosphere-correction policy is a major additional gate on P05023.** In the one P05023 pair where Target 4 has 10 cells after coherence + connected components + the low-byte land/subswath gate, the default reader excludes pixels flagged as ionosphere-screen interpolated, leaving only 3/169 local cells. The `--no-iono` diagnostic-only sensitivity preserves those 10 cells, but that is a different physical/noise model and is not evidence that the default gate should be removed.
4. **Target 9 has usable raw phase support in three P05023 ascending GUNW pairs when assessed inside `langtang`.** For frame 016_024: 80 cells pass connected components, 29 remain after the default ionosphere-screen gate, and 27 survive the Langtang polygon. For frame 016_025: 80 pass connected components, 31 remain after the default screen gate, and 29 survive the Langtang polygon. For frame 016_026: 83 pass the earlier gates, 37 remain after the screen gate, and 33 survive the Langtang polygon. Other tested GUNW pairs at Target 9 fail at the connected-component gate.
5. **GOFF local sample availability is not the bottleneck in `langtang`.** On the layer-2 offset window, the six tested products provide 100–160 cells for Target 4 and 52–99 for Target 9 after finite range/azimuth, correlation ≥0.3, SNR ≥3, and AOI clipping. This is sample availability, not proof that the raw/unculled offsets are accurate; outlier handling and independent displacement validation remain necessary.

## Target-local GUNW counts (13×13 cells, radius 6)

Each row is an exact raw HDF5 input; counts are cumulative at each gate. A 13×13 window has 169 grid cells. `Connected components > 0` is evaluated after finite non-zero phase and coherence ≥0.3. `Iono` is the additional screen-interpolation mask applied by the current reader when ionosphere correction is enabled.

| Pair / CRID | Target | Finite phase | Coherence ≥0.3 | Connected components >0 | Land/subswath + no anomaly | After ionosphere-screen gate | After chosen AOI |
|---|---:|---:|---:|---:|---:|---:|---:|
| GUNW BETA X05010, frame 016_010 | 4 | 169 | 40 | 0 | 0 | 0 | 0 |
| GUNW BETA X05010, frame 016_010 | 9 | 169 | 51 | 0 | 0 | 0 | 0 |
| GUNW P05023, frame 016_024 | 4 | 169 | 17 | 0 | 0 | 0 | 0 |
| GUNW P05023, frame 016_024 | 9 | 169 | 84 | 80 | 80 | 29 | 27 in Langtang; 0 in source |
| GUNW P05023, frame 016_025 | 4 | 169 | 22 | 0 | 0 | 0 | 0 |
| GUNW P05023, frame 016_025 | 9 | 169 | 80 | 80 | 80 | 31 | 29 in Langtang; 0 in source |
| GUNW P05023, frame 074_026 (descending) | 4 | 0 | 0 | 0 | 0 | 0 | 0 |
| GUNW P05023, frame 074_026 (descending) | 9 | 169 | 16 | 0 | 0 | 0 | 0 |
| GUNW P05023, frame 016_026 | 4 | 169 | 31 | 10 | 10 | 3 | 3 |
| GUNW P05023, frame 016_026 | 9 | 169 | 89 | 83 | 83 | 37 | 33 in Langtang; 0 in source |
| GUNW P05023, frame 074_028 (descending) | 4 | 169 | 66 | 4 | 4 | 3 | 3 |
| GUNW P05023, frame 074_028 (descending) | 9 | 169 | 4 | 0 | 0 | 0 | 0 |
| GUNW P05023, frame 016_028 | 4 | 169 | 35 | 0 | 0 | 0 | 0 |
| GUNW P05023, frame 016_028 | 9 | 169 | 79 | 0 | 0 | 0 | 0 |

**Important interpretation:** In the descending frame 074_026, `inside_grid_extent=true` describes the rectangular coordinate-axis extent; it does not mean the target window contains non-nodata phase. All 169 cells have no finite non-zero phase for Target 4 in this product, so it cannot support this target. Do not conflate grid-bounds inclusion with valid SAR footprint.

## Why the 32-bit P05023 GUNW mask removes so many cells

The raw HDF5 attributes identify the mask at `science/LSAR/GUNW/grids/frequencyA/unwrappedInterferogram/mask` as `uint32` for P05023. Its description states:

- Bits 0–7 encode water plus reference/secondary subswath IDs.
- Bits 8–23 encode secondary/reference anomaly flags.
- Bit 24 marks pixels masked out and filled by interpolation in the ionosphere phase screen.
- Bits 25–31 are reserved.

A frequent raw mask value is decimal `16777228` = hex `0x0100000C`: bit 24 is set while the low byte `0x0C` encodes land and valid nonzero subswaths. This is not a decimal three-digit mask value. The reader correctly must not decode the whole P05023 word as a decimal number.

| P05023 GUNW frame | After connected components | After low-byte land/subswath + anomaly | After default ionosphere-screen interpolation gate | Percent removed at this last gate |
|---|---:|---:|---:|---:|
| 016_024 | 7,372,633 | 7,350,651 | 77,074 | 98.95% |
| 016_025 | 7,258,483 | 7,242,135 | 79,801 | 98.90% |
| 074_026 descending | 6,487,623 | 3,739,028 | 146,004 | 96.09% |
| 016_026 | 7,153,109 | 7,141,049 | 80,156 | 98.88% |
| 074_028 descending | 7,455,180 | 7,445,169 | 71,379 | 99.04% |
| 016_028 | 7,183,273 | 7,173,834 | 112,518 | 98.43% |

This is an important policy/sensitivity result, **not automatic proof of a coding bug**. The bit labels interpolated ionosphere-screen pixels, not necessarily invalid raw phase measurements. The project subtracts the ionosphere screen by default and excludes these pixels because their correction is interpolated. That is a defensible conservative choice, but it can leave insufficient sample support. Before changing this behavior, quantify correction-screen uncertainty or residuals and test a clear alternative model; do not merely drop the mask gate or lower `MIN_PX`.

For contrast, the sampled BETA X05010 GUNW mask is `uint8` using the documented decimal digit layout; it does not carry bit 24. This is why maturity/CRID and actual dataset attributes need to be recorded with every science result.

## GOFF support in the Langtang AOI

Layer-2 target-window cells remaining after finite range + azimuth offsets, correlation ≥0.3, SNR ≥3 and Langtang AOI clipping:

| GOFF frame / CRID | Target 4 | Target 9 |
|---|---:|---:|
| 016_010 / X05010 BETA | 138 | 59 |
| 016_025 / P05023 | 111 | 70 |
| 074_026 descending / P05023 | 136 | 58 |
| 016_026 / P05023 | 100 | 99 |
| 074_028 descending / P05023 | 160 | 52 |
| 016_028 / P05023 | 107 | 91 |

ASF's GOFF documentation says offsets are delivered unculled and not low-pass filtered and may contain outliers. These counts prove local data availability only. The next validation must still test range and azimuth outliers separately, offset units/sign, the selected layer, temporal common-reference behavior and error calibration. The BETA item is not a same-frame matched control for every P05023 item; do not infer that maturity alone caused differences from these mixed frames.

## Inputs downloaded and identity evidence

The 13 HDF5 files were downloaded sequentially with the backend's `tools/download_nisar_product.py`, which records expected response size and SHA-256. Combined bytes downloaded: **23,142,072,320 bytes (21.553 GiB)**. Each file was inspected for both AOIs; each GUNW also has `--no-iono` sensitivity summaries for source and Langtang. Raw HDF5 was deleted from runner scratch after each product was diagnosed. The final artifact contains per-file manifests and JSON/CSV reports, not any source HDF5 bytes. The companion summary has all 128 target/layer/mode rows, including dataset paths/shapes, raw mask histograms, global gates and per-target stage counts.

Source artifact: [authenticated raw-data diagnostic run #5](https://github.com/moshiour0/SARGuardian-backend/actions/runs/39712719257), artifact ID 12398535967. The first isolated downloader attempt used an incomplete ASF exact-name catalogue lookup; the next direct-downloader run exposed and then fixed a shell variable expansion error. The final run completed **13/13 products, 40/40 reports, 128/128 rows, zero validation errors**, and uploaded diagnostics only. These workflow changes remain isolated on `science/target4-9-raw-investigation`; no production workflow or `main` was changed.

## Cross-run source identity against the successful Run #55

I compared each new raw-diagnostic download manifest with the corresponding Run #55 per-product artifact's `download-manifest.json`. **All 13/13 match on exact basename, observed byte count, and SHA-256**. Therefore the gate diagnostics are tied to the exact raw bytes that the successful processing run consumed, rather than merely to a similar granule name or a newer catalog result.

| Product | Observed bytes | SHA-256 |
|---|---:|---|
| GUNW-006 | 2,365,587,456 | `e8b7571c0cb778c9a9e0d8e13d441a685b9bceb508e06cc2e3e16229e1c06392` |
| GUNW-007 | 2,378,170,368 | `dd0e29376f9113060404b6a1b799a3a614f194cca02c8ff80736426699b7a6d2` |
| GUNW-009 | 2,353,004,544 | `35b796f959b0ac76ed5a9b853976530e656aab43d8155274f8fc135825c3af8a` |
| GUNW-010 | 2,336,227,328 | `05dbd5e131a236cf3d141dfadec0a586915400e21862e75e82d4fb5913ea2fa0` |
| GUNW-011 | 2,344,615,936 | `3def328bd568c921486e4c1ed5a6ef7e359db43e99b63eaab083c0393029a499` |
| GUNW-012 | 2,373,976,064 | `2f7b50e7029651cbb96ba58eb4a6c89e8a75106f81abddc2b3a1590f8f609e13` |
| GUNW-013 | 2,390,753,280 | `cae4903019d10bc5568da3e64b7ad9f60bb050c67c947841c3b1104082d949d3` |
| GOFF-023 | 1,115,684,864 | `57818f28c8564372d614eedc373e36fbb917f6fdd042fbe61b7d4eff38db28d8` |
| GOFF-025 | 1,090,519,040 | `8a947f7fa4458afb564431514f36950a89e004d79e6f0e6bfb7461f7e867123c` |
| GOFF-026 | 1,080,033,280 | `541dd2664ad23a18976f8996ee3c230483df6c73d38c6fba0a10bd62b80bf4b7` |
| GOFF-027 | 1,092,616,192 | `ac3eb50eb5c35cb70769b44492398d6b212e880b70558d4d7d2376c5d8e49313` |
| GOFF-028 | 1,119,879,168 | `02af3e9c06b5946dc4c7d42566e0512103e25f7d38e9e782a2f5848b2726eec5` |
| GOFF-029 | 1,101,004,800 | `45e4641da04f235941b1e88262961e9e62ca05ce3ceb35fdac452649b13e416f` |

**Size-field caveat:** the original product manifests sometimes have an `expected_size_bytes` value different from the actual `Content-Length`/downloaded byte count (for example GUNW-006's expected size is 2,400,000,000 bytes but the observed file is 2,365,587,456 bytes). The cross-run identity check deliberately compares the observed `source_bytes`, the exact basename and SHA-256; it does not mistake a catalog/expected-size estimate for the downloaded object's actual size.

## Decision and next steps

- **Target 9:** exact AOI guards are now implemented on the science branch: `scar_map.py` masks NDSI, coverage denominators and clusters with `SOURCE_RING`; `local_floor.py --aoi source` refuses the target centre outside that polygon and records `TARGET_OUTSIDE_AOI`, never a floor. In Langtang, three tested P05023 ascending GUNW pairs have supported phase windows after the current gate, while all six tested GOFF products have layer-2 sample support. Temporal/geometry requirements and stability of the proposed reference still need validation.
- **Target 4:** keep `MIN_PX=6`. Do not claim a GUNW floor from 0–3 pixels. One P05023 pair has 10 pre-ionosphere-gate samples and 3 afterward; another has only 4 connected-component samples; five have zero connected-component support, including one with no finite nonzero phase. Investigate the ionosphere-screen uncertainty/interpolation semantics before altering masks; otherwise report no GUNW floor for the unsupported pairs.
- **Scripts:** AOI eligibility guards and stage-wise HDF5 gate counts have been implemented on the isolated branch. Remaining tasks: add a dedicated GOFF per-axis outlier test and known-offset validation; propagate uncertainty into detectability bounds; pin source/DEM samples and CRS for impoundment outputs; add semantic regression gates for the changed GeoJSON values.
- **Do not merge the scientific branch or change Run #55 output/threshold behavior based on counts alone.** The successful 13-file diagnostic proves the observed counts and locations; it does not certify geophysical truth or settle how interpolated ionosphere-screen pixels should be treated.

## Primary technical references

- [ASF/NASA NISAR GUNW guide](https://nisar-docs.asf.alaska.edu/gunw/) — product layers and ionosphere-screen context.
- [ASF/NASA NISAR GOFF guide](https://nisar-docs.asf.alaska.edu/goff/) — unculled offset caveats.
- [NISAR product availability and CRID maturity](https://nisar-docs.asf.alaska.edu/availability-overview/) — use CRID rather than acquisition dates alone.
- [ASF/NASA Provisional known issues (updated 2026-07-30)](https://nisar-docs.asf.alaska.edu/provisional-known-issues/).
- [ASF/NASA BETA known issues](https://nisar-docs.asf.alaska.edu/product-known-issues/) — BETA valid-data-mask issue is maturity-specific; do not assume it automatically applies to P05023.
- [ISRO/NISAR product-format document](https://bhoonidhi.nrsc.gov.in/NISAR/NISAR_Data_Product_Format_Document_V1.2.1_digisigned.pdf) — authoritative mask bit layout and dataset-level description.
- [ASF authenticated download docs](https://docs.asf.alaska.edu/asf_search/downloading/).
