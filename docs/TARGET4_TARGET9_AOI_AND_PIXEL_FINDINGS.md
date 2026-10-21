# Target 4 / Target 9: AOI and pixel-support findings

Date: 2026-10-22  
Evidence: SARGuardian source snapshot + SARGuardian-backend Run #55 (39657764262) logs and output TIFFs.  
Status: Target 9 AOI mismatch confirmed; Target 4 GUNW support unresolved; raw HDF5 diagnosis pending.

## Main correction

The earlier quick pixel-count interpretation treated source-AOI-clipped TIFFs as if they measured raw source-product validity. That was too strong. Those TIFFs are already clipped by the source AOI. Target 9 is outside that AOI, so its 0/1 finite-cell counts are primarily a spatial eligibility problem and do not establish that the original HDF5 contains no useful values.

A polygon check was performed by projecting each AOI and candidate from EPSG:4326 to EPSG:32645 and measuring distance to polygon boundaries:

| Candidate | Coordinate | source AOI | Distance to source boundary | langtang AOI | Distance to Langtang boundary |
|---|---|---|---:|---|---:|
| 4 | 85.52983E, 28.27799N | inside | 2,237.9 m inside | inside | 3,307.3 m inside |
| 9 | 85.47591E, 28.25460N | outside | 670.0 m outside | inside | 183.1 m inside |

The 13×13 window radius is six 80-m postings, roughly 480 m in each axis. A target centre outside an AOI can still have a few diagonal edge cells cross the polygon boundary. Thus Target 9's rare 1/169 GOFF values after clipping are not point support. The langtang AOI contains event/runout-affected ground and must not be assumed a stable control.

## Current clipped-output counts

Counts below describe existing derived, quality-gated, source-AOI-clipped TIFFs only:

| Target | Family | Outputs | Zero-cell windows | Windows below six pixels | Windows with at least six finite cells | Min–max finite cells per 13×13 window |
|---|---|---:|---:|---:|---:|---:|
| 4 | GUNW | 15 | 13 | 15 | 0 | 0–3 |
| 4 | GOFF | 18 | 0 | 0 | 18 | 23–119 |
| 9 | GUNW | 15 | 15 | 15 | 0 | 0 |
| 9 | GOFF | 18 | 10 | 18 | 0 | 0–1 |

These summaries are not raw-HDF5 validity. The Target 4 GUNW issue remains real within the chosen source polygon, although footprint coverage and mask effects still need separation. For Target 9, rerun with --aoi langtang or a defined target-specific AOI before considering local support.

## Relevant Run #55 GUNW stage counts

| Product | CRID/maturity | Finite, nonzero phase | Coherence >= 0.3 | Connected component > 0 | After current water/subswath mask | After source AOI |
|---|---|---:|---:|---:|---:|---:|
| GUNW-006 | X05010 / BETA | 9,421,891 | 8,284,468 | 8,157,627 | 8,137,813 | 6,314 |
| GUNW-007 | P05023 / PROVISIONAL | 9,061,090 | 7,460,176 | 7,372,633 | 77,074 | 174 |
| GUNW-009 | P05023 / PROVISIONAL | 8,896,083 | 7,356,248 | 7,258,483 | 79,801 | 363 |

These are different frame IDs (016_010, 016_024, 016_025), not same-footprint controls. Whole-product counts alone cannot prove a CRID-specific mask bug. The ISRO/NISAR product format document defines a three-digit UByte mask: the hundreds digit is the reference-RSLC water flag, the tens digit the reference subswath, and the units digit the secondary subswath. The raw HDF5 audit must record the exact selected mask path, shape and histogram before any change to mask logic. ASF documents valid-data mask alignment as a BETA issue; that does not prove the same issue exists in the P05023 products.

## What is and is not proved

Proved from delivered derived TIFFs and logs:
- Target 9 is outside the source AOI but inside langtang.
- Target 4's GUNW local window never reaches MIN_PX=6 in the 15 delivered source-clipped TIFFs; two files have 3 finite cells each.
- Target 4 has sufficient finite GOFF cells in the current derived TIFFs.
- GUNW-007/009 show far fewer cells after the current mask gate than GUNW-006; products have different frame IDs, so this is a warning signal, not proof of mask interpretation error.

Not yet proved:
- Whether raw GUNW HDF5 has usable Target 9 phase/quality values before source-AOI clipping.
- Whether Target 4 is outside the footprint for any pair, lacks local coherence/connected components, is removed by masks, or is affected by other quality filters.
- Whether the very different GUNW global mask counts are due to frame coverage, product mask quality/encoding, sibling dataset selection or combinations.
- Whether Target 4's GOFF offsets are geophysically accurate; finite cells alone are not an accuracy validation.

## Next experiment

1. Download exact Run #55 granules using the new downloader on a host with internet and Earthdata Login configured using ASF's supported .netrc or ASFSession methods. Keep credentials out of the repository and output manifest.
2. Run the raw-HDF5 diagnostic once with --aoi source and once with --aoi langtang. Preserve SHA-256, CRID, selected dataset paths, shapes, EPSG/grid extent, raw mask histograms, local stage counts and AOI membership.
3. Target 4: establish whether its native coordinate lies inside each granule's grid extent before interpreting local stage counts. Compare finite phase, coherence, connected components, mask, then AOI.
4. Target 9: compare pre-AOI native-grid support, then separately show why a source-AOI clip rejects it. Do not pass Target 9's clipped 0/1-cell window to local_floor.py as an instrument floor.
5. Compare X05010 and P05023 only where spatially overlapping, on matched grids/frames. Do not globally compare counts from unmatched frame IDs.
6. Do not change thresholds or scientific outputs until raw HDF5 evidence and the official product format are checked.

## Official technical references

- ASF/NASA NISAR availability and CRID maturity: https://nisar-docs.asf.alaska.edu/availability-overview/
- ASF/NASA Provisional known issues (updated July 30, 2026): https://nisar-docs.asf.alaska.edu/provisional-known-issues/
- ASF/NASA BETA known issues: https://nisar-docs.asf.alaska.edu/product-known-issues/
- GUNW overview: https://nisar-docs.asf.alaska.edu/gunw/
- GOFF overview: https://nisar-docs.asf.alaska.edu/goff/
- ISRO/NISAR product format document, GUNW mask definition: https://bhoonidhi.nrsc.gov.in/NISAR/NISAR_Data_Product_Format_Document_V1.2.1_digisigned.pdf
- ASF download/authentication manual: https://docs.asf.alaska.edu/asf_search/downloading/
