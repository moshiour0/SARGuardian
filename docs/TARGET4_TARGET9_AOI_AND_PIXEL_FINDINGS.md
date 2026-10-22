# Target 4 / Target 9: AOI and pixel-support findings

Date: 2026-10-22  
Evidence: SARGuardian source snapshot + SARGuardian-backend Run #55 (39657764262) logs and output TIFFs.  
Status: raw HDF5 diagnosis completed for 12 exact Run #55 products (6 GUNW, 6 GOFF); Target 9 source-AOI exclusion confirmed; Target 4 GUNW sample support remains below the six-cell floor. See [the full raw HDF5 findings](TARGET4_TARGET9_RAW_HDF5_FINDINGS_2026-10-22.md).

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

These summaries are from derived, source-AOI-clipped TIFFs. The raw HDF5 run now separates their causes: Target 4 GUNW support remains scarce after coherence, connected-components and the current ionosphere-screen mask; Target 9's source-AOI outputs cannot represent its point because the point is outside that polygon. Raw per-stage counts are in [the full findings report](TARGET4_TARGET9_RAW_HDF5_FINDINGS_2026-10-22.md).

## Relevant Run #55 GUNW stage counts

| Product | CRID/maturity | Finite, nonzero phase | Coherence >= 0.3 | Connected component > 0 | After current water/subswath mask | After source AOI |
|---|---|---:|---:|---:|---:|---:|
| GUNW-006 | X05010 / BETA | 9,421,891 | 8,284,468 | 8,157,627 | 8,137,813 | 6,314 |
| GUNW-007 | P05023 / PROVISIONAL | 9,061,090 | 7,460,176 | 7,372,633 | 77,074 | 174 |
| GUNW-009 | P05023 / PROVISIONAL | 8,896,083 | 7,356,248 | 7,258,483 | 79,801 | 363 |

These are different frame IDs (016_010, 016_024, 016_025), not same-footprint controls. Raw HDF5 attributes confirm that BETA X05010 uses a uint8 decimal-digit mask while P05023 uses a uint32 packed bit mask: bits 0–7 water/subswath, bits 8–23 anomaly, bit 24 interpolated ionosphere screen, and bits 25–31 reserved. The decoder must use CRID/dtype-aware bit logic. The current major P05023 count drop occurs at the bit-24 ionosphere-screen gate, not the reserved-bit gate; see the full raw findings report before making any mask change.

## What is and is not proved

Proved by the 12-file raw HDF5 audit (workflow run 39693169804):
- Target 9 is outside source AOI by about 670 m but inside Langtang by about 183 m. In source runs, AOI clipping removes target support; in Langtang, two P05023 GUNW pairs retain 29 and 33 cells after the current default gates.
- Target 4 lies in both polygons. Among six tested GUNW files, four local windows have zero connected-component support or no finite phase; another has only four connected-component pixels. Two P05023 pairs have 10 pixels before the interpolated-ionosphere-screen gate and only 3 after it.
- The P05023 GUNW bit-24 screen mask removes roughly 96.1–99.0% of the phase samples remaining after prior gates in the five tested P05023 files. This is a major scientific policy decision to validate; it is not automatic proof of a software defect.
- Six raw GOFF products have many layer-2 cells after correlation/SNR gates in Langtang: Target 4 has 100–160 cells and Target 9 has 52–99 cells. That does not establish offset accuracy.

Still not proved:
- Whether the interpolated ionosphere-screen locations should be hard-excluded from corrected phase, retained with a defensible uncertainty model, or treated with a different estimator.
- Whether GOFF-derived local displacement noise floors are well-calibrated, because finite/correlation/SNR-valid offsets can still contain outliers and spatially correlated errors.
- Geophysical validity against independent reference displacement or stable control data; BETA and P05023 samples are not all matched frame/track controls.

## Next experiments

1. Investigate the P05023 ionosphere screen's spatial interpolation and uncertainty/residuals; compare corrected/default and uncorrected diagnostic models without treating them as interchangeable measurements.
2. Keep `MIN_PX=6`; do not report a GUNW floor where 0–3 valid samples remain. Add candidate-centre AOI membership as a hard precondition before any local-floor calculation.
3. For Target 9, use a scientifically justified candidate AOI rather than the source polygon. Document whether that ground is stable or event-affected before using it as a control.
4. For GOFF, independently test range and azimuth outlier handling, layer selection, offset sign/units, temporal reference behavior, and noise calibration.
5. Extend the sample to matched overlapping frames for maturity comparisons. Do not attribute global differences to CRID alone if the footprints differ.
6. Do not change published science outputs or merge this branch until these checks are complete.

## Official technical references

- ASF/NASA NISAR availability and CRID maturity: https://nisar-docs.asf.alaska.edu/availability-overview/
- ASF/NASA Provisional known issues (updated July 30, 2026): https://nisar-docs.asf.alaska.edu/provisional-known-issues/
- ASF/NASA BETA known issues: https://nisar-docs.asf.alaska.edu/product-known-issues/
- GUNW overview: https://nisar-docs.asf.alaska.edu/gunw/
- GOFF overview: https://nisar-docs.asf.alaska.edu/goff/
- ISRO/NISAR product format document, GUNW mask definition: https://bhoonidhi.nrsc.gov.in/NISAR/NISAR_Data_Product_Format_Document_V1.2.1_digisigned.pdf
- ASF download/authentication manual: https://docs.asf.alaska.edu/asf_search/downloading/
