# Find Displacements — target-free AOI discovery boundary

## Product boundary

SARGuardian has two different user workflows. They must not be merged into one pipeline.

### A. Find Displacements (this prototype's intended home)

Purpose: a user submits an arbitrary AOI polygon and date window **without choosing a target point or crown polygon**. The system screens the whole usable AOI and returns spatially explicit candidate patches with evidence and uncertainty.

Expected inputs:
- AOI WKT/polygon and time window.
- NISAR GOFF pairs and, as a separate corroboration layer, available GUNW pairs.
- Explicit data/quality/geometry configuration and product provenance.

Expected outputs:
- Per-pair usable-pixel support and exclusion reasons.
- AOI-wide candidate patches, polygons, and per-patch displacement/anomaly statistics.
- Persistence across independent acquisition intervals, matched control/null statistics, and candidate false-alarm context.
- Separate GOFF offset evidence and GUNW LOS-phase evidence. A comparison view may show both, but the measurements must not be merged as if they were the same quantity.
- No user-specified target is required for discovery.

The current `src/goff_spatial_screen.py` is only an **exploratory GOFF prototype** for this workflow. It does not yet constitute a scientifically validated AOI-wide detector, does not yet have calibrated spatial-null/multiple-testing control, and is not wired into production actions or the dashboard. It must not be described as a validated detection engine.

### B. Know Place Analysis (separate workflow)

Purpose: reproduce and audit a known, externally defined scientific case: the exact paper/site/crown polygon, author-selected interval, processing release, estimator and control/null distribution.

This workflow may use a target crown polygon because the research question is defined around that known feature. It must preserve the author's exact product identity (including P05023 versus X05026 distinctions), interval and licensing boundaries.

Know Place Analysis is **not** the implementation target of the target-free prototype. Do not import its target polygons, event-specific priors, benchmark thresholds or paper-specific assumptions into Find Displacements. Do not use the target-free screen as a substitute for exact paper reproduction.

## Current implementation status

- Synthetic tests check a coherent patch, noise-only data, AOI masking, product identity and explicit DEM datum declaration.
- The raw-HDF5 path uses per-product GOFF offsets, SNR/correlation and radar-grid geometry.
- Pair processing is streamed one HDF5 product at a time.
- Terrain direction is reported as post-hoc review metadata, not used as proof or as a hazard alert.
- A result is exploratory ranking only; an empty result is not evidence of no movement.
- No raw replay/download has been started by this code branch.
- Main has not been modified by this prototype.

## Required gates before any production integration

1. Verify range/along-track to East/North sign, axes, units and geometry interpolation against official NISAR examples and known synthetic vectors.
2. Verify DEM vertical conversion with a controlled test and log the actual PROJ operation/grid used.
3. Add matched spatial-null distributions and multiple-testing correction for the AOI.
4. Validate on independent positive and negative AOIs, reporting false positives, missed detections and sensitivity to thresholds.
5. Check temporal independence: do not count same-calendar-interval subswaths as separate evidence.
6. Require enough compatible pairs and meaningful common valid support before persistence claims.
7. Add GUNW corroboration only as a separate evidence family with phase/coherence/connected-component/ionosphere sensitivity checks.
8. Add dashboard wiring only after the above scientific gates pass; keep family-specific playback and measurement units explicit.

Until then, do not label candidates as landslide/glacier collapse detections, precursors, forecasts or warnings.
