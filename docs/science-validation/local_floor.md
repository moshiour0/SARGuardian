# Scientific validation: `src/local_floor.py`

**Code pin:** `f038e4dfdb53c8d9793aaaefe173570a1032026c`  
**Blob:** `b37816d877a37bf5eb83df4ca0738d94227d07f0`  
**Find Displacements use:** skipped without explicit target coordinates (`NOT_RUN_NO_IN_AOI_TARGET`). Keep this target/window estimator in Known Place Analysis.

## Method

Uses robust MAD scatter (1.4826 × MAD) and a three-sigma pairwise floor divided by pair span. It provides a percentile-bootstrap interval and also estimates a block-median floor to expose correlated errors between neighbouring offset estimates. It refuses zero-MAD “perfect instrument” floors and requires minimum finite pixel support.

## Scientific assessment

This is a point/window observability floor, not a per-pixel velocity map, nor a detector for unknown locations. Pixel correlation reduces effective sample size; the code itself correctly warns that the ordinary bootstrap interval is optimistic when neighbouring offset measurements share correlation windows. The block-median estimate helps expose shared error, but its tile partition and minimum six blocks are still a model choice requiring validation.

Do not use AOI centroid as a substitute target and do not turn the point floor into an AOI-wide velocity threshold. A new target-free component should estimate spatially varying uncertainty over the AOI and preserve per-pixel support.

## Tests and evidence

`tests/test_local_floor.py` tests 3-sigma/span units, quiet/noisy point versus AOI, minimum support, edge windows, bootstrap interval behaviour, AOI eligibility, provenance/CRS checks and refusal to publish zero noise.

## Required additions

- Spatial block/bootstrap calibration against offset fields with known correlation lengths.
- Sensitivity analysis for block size, window size and AOI mask geometry.
- Keep target/AOI, source hash, pair dates, sample counts and floor method attached to each result.
- Do not modify this module as part of Find Displacements without a separately approved Known Place Analysis change.

**Verdict:** PASS for its guarded point/window floor calculation; intentionally NOT RUN in target-free Find Displacements.
