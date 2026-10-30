# Scientific validation: `src/goff_reader.py`

**Code pin:** `f038e4dfdb53c8d9793aaaefe173570a1032026c`  
**Blob:** `fc672509f831d3982cf997aa879e909c5980dc71`  
**Find Displacements use:** executed per GOFF pair with `--layer layer2`; default reader policy also enables auto-reference and deramping.

## Measurement model

The reader loads NISAR GOFF slant-range offset, along-track offset, correlation surface peak, SNR and range-offset variance. It requires finite offset and quality arrays, valid shapes, projected coordinates for AOI clipping, correlation and SNR gates, then performs per-axis MAD outlier treatment. It preserves a pre-global-MAD `screen_valid` mask so an AOI-local coherent outlier patch remains available to the separate spatial candidate screen, while the stricter `valid` mask is used for standard robust summaries.

The default gates in this code are correlation >= 0.3, SNR >= 3, and 4-MAD culling. When enabled, a planar ramp is fitted outside the submitted AOI using pixels that are joint range/azimuth inliers. This is intended to remove orbit/coregistration residual ramps without fitting away the AOI signal.

## Validation findings

The axis-specific outlier masks and joint plane-fit sample selection are tested. The reader has strict quality-array presence/shape requirements and a pixel-attrition ledger. The workflow fixes the product to `layer2`; it does not compare the finer/multiple GOFF layers for scale consistency. It exports range and along-track offsets, correlation and robust-inlier flags; variance/SNR influence filtering but are not currently retained as full per-pixel uncertainty layers for the later AOI-wide screen.

ASF describes GOFF as geocoded pairwise pixel offsets, with spacing/window/search settings that vary by bandwidth and observation area, and explicitly warns that delivered offset layers are not post-filtered and may contain outliers. The current quality gates are therefore necessary, but they are not proof of displacement accuracy. The outside-AOI ramp-fit region is not automatically independently stable ground.

## Tests and evidence

`tests/test_goff_reader.py` covers range-only and azimuth-only outliers, upstream invalid-mask preservation, shape mismatch, active AOI mask, missing quality layers, unknown layers, and retention of coherent high-offset patches for spatial review.

## Required additions

- Export `slantRangeOffsetVariance`, `alongTrackOffsetVariance` where present, cross-covariance if available, and a provenance-linked SNR layer to the AOI-screen stage.
- Use component-specific uncertainty or covariance-aware vector scoring; do not apply a range-only detection floor to a joint range+azimuth magnitude.
- Compare layers 1 and 2 at compatible locations and report disagreement/scale-dependent support.
- Test deramp sensitivity to the size and stability of the outside-AOI fitting region.
- Calibrate correlation/SNR thresholds using held-out matched pairs and reference offsets.

**Verdict:** PASS for tested structural and reader-quality invariants; PARTIAL for measurement uncertainty and multi-resolution scientific validation.
