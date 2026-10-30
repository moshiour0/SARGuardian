# Scientific validation: `src/pixel_stack.py`

**Code pin:** `f038e4dfdb53c8d9793aaaefe173570a1032026c`  
**Blob:** `758be2a93c6c64d0edc4c13fdb31f26af7df6e4a`  
**Find Displacements use:** conditionally executed for compatible PR GUNW exports; diagnostics run both default common-mode-removed and common-mode-retained variants.

## Algorithm

1. Aligns rasters by CRS, resolution, affine transform and common pixel lattice; refuses half-pixel shifts, rotated/sheared grids and incompatible grids.
2. Groups consecutive pairs within one geometry only; it does not bridge gaps or mix ascending and descending tracks.
3. Requires at least two consecutive intervals (three acquisition epochs).
4. Optionally subtracts each pair's spatial median. The median removed per pair is stored; a second run retains common mode for sensitivity comparison.
5. Builds cumulative relative LOS values only where the full contiguous chain is observed at that pixel. Partial and absent pair support are separately recorded.
6. Fits a linear OLS trend against actual elapsed days; outputs LOS velocity, residual-based sigma, observed-pair count and longest contiguous support.

## Scientific assessment

The grid alignment, geometry separation and explicit missing-support handling are strong safeguards. Pair-median removal can suppress genuine broad AOI-wide deformation; retaining the unremoved sensitivity product is appropriate, but the retained series can be dominated by atmosphere/orbit residuals. Neither branch is an absolute displacement reference.

The residual sigma is conditional on the linear fit and its residual degrees of freedom; it does not automatically include atmospheric, orbit, reference, unwrapping, spatial correlation or systematic errors. The code's synthetic spatial-noise null is useful as an exploratory diagnostic but uses only 20 Monte Carlo trials, so its 95th percentile is coarse and should not be treated as a calibrated false-alarm probability. The cluster cross-geometry check currently compares geographic bounding boxes of cluster extents, not exact cluster polygon overlap; a bounding-box overlap can falsely suggest same-ground corroboration.

## Tests and evidence

`tests/test_pixel_stack.py` covers half-pixel rejection, one-pair refusal, partial temporal support, support-band export and explicit common-mode removal. Backend `tools/tests/test_run_export_science_diagnostics.py` checks both common-mode variants are run and that output is not presented as hazard confirmation.

## Required additions

- Replace bounding-box cross-geometry corroboration with exact polygon/raster-mask overlap and document temporal coincidence.
- Increase and seed-control null simulations; validate block/correlation structure and multiple-testing behaviour on independent null scenes.
- Inject broad, medium and compact motion patches to measure how pair-median subtraction and local detrending change recoverability.
- Propagate per-pixel measurement uncertainty rather than relying solely on OLS residual scatter.
- Map the velocity and sigma GeoTIFFs in the dashboard as distinct layers with support masks and common-mode policy.

**Verdict:** PASS for grid/network/support invariants; PARTIAL for statistical calibration and cross-geometry corroboration; output is relative LOS velocity only.
