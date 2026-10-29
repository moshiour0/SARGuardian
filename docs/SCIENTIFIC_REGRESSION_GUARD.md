# SARGuardian Scientific Regression Guard

**Document role:** mandatory engineering and scientific acceptance contract for every AI, contributor, reviewer, and automation changing SARGuardian science or its Find Displacements pipeline.

**Scope:** first-class priority is **Find Displacements** (target-free screening over the submitted AOI and date window). **Known Place Analysis is a separate workflow** and must not be mixed into Find Displacements or treated as its prerequisite. Shared reader fixes may serve both workflows, but their input contracts, outputs, and validation must remain explicit.

**Safety rule:** a green software test, successful download, structurally valid GeoTIFF, or successful dashboard build is not proof of geophysical accuracy. Never claim confirmed movement, a validated velocity, a reproduced paper result, or a hazard alert without the evidence defined below.

## 0. Mandatory behavior for every AI change

Before editing:
1. Read this entire document and the current branch's workflow(s), pinned source commit(s), module code, tests, and recent relevant audit reports.
2. Record repository, branch, base SHA, exact source SHA used by the workflow, and the actual entry points. Do not infer that a script runs because it exists.
3. Check whether a proposed fix already exists on the current branch. Do not apply patches derived from an old pinned commit to newer code without re-auditing the current source.
4. Work only on an explicitly named feature/fix branch. Do not commit, merge, cherry-pick, or push changes to `main` unless the repository owner separately authorizes it.
5. Keep **Known Place Analysis** independent of Find Displacements. Do not require the user to nominate a target/crown polygon for target-free AOI-wide screening.

After editing:
1. Add or strengthen a regression test that fails on the original defect and passes only with the fix. Do not rely only on code inspection or a test that asserts the code repeats its prior output.
2. Run focused tests, the entire relevant test suite, static checks, and workflow/integration tests that are available. Record exact commands, commits, run IDs, results, and failures.
3. Re-open the changed files and inspect the diff for scientific, numerical, geospatial, and integration errors. Check the change against this guard again.
4. Report remaining uncertainty and unrun tests plainly. Do not label a partial or simulated validation as an independent real-world validation.
5. Never silently change scientific thresholds to increase pixel counts. A more permissive mask must be justified by signal-recovery **and** false-positive/null tests.

## 1. Required output/status vocabulary

Every stage must return a machine-readable, explicit status and a reason. Use semantically distinct states such as:

- `PASS_STRUCTURAL_ONLY`: files, shapes, metadata, provenance, and declared structural checks passed; no claim of physical validity.
- `MEASUREMENT_AVAILABLE_DIAGNOSTIC`: a measurement exists but its reference, uncertainty, or validation is not sufficient for a verified scientific claim.
- `CANDIDATE_FOR_REVIEW`: an anomaly met stated screening rules; not a confirmed landslide or hazard.
- `INSUFFICIENT_DATA`: missing coverage, too few valid samples/epochs, weak quality, or unsupported product/schema.
- `INCOMPATIBLE_GRID`, `INCOMPATIBLE_SERIES`, `REFERENCE_UNVALIDATED`, `NETWORK_DISCONNECTED`, `RANK_DEFICIENT`, `UNCERTAINTY_UNAVAILABLE`.
- `WITHHELD`: do not publish the requested metric because its prerequisites did not pass.
- `FAIL`: a required structural, provenance, invariant, or runtime check failed.

A missing measurement is **never** a numeric zero. A zero candidate count is not evidence that no displacement or hazard exists. No downstream stage may convert `INSUFFICIENT_DATA`, NaN, missing values, or exceptions into a zero measurement or an unqualified `PASS`.

## 2. Pixel attrition audit: mandatory for every reader and detector

For every product, layer, and processing variant, persist a reproducible pixel-gate ledger. At minimum report:
- native input dimensions, dtype, nodata/fill definition, CRS, transform/posting, pixel centres, product/version/CRID, processing variant, and source SHA-256;
- total source cells;
- finite / source-valid cells;
- count after **each named gate**, in exact execution order;
- count rejected by each gate and why (mutually exclusive first-rejection counts where possible; otherwise label overlapping reasons);
- cells intersecting/inside AOI under the explicitly declared pixel-centre or all-touched rule;
- final reader-valid pixels, exported valid pixels, and candidate pixels;
- percentages using a relevant denominator (AOI-covered cells for AOI coverage, not the entire product frame);
- per-layer and spatial maps of gate rejection, including the target/candidate/feature locations when relevant;
- gate configuration, software commit, and product identity in the output provenance.

**Non-negotiable mask invariant:** never overwrite the only copy of an earlier validity mask. Keep named masks such as `source_valid`, `quality_valid`, `aoi_valid`, `analysis_valid`, and `candidate_mask` separately where applicable. A later scientific stage must be able to distinguish source-invalid, quality-rejected, AOI-excluded, statistically filtered, and merely not-selected pixels.

Every filter must have a test proving (a) invalid pixels stay invalid, (b) valid pixels are not silently lost for an unrelated reason, (c) rejection counts reconcile, and (d) a physical/synthetic anomaly of known truth can pass the intended path.

## 3. GUNW (unwrapped interferometric phase) safeguards

### 3.1 Product schema and masks
- Decode the mask according to the actual product release/CRID, dtype, dataset metadata/attributes, and official schema. Do not infer a packed mask's meaning from its decimal representation or assume every release has the same layout.
- For a supported schema, required quality datasets must be present and be compatible with the phase grid. Missing mandatory mask/quality data, unknown dtype, conflicting shapes, or ambiguous dataset paths must be explicit errors or a separately justified degraded mode. Never silently skip a required quality gate.
- Report each of: finite/nodata, coherence, connected-component, water/subswath, anomaly, ionosphere-screen interpolation, reserved/unknown bits, and AOI masks separately.
- Respect the product's documented fill semantics. Do not classify every numerical phase value of zero as invalid unless the particular product schema establishes zero-fill and the other mask evidence is consistent. Include a test for a valid zero-phase / zero-displacement observation and for a true fill pixel.
- Confirm phase, coherence, connected-component, mask, ionosphere screen, coordinates, and projection are on the same declared grid. Fail closed on an unexplained mismatch.

### 3.2 Ionosphere screen
- Log whether the correction is present, shape-compatible, applied, and subtracted; report its units and the declared sign convention.
- Treat a pixel flagged as interpolated in the correction screen as a **quality/uncertainty issue**, not as automatically proven physical garbage or automatically safe data.
- Preserve a conservative production policy until independent evidence justifies changing it. Run correction-on and correction-off sensitivity tests, but label correction-off output as a diagnostic; never silently substitute it for the production product.
- Before changing an interpolation gate, quantify attrition, residual/uncertainty, detection sensitivity, false positives, and performance on independent controls/known-truth simulations. Never weaken the gate solely to increase retained pixels.

### 3.3 Reference and connected components
- GUNW phase is relative. A reference sample must be valid and belong to the same unwrapped connected component as the pixel/region whose value is reported.
- Do not use a single scalar reference to combine disconnected phase components. Keep components separate, or explicitly mask cross-component values from comparison/inversion.
- The selected AOI's dominant component is not automatically the component of every candidate. Candidate-specific/component-specific support must be evaluated.
- Record the chosen reference location/component, reference sample count, reference statistic, in-AOI/out-of-AOI status, and stability evidence.
- If a stable reference cannot be established, emit relative/reference-unvalidated diagnostics or withhold a stronger claim; do not manufacture an absolute displacement.

## 4. GOFF (pixel offsets) safeguards

- Verify the correct data layer, offset arrays, correlation, SNR, coordinates, CRS, shape, and product identity. Missing required quality layers or mismatched shapes must not silently disable gates.
- Keep range/slant-range and along-track offsets as distinct measured components with explicit units, sign, geometry, and independent validity/uncertainty. They are not by themselves a complete 3-D ground displacement vector and must not be mixed with GUNW phase values as though they were identical observations.
- Apply AOI clipping on the documented geospatial grid. Keep AOI denominator and pixel-centre/all-touched semantics explicit.
- If using an auto-reference, plane deramp, or background fit, record the fit area, number/distribution of inliers, residuals, condition/rank, and whether the fit has adequate support outside the signal AOI. Do not call the outside-AOI fit area “stable ground” without evidence of stability.
- Do not allow a deramp failure or insufficient background support to masquerade as a normal measurement. Mark output quality as degraded or withhold dependent measures when residual ramps could dominate.
- MAD/outlier culling is not proof that rejected pixels are noise. Preserve the upstream quality-valid offset field and record culling masks independently. Validate thresholds on (i) isolated gross noise outliers, (ii) compact coherent displacement patches, (iii) narrow/elongated movement features, and (iv) background trends. A spatially coherent, quality-supported anomaly must remain reviewable even if its magnitude is a global statistical outlier.
- Do not require both range and azimuth components to be finite or inlier if that requirement would make a valid one-axis measurement unusable without an explicit scientific reason. When vectors are analysed, declare whether both axes are required and why.
- Do not present a GOFF pair's offset divided by elapsed days as a validated ground velocity without a validated common datum, geometry, uncertainty, and temporal model.

## 5. AOI and grid invariants

- Parse only the declared supported polygon contract; validate WGS84 coordinate order, finiteness/ranges, ring closure, polygon validity, and unsupported holes/multipolygons explicitly.
- Never silently substitute a named/historical AOI, target, or fallback polygon for the submitted AOI. Never replace exact author granule/product identity with another product merely because it is available.
- Confirm CRS transformations use longitude/latitude axis order explicitly and fail closed if required CRS metadata is missing.
- Alignment must be based on geospatial transforms and grids, not array index or shape alone. Record CRS, origin, pixel size, shape, transform, and grid policy.
- Do not silently resample between incompatible grids. If resampling is introduced, state interpolation, its effect on uncertainty, and an independent test; otherwise mark the comparison incompatible.
- AOI boundary decisions must be explicit. Test pixels fully inside/outside and on/near polygon boundaries. Never report an outside-AOI target or pixel as an in-AOI measurement.
- A valid grid cell, finite measurement, or coverage fraction is not evidence of deformation accuracy.

## 6. Temporal network and velocity invariants

- Parse only exact product identities and acquisition dates. Preserve processing variant, CRID/maturity, track/direction, frame/subswath, polarization, frequency/bandwidth, and spatial footprint. Do not merge distinct footprints/releases/processing variants merely because date pairs match.
- Use **measured pair edges**, not catalog-only edges, to determine numeric network topology. Separate graph-topology diagnostics from numerical inversion.
- Before inversion, require a connected epoch component, compatible measurement type/geometry/grid, a common validated datum (or explicitly named diagnostic datum), enough independent observations, full column rank, adequate temporal span, and appropriate residual/error degrees of freedom.
- A graph can be connected yet exactly determined with zero redundancy. Such a graph cannot estimate residual variance from its own closure. Identify loops/redundancy and do not describe an exactly determined trend as independently tested.
- Never invoke least-squares minimum-norm output from a rank-deficient graph as a valid cumulative series. Check rank/condition before publishing any values.
- Never sum pair medians from separately re-referenced GUNW/GOFF observations as though they share a zero. If common datum is AOI-derived rather than independently verified stable ground, label every output diagnostic/reference-unvalidated and test sensitivity to the choice of datum.
- Account for irregular observation intervals using actual acquisition dates, not pair order or assumed fixed revisit cadence. Report missing epochs and disconnected time blocks separately.
- Never replace missing per-pixel measurements with zero. If complete-case stacking is used, quantify and map the support loss and test for bias when missingness correlates with deformation/decorrelation. Prefer an explicitly justified estimator and minimum-support policy.
- Do not combine GUNW and GOFF or ascending and descending observations in one inversion unless the measurement model explicitly harmonizes units, sign, datum, geometry, acquisition support, and covariance. A 2-D decomposition with a north-motion assumption is not a full 3-D solution.
- Report uncertainty, covariance assumptions, residuals, fit degrees of freedom, network redundancy, leave-one-pair-out sensitivity, and the limits of interpretation. If these cannot be estimated, say so.

## 7. Pixel stack and anomaly screen safeguards

- A two-epoch/one-pair stack cannot estimate residual-based slope uncertainty from its own fit. Do not substitute an artificial `max(dof, 1)` for a genuine degree of freedom and then publish the result as a statistical uncertainty. Withhold unsupported significance claims.
- Test common-mode removal both enabled and disabled. It can suppress broad or AOI-wide real motion as well as shared atmospheric/reference error. Output both as sensitivity diagnostics when useful; do not choose whichever supports a preferred conclusion.
- Spatial pixels and neighbouring GOFF estimates are correlated. Do not use independent-pixel false-positive rates without a spatial-null model.
- Use positive and null controls, held-out spatial/temporal data, multiple feature sizes and amplitudes. Record threshold sensitivity and calibration scope.
- If temporal persistence is required, also preserve a separate one-pair/low-persistence review layer. A catastrophic event may appear in one acquisition pair only or may destroy coherence. A persistence gate must not erase the evidence without a trace.
- A zero-candidate output is “no candidate met these filters on the usable data,” not “no movement”.
- A candidate region is a screening result, not a confirmed event, landslide, glacier collapse, forecast, or public warning.

## 8. Local floor, inverse velocity, geometry, and bound

- A local detection/noise floor is not measured motion, a stable reference, a velocity, or a forecast. Require enough spatial support and report how spatial autocorrelation limits effective sample size.
- Never emit a numeric noise floor from zero/nonfinite MAD; zero MAD does not prove zero noise.
- Inverse-velocity/Fukuzono analysis is an optional downstream analysis, not a general movement detector. Require a defensible time series, consistent direction/sign interpretation, matched per-pair floors, enough above-floor observations, an acceleration fit, a bounded/appropriately reported confidence interval, and strict exclusion of post-event observations for retrospective forecasting tests.
- Geometry inversion must use product look vectors when available and report matrix rank/condition number and propagated uncertainties. Mark weak-sensitivity tracks unresolved. State any imposed zero component; never call an east/vertical two-component solution full 3-D motion.
- Candidate/event-specific bound calculations must not run with generic AOI priors, stale candidates, an out-of-AOI target, an invalid floor, or incompatible LOS/down-slope units. A bound is conditional detectability, not proof that a precursor did not occur.

## 9. Synthetic and real-data test matrix

Every material science change must add or preserve test coverage for applicable cases below.

### GUNW
- Known positive, negative, and zero displacement sign/scale/wavelength round-trips.
- Valid zero-phase pixel versus true fill/nodata.
- Every supported mask dtype/schema/release; water, missing subswath, anomaly, reserved/unknown bits, and ionosphere interpolation bit.
- Missing/shape-mismatched mask, coherence, component, coordinates and ionosphere layers fail closed or enter an explicitly justified degraded status.
- Reference in the same component; disconnected components with whole-fringe offsets; an AOI that contains multiple components.
- Sparse coherence, partial validity, AOI edge pixels and no stable outside-AOI reference.
- Correction-on/off comparison with explicit diagnostic-only status for the non-production branch.

### GOFF
- Known range-only and azimuth-only signals, offset sign/units, correlation/SNR invalidity, absent required layers and mismatched shapes.
- Isolated gross outlier and compact coherent high-amplitude movement patch survive as separate cases; assert which gates reject each and why.
- Deramp inliers are spatially distributed and from the declared fit region; ill-conditioned/too-small background fails explicitly.
- Culling threshold/local-window/direction-coherence/persistence sensitivity and AOI-boundary tests.
- Candidate map maintains a traceable path from source pair and pre-filter quality-valid pixels to every final candidate.

### Network/stack/integration
- Chain, loop/redundant, disconnected and rank-deficient networks; missing measurements; duplicate date pairs across tracks/CRIDs; irregular dates; incompatible grids/datum.
- Exactly-determined network has no false formal residual uncertainty; leave-one-out breaks are not bridged with arbitrary zeros.
- Two-epoch stack cannot report a significance test; one missing epoch cannot silently become zero.
- Broad AOI-wide motion is not hidden without an explicit common-mode warning.
- Candidate exists in one pair but not two; a zero-candidate result is explicitly scoped to filters/support.
- The workflow actually invokes the intended module; output summary, raster/GeoJSON/CSV, provenance and map assets are produced and validated. A green reader-only or dashboard-only job is not an E2E science pass.

## 10. Workflow and UI acceptance contract

Before describing Find Displacements as integrated:
1. Trace the UI/API dispatch fields to the exact GitHub workflow, branch/ref, inputs, product discovery, per-product reader, artifact aggregator, candidate screen, and downloadable/map outputs.
2. Verify by logs and source SHA that the run used the intended workflow and immutable source pin.
3. The workflow must run target-free screening without inventing a target coordinate or requiring a user-selected crown/landslide polygon.
4. The map must expose AOI-wide review candidates and their provenance/status, not merely the submitted AOI outline or selected target marker.
5. Velocity is shown only if a valid same-family/geometry network, common datum, rank/redundancy, temporal support, and uncertainty gates pass. Otherwise show why it is withheld.
6. Candidate screen outputs and velocity products must be linked into the **same expected aggregator build and UI contract**. Code on an experimental branch is not integrated merely because a test can run against a retained artifact.
7. Keep Known Place Analysis as a separate deeper-analysis workflow. Do not let its target-local floor/measurement prerequisites block the target-free Find Displacements feature.
8. Keep raw-data retention, checksums/provenance, Google Drive upload/receipt contract, and explicit raw-replay authorization gates intact. Do not launch an expensive raw replay without the required explicit authorization.

## 11. Required change record for every fix

Every PR/commit that changes science should provide:
- observed defect and the exact code path/line(s);
- intended scientific invariant and product/schema assumptions;
- affected files/workflows and exact pinned SHA changes;
- before/after pixel attrition table and representative spatial masks when the change affects filtering;
- positive truth test and null/negative control;
- focused test command/results and full-suite results;
- live workflow run ID/URL only when it was actually run;
- effect on known historical/target test outputs and whether that comparison is scientifically comparable;
- remaining uncertainty and restrictions;
- explicit statement that main has not been modified and no merge occurred, if applicable.

## 12. Final release blockers (must all be cleared or explicitly unresolved)

- [ ] All pixel-gate counts reconcile and the rejection maps/ledger are present.
- [ ] Every mask/schema and grid is verified; missing required quality input cannot silently bypass a gate.
- [ ] Disconnected GUNW components cannot share an invalid scalar reference.
- [ ] Coherent high-amplitude GOFF anomalies remain reviewable through the outlier policy.
- [ ] Numeric inversion refuses rank-deficient, disconnected, or incompatible measurements.
- [ ] Unsupported degrees of freedom never produce misleading uncertainty or significance.
- [ ] Target-free AOI screening runs in the same branch/ref that the expected aggregator/UI actually consumes.
- [ ] Candidate outputs are visible, georeferenced, provenance-bound, and explicitly non-alert unless separately validated.
- [ ] Unit, regression, synthetic truth, null-control, and integration tests are green; unresolved items are labelled as such.
- [ ] No main-branch commit or merge has occurred without separate explicit authorization.

**Scientific principle:** maximize recoverable valid evidence, not the raw count of surviving pixels. Preserve alternative masks and diagnostics; reject only for a documented, testable reason; never turn missing, uncertain, or unvalidated measurements into certainty.
