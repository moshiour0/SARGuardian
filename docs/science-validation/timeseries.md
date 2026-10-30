# Scientific validation: `src/timeseries.py`

**Code pin:** `f038e4dfdb53c8d9793aaaefe173570a1032026c`  
**Blob:** `4f2c83af18a79110d255511e9609e39dbd6fd60c`  
**Find Displacements use:** the time-series module supplies network identity/topology/inversion primitives to the target-gated common-reference diagnostic; the numerical target-free common-datum inversion is currently withheld.

## Algorithm and contracts

- Builds pair edges from actual NISAR metadata and acquisition dates.
- Keeps product family, processing type, track, frame, subswath, bandwidth, polarization and CRID in network identity; conflicting processing releases are not silently guessed.
- Computes connected temporal components. Separate components keep separate zero references and cannot be merged.
- Inverts pair differences using least squares with the first epoch fixed to zero.
- Rejects rank-deficient/ill-conditioned systems; reports redundancy and degrees of freedom.
- Leave-one-out jackknife refuses pair removal that disconnects the network.

## Scientific assessment

This follows the essential SBAS graph principle: use a compatible network of interferograms, do not manufacture links between disconnected temporal components, and distinguish exactly determined solutions from redundant solutions with residual degrees of freedom. Berardino et al. (2002) describe small-baseline networks linked via SVD; atmospheric phase artefacts are an important part of time-series interpretation.

A pair graph inversion cannot cure inconsistent per-pair phase datums. All edges must represent comparable measurements on a common datum or explicitly modeled pair offsets. Its reported residual errors also do not automatically contain atmospheric/reference systematic uncertainty.

In current Find Displacements, `dashboard/common_reference.py` exits before importing the science solver unless a validated target is present. Therefore no claim should be made that this complete solver currently produces a target-free AOI velocity map.

## Tests and evidence

`tests/test_timeseries.py` covers known cumulative-series recovery, disconnected epochs, redundancy, loops/residuals, duplicate-product identity, processing/CRID ambiguity, rank deficiency and exact-determined zero-DOF behaviour.

## Required additions

- Build a target-free reference-selection adapter independent of the Known Place target contract.
- Require a redundant, compatible pair graph and stable-reference sensitivity evidence; preserve disconnected components separately.
- Propagate measurement covariance and atmospheric/reference uncertainty into the inversion.
- Add synthetic network tests with missing observations, cycle closure errors, datum offsets and spatially varying masks.

**Verdict:** PASS for graph construction and least-squares safeguards; NOT YET AVAILABLE as a target-free, independently referenced AOI-wide time series.
