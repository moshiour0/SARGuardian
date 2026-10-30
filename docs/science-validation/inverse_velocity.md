# Scientific validation: `src/inverse_velocity.py`

**Code pin:** `f038e4dfdb53c8d9793aaaefe173570a1032026c`  
**Blob:** `02d4fc49833af09057f87734ff94e444d7aad67f`  
**Find Displacements use:** withheld (`BLOCKED_NO_VERIFIED_COMMON_DATUM_SERIES`).

## Method

Implements the Fukuzono/Voight inverse-velocity construction for accelerating creep: under the specific `alpha = 2` form of `dv/dt = A v^alpha`, inverse speed is linear in time and extrapolates to a finite failure-time intercept. The code groups by geometry/network component, supports irregular intervals, applies measured pair floors, refuses zero speed/direction reversal, handles event-date cutoffs and withholds forecasts when uncertainty is unbounded.

## Scientific assessment

The model is conditional, not universal. Fukuzono (1985) describes controlled slope-model experiments; the inverse-velocity straight-line relation follows specifically from the acceleration-law assumptions. A successful fit is not proof of a physical failure mechanism or a reliable forecast. It requires a valid measured displacement/velocity sequence, consistent geometry and datum, adequate sampling, a compatible noise floor and independent event/cutoff semantics. It is unsuitable as a target-free AOI scan across millions of pixels.

## Tests and evidence

`tests/test_inverse_velocity.py` covers textbook synthetic failure, sign invariance, direction reversal, pair-specific floors, post-event exclusion, irregular intervals, small-sample t intervals, zero-speed refusal and unbounded failure-time uncertainty.

## Required additions

- Do not invoke until a candidate-specific or pixel-cluster series has passed reference, temporal-network and uncertainty gates.
- Evaluate alpha/model misspecification, non-accelerating and decelerating controls, autocorrelated noise and false alarms on held-out scenes.
- Any forecast must report fitted interval, sample count, noise-floor provenance, residuals, uncertainty interval and the model assumptions.

**Verdict:** PASS for tested synthetic/invariant safeguards; not an AOI-wide discovery method and correctly withheld in current Find Displacements.
