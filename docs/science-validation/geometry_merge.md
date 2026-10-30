# Scientific validation: `src/geometry_merge.py`

**Code pin:** `f038e4dfdb53c8d9793aaaefe173570a1032026c`  
**Blob:** `a95bdb3b5e3522a0c61a237ba12cfc792522b591`  
**Find Displacements use:** withheld (`BLOCKED_NO_COMPARABLE_ASC_DESC_SERIES`).

## Method

Computes track/terrain sensitivity and, when compatible per-geometry LOS series exist, can resolve constrained downslope or east/up components. It records sensor look side explicitly and prefers per-product LOS vectors when available. The code refuses poorly conditioned geometry and low-slope terrain where aspect uncertainty dominates.

## Scientific assessment

One interferogram measures only the displacement component along its LOS. Ascending and descending LOS measurements can constrain selected components only under assumptions and sufficient independent geometry; two near-polar looks do not generally determine unconstrained east, north and up displacement. The code's east/up reduction relies on assumptions that must be explicit. A generic user AOI cannot inherit the fixed named-track geometry/DEM context without checking the actual product vectors and terrain.

Tests exercise mission look-side sign, geometry sensitivity, blind tracks, slope floor and ill-conditioned east/up decomposition. These prove the implementation's selected geometry cases, not that a particular live NISAR target has reliable 3-D motion.

## Required additions

- Require actual per-pixel look vectors and incidence layers from each product where available.
- Match spatial support, acquisition intervals, CRS, reference and compatible measurement semantics before decomposition.
- Propagate LOS covariance, DEM slope/aspect error and ill-conditioning into component uncertainty.
- Publish the condition number and withheld status; never label a two-look product as full 3-D without a third constraint or explicit model assumption.

**Verdict:** PASS for tested geometry mathematics; not currently eligible for Find Displacements 3-D velocity and correctly withheld.
