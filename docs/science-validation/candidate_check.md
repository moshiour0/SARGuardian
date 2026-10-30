# Scientific validation: `src/candidate_check.py`

**Code pin:** `f038e4dfdb53c8d9793aaaefe173570a1032026c`  
**Blob:** `8a1d9bed67e86737e2f72f9928a9d15d2285f06c`  
**Find Displacements use:** skipped for arbitrary WKT AOIs (`SKIPPED_UNSUPPORTED_CUSTOM_AOI`).

## Purpose

Stress-tests one already nominated candidate location by comparing target-window and AOI measurements, joining by exact granule identity, checking phase-aliasing risk and comparing geometry/seasonal hypotheses. It requires target latitude/longitude and one of the named source/langtang/lhende AOIs.

## Scientific assessment

This is not a candidate discovery algorithm: it expects a candidate point and target log as inputs. It therefore must not be repurposed by supplying an AOI centre or arbitrary first candidate. The exact-granule join is important because dates alone can collide across tracks, frames and CRIDs. Geometry/snow/atmospheric comparisons are screening evidence, not definitive physical attribution, especially without propagated covariance and independent weather/field observations.

## Tests and evidence

`tests/test_candidate_check.py` covers sensitivity sign hypotheses, blind tracks, inconclusive ratios, gentle terrain, exact granule joins, phase aliasing and non-finite geometry inputs. It remains an AOI-specific point-check, not a general WKT detector.

## Required additions

A generic AOI discovery path must produce candidate regions first and then apply a separately designed candidate-validation adapter to each region, with explicit provenance. Do not remove the named-AOI gate until polygon, candidate, DEM and reference assumptions are all parameterized and tested.

**Verdict:** PASS for its restricted target-centred use; correctly skipped in target-free custom AOI mode.
