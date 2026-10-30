# Scientific validation: `src/gunw_reader.py`

**Code pin:** `f038e4dfdb53c8d9793aaaefe173570a1032026c`  
**Blob:** `444386c2004dba07ba7507b8782b39c4c725f600`  
**Find Displacements use:** executed per GUNW pair using `--read ... --aoi custom --auto-ref --csv ... --export ...`.

## Measurement model

The reader converts unwrapped interferometric phase `phi` to LOS displacement using `d_LOS = -lambda * phi / (4*pi)`, taking the product wavelength/frequency when present. Positive is defined as away from the satellite. It uses the unwrapped-phase group's coherence, connected-components and packed mask datasets, AOI grid, and ionosphere phase screen. With default policy, the ionosphere screen is subtracted and interpolated-screen pixels are conservatively excluded. The reference is pair-specific and selected from a fully valid high-coherence block outside the AOI when possible; it may fall back inside the dominant connected component and records that the datum is unvalidated.

## Validation findings

- Dataset discovery is schema-aware and refuses to silently skip required mask/connected-component/ionosphere layers.
- P05023+ packed uint32 mask is decoded by bit fields rather than interpreted as decimal digits; legacy uint8 schema is separately supported.
- Reference values are not allowed to leak across unwrapping connected components.
- A zero phase can be a valid measurement; it is not discarded merely because zero is also a fill convention.
- AOI gate attrition is reconciled and reported.

The formula is physically appropriate for interferometric LOS range change. It does not yield downslope or 3-D motion. A per-pair auto-reference does not guarantee that separate interferograms share a common stable datum. The reader's ionosphere subtraction is a deliberate project policy; ASF says the distributed GUNW ionosphere screen is included but is not applied to the product data layers by default. Its interpolation/uncertainty treatment therefore needs a sensitivity audit against the screen uncertainty layer. ASF also states available external tropospheric corrections are not removed from the interferometric data by default.

## Tests and evidence

`tests/test_gunw_reader.py` covers known-displacement round trips, product wavelength, coherence-grid identity, ionosphere subtraction, mask semantics, connected-component-safe reference, fail-closed missing datasets, zero phase, and AOI gate-ledger reconciliation.

## Required additions

1. For every supported CRID/schema, validate the packed-mask bit definitions against product documentation and raw HDF5 examples.
2. Preserve both ionosphere-corrected and uncorrected sensitivity outputs with screen uncertainty; never silently switch policy.
3. Quantify pair-reference drift over the AOI and build a target-free stable-reference method before temporal AOI-wide LOS inversion.
4. Add atmosphere/topography sensitivity or external weather/reanalysis support; publish residual uncertainty rather than imply it is included.
5. Validate sign and scale with an independent known-displacement source, not only a synthetic array round trip.

**Verdict:** PASS for reader arithmetic/schema guards under tested schemas; NOT YET VALIDATED as an AOI-wide common-datum deformation series.
