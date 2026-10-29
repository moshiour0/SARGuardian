# Custom-AOI GOFF spatial screen — research prototype

## Status and intended use

This module is a **target-free exploratory ranking tool**, not an operational landslide/glacier-collapse detector. It searches for repeatable, spatially coherent anomalous offset fields inside an arbitrary single-ring WGS84 polygon. A candidate means “review these pixels against geology, terrain, optics, and independent controls”; it does not mean “a landslide was detected”.

The module does not require `target_latitude` / `target_longitude`. That is deliberate: the target is not known at discovery time. This is separate from `local_floor.py` and target-local measurements, which do require an explicit measurement point.

## Why it reads raw GOFF HDF5

A clipped 3-band display GeoTIFF is not a sufficient source for a defensible offset screen. The raw GOFF source contains per-pixel slant-range and along-track offsets, SNR, correlation, coordinates and the product's radar-grid geometry. The screen requires those fields to apply quality gates and convert the two offset components to horizontal East/North at the local imaging geometry.

The implementation processes **one HDF5 pair at a time** and accumulates only the spatial candidate counts and score rasters. It intentionally does not keep all pair-sized arrays in RAM. Every input pair is recorded with source SHA-256 and byte count.

## Input contract

- Full, unmodified NISAR L2 GOFF HDF5 product filenames so the identity can be checked.
- One single-ring `POLYGON` WKT in WGS84 longitude/latitude order.
- A georeferenced DEM that covers the AOI. Its vertical datum must be verified against `radarGrid/heightAboveEllipsoid`; missing elevation is not silently imputed.
- At least two compatible acquisition pairs for the persistence screen.
- The selected polarization, frequency and offset layer must exist in every product.

The screen **fails closed** if track, direction, frame, bandwidth, polarization, processing type or CRID differs; if duplicate calendar-day acquisition intervals are included; if a pair is reversed; if a requested layer is missing; or if the native pixel grids do not match exactly. It does not silently reproject one pair to make it fit another.

It is acceptable for the subswath ID to change when the source products still map to the exact same exported grid and series identity. The duplicate-interval gate prevents two subswaths of the same acquisition pair from being counted as independent temporal evidence.

## Run locally after the environment and DEM datum have been checked

```bash
python src/goff_spatial_screen.py \
  --read /path/to/NISAR_L2_PR_GOFF_...P05023_N_F_J_001.h5 \
         /path/to/next-compatible-GOFF-pair.h5 \
  --dem /path/to/verified-dem.tif \
  --dem-vertical-crs EPSG:3855 \
  --correlation-min none \
  --aoi-wkt 'POLYGON ((lon lat, lon lat, lon lat, lon lat))' \
  --out outputs/custom-aoi-goff-screen
```

This command does not download data and is not wired into the production backend workflow. The caller is responsible for selecting the correct input pairs and not mixing processing releases or acquisition geometries.

## What the output contains

- `pair_persistence_count.tif`: number of supplied pairs in which a pixel passes the screen.
- `mean_anomaly_m.tif`: mean positive local-magnitude anomaly over pairs that flagged the pixel.
- `mean_directional_coherence.tif`: mean coherence of local residual vectors in flagged pairs.
- `persistent_candidate_mask.tif`: connected components that persist in at least the configured number of pairs and meet the minimum component-size gate.
- `candidates.geojson` (WGS84 coordinates) and `candidates.csv`: review polygons and coordinates.
- `screen-summary.json`: AOI, grid, pair identity, source hashes, parameters, candidate counts, output hashes, and limitations.

Empty output means no candidate met these exploratory gates on the usable pixels; it does **not** mean no motion or no hazard exists. A product with no usable AOI pixels is a data-support limitation, not a measured zero displacement.

## Algorithm and provenance

For each pair, the screen:

1. Applies finite-value, SNR, correlation, DEM and imaging-geometry gates.
2. Converts range/along-track offsets to horizontal East/North using per-pixel `radarGrid` geometry interpolated at DEM elevations; it rejects geometry extrapolation beyond the metadata support.
3. Removes the local vector mean and computes a local magnitude anomaly relative to robust local scatter.
4. Requires a sufficiently large anomaly and high directional coherence of the residual vectors.
5. Counts spatially coincident hits across independent pairs and retains only persistent connected components.

The vertical conversion uses the public [PROJ EGM2008 geoid grid](https://cdn.proj.org/us_nga_egm08_25.tif). The Copernicus DEM height reference is EGM2008; see the [Copernicus DEM documentation](https://documentation.dataspace.copernicus.eu/APIs/SentinelHub/Data/DEM.html). Conversion is applied only when the caller explicitly declares `EPSG:3855`; otherwise the caller must declare that heights are already ellipsoidal.

The initial thresholds are adapted from Meena (2026), *Bhote Koshi 2026 ice–rock avalanche: NISAR L2 GOFF site series and footprint-wide offset-tracking screen — data and code*, Zenodo DOI [10.5281/zenodo.22730419](https://doi.org/10.5281/zenodo.22730419). Use `--correlation-min none` for an author-like SNR-only quality gate, or leave the default `0.30` to add SARGuardian's stricter correlation filter. The latter is deliberately not a paper-equivalent gate. Component extraction uses 4-connectivity to match the author's default connected-component labelling. Even with `none`, this remains a reproduction-oriented implementation, not a validated exact paper reproduction: the geometry interpolation, vertical-datum correction, AOI edge support and dependency versions still must be cross-checked against the author's source implementation. Those thresholds are not calibrated for arbitrary AOIs. The output includes this limitation explicitly.

## External benchmark: what not to claim

The supplied Zenodo archive contains five site-area derived offset-field GeoTIFFs, author-produced footprint-wide candidate rasters/inventory, terrain layers and tables. It does not contain the full raw GOFF HDF5 products needed to run this module end to end. These derived fields are useful to check published summary values and map plumbing, but they are not interchangeable with raw GOFF because the full per-pixel SNR and original radar-grid geometry are not retained.

The author's five site intervals report crown median displacement magnitudes of approximately 0.99, 0.54, 0.30, 0.96 and 0.33 m, over 12, 12, 12, 24 and 12 days respectively. The accompanying size-matched polygon null shows that only two of the four pre-event intervals were resolved by one estimator, and neither of those crown coherence statistics crossed the 95th percentile of the reported spatial null (77.2% and 85.2%). Therefore reproducing a large median is not sufficient to claim a statistically credible precursor. Compare the exact zone, interval, estimator, null distribution and processing release.

This prototype has **not** rerun the source HDF5 and has not independently reproduced those numbers. It has only passed synthetic controls and package regression tests.

## Required validation before any operational integration

- Verify the DEM vertical datum against the product's ellipsoidal radar-grid heights, and document the vertical-datum conversion.
- Run first with `--correlation-min none` to align the author's SNR-only screen, then repeat with `--correlation-min 0.30`; compare hit masks, valid-pixel denominators and components.
- Use matched spatial nulls and multiple-testing correction for the actual AOI, not just a z-score or a “persistent” label.
- Test against held-out known positive and negative controls; quantify false alarms and missed detections.
- Independently validate candidate locations, terrain orientation, glacier/debris-cover context and source-product maturity.
- Test threshold sensitivity, missing pairs, different CRIDs, no-data and mismatched grids.
- Keep GUNW phase-derived LOS measurements and GOFF offset measurements as separate families. Do not derive cumulative displacement or inverse velocity from independently referenced pair medians.

No hazard alert, failure probability, precursor claim or operational warning should be produced by this prototype.
