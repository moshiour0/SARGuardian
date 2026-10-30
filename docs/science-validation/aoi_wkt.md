# Scientific validation: `src/aoi_wkt.py`

**Code pin:** `f038e4dfdb53c8d9793aaaefe173570a1032026c`  
**Blob:** `36ac1d2a50fd3dcf642cf7d1709f34962ac98136`  
**Find Displacements use:** executed at AOI discovery/validation, before science workers.

## Function and contract

Parses a user-supplied WGS84 polygon WKT. It accepts a single exterior ring, longitude-latitude coordinate order, optional SRID 4326, up to 500 vertices; it closes an open ring and rejects invalid coordinate bounds, non-finite values, repeated consecutive points, self-intersection and zero area. Holes and MULTIPOLYGON are deliberately rejected instead of silently flattened.

## Scientific/geospatial assessment

This is a safety-critical input contract, not a displacement estimator. Strictly rejecting unsupported geometry is preferable to silently changing the user's AOI. The polygon hash must be carried unchanged through discovery, reader clipping, candidate screening and dashboard provenance. A single-ring restriction is an explicit product limitation; do not claim holes or multipart polygons are supported.

The parser's zero-area check is performed in longitude/latitude coordinates. For very large or dateline-crossing polygons, planar lon/lat area and simple ring logic are not a substitute for geodesic topology validation. Such AOIs should be rejected or routed through a tested geospatial library path, not approximated.

## Tests and evidence

`tests/test_aoi_wkt.py`: open/closed ring normalization, SRID 4326, invalid WKT cases and zero-area rejection. This is parser regression coverage, not independent validation that all valid AOIs are supported.

## Required additions

- Add fixtures for antimeridian crossing, polar/high-latitude polygons, very small AOIs, long thin AOIs, near-collinear edges and 500/501 vertices.
- Validate projected AOI area and raster-mask coverage after transformation; record the exact geometry engine/library version.
- If holes/MULTIPOLYGON are later supported, version the schema and ensure every downstream rasterizer preserves topology exactly.

**Verdict:** PASS for the declared single-ring EPSG:4326 contract; unsupported geometries must fail explicitly.
