# SARGuardian

Multi-hazard Earth intelligence from SAR. Current focus: **glacier and landslide
deformation in the Nepal Himalaya**, built on NASA NISAR L-band interferometry.

---

## Why NISAR L2

NASA publishes NISAR **Level-2** products that are already coregistered,
unwrapped and geocoded:

| Product | What it gives you |
|---------|-------------------|
| `GUNW`  | Geocoded unwrapped interferogram → LOS displacement |
| `GOFF`  | Geocoded pixel offsets → large/fast motion, needs no coherence |
| `GCOV`  | Geocoded covariance → backscatter change detection |

So the first displacement time series needs **no SNAP, no ISCE2, no DEM, no
orbit files, no burst handling**. Drop to `RSLC` + ISCE2 only when you need a
custom pair network.

L-band (λ = 23.8 cm) also penetrates vegetation far better than Sentinel-1's
C-band (λ = 5.6 cm) and raises the unwrapping ceiling roughly fourfold.

**Constraint:** the NISAR archive starts mid-2025. For anything earlier, use
Sentinel-1 C-band.

---

## Areas of interest

| Name | Extent | Notes |
|------|--------|-------|
| **Langtang** | 80.8 km², 28.245–28.330 N | Langtang Lirung massif, 3,048–7,188 m, 80% above 4,000 m. Glacier monitoring; heavily studied, good validation literature. |
| **Lhende Khola** | 85.44–85.62 E, 28.34–28.47 N | Reported source zone of the 26 Aug 2026 collapse, ~9 km north of the Langtang box. |

Both are covered by NISAR paths **48 (descending)** and **98 (ascending)**.
Switch with `AOI_RING` at the top of each script.

Sentinel-1 tracks over the same ground: ASC 85 (frame 88), DESC 19 (frame 497),
DESC 121 (frames 498–499). All 12-day repeat. **Track 19 descending** is the one
to use — single frame, best geometry.

---

## Setup

```bash
python -m venv .venv
.venv\Scripts\activate        # Windows;  source .venv/bin/activate elsewhere
pip install -r requirements.txt
```

Credentials — either a `.env` from `.env.example`, or `~/.netrc`
(`~/_netrc` on Windows), which most geospatial tools already expect:

```
machine urs.earthdata.nasa.gov login YOUR_USER password YOUR_PASS
```

Register at <https://urs.earthdata.nasa.gov>. **Never commit either file.**

---

## Usage

Catalogue search is stdlib-only and needs no credentials, so recon runs anywhere.

```bash
# What exists over the AOI
python src/nisar_acquisition.py --recon --sentinel1

# Fetch the interferograms
python src/nisar_acquisition.py --download GUNW GOFF

# Confirm the HDF5 layout of your first real file
python src/gunw_reader.py --inspect data/nisar_l2/NISAR_L2_PR_GUNW_*.h5

# One pair → displacement, referenced to stable ground
python src/gunw_reader.py --read FILE.h5 --ref-lat 28.21 --ref-lon 85.47 \
    --geotiff outputs/d.tif --quicklook outputs/d.png

# Every pair → statistics table
python src/gunw_reader.py --batch data/nisar_l2 --ref-lat 28.21 --ref-lon 85.47 \
    --csv outputs/timeseries.csv
```

### Watching for the co-event pair

The 26 Aug 2026 collapse happened between acquisitions. The interferogram
spanning it publishes on the next 12-day repeat. `--watch` exits **10** when
something new appears, so cron can alert:

```bash
0 */6 * * * cd /path/to/sarguardian && python src/nisar_acquisition.py \
    --watch --new-since 2026-08-27 || notify-send "NISAR co-event pair available"
```

---

## Two things to verify on real data

**Sign convention.** `d_los = -(λ/4π)·φ`, positive = away from satellite.
Conventions differ between products and versions. Check against a signal whose
direction you already know before interpreting anything; `--flip-sign` inverts.

**Layover/shadow polarity.** The reader assumes `0 = good`. If the mask is
inverted you will gate out everything usable — obvious immediately from the
gate table.

Also: **unwrapped phase is relative.** Each connected component carries an
arbitrary constant, so absolute displacement is meaningless without
`--ref-lat/--ref-lon` on ground you believe is stable. Pick bedrock outside the
AOI, away from ice and the valley floor, and check its coherence first.

---

## Detectability: how much revisit do you need?

`src/detectability.py` simulates accelerating creep to failure (Voight,
alpha = 2), runs an inverse-velocity detector (Fukuzono) over samples taken at
a fixed revisit interval, and measures how much warning survives. It needs no
data.

```bash
python src/detectability.py --sweep --plot outputs/detectability.png
python src/detectability.py --demo --precursor 10 --revisit 4
```

Result, across 600 noise realisations per cell (5 mm LOS noise, 300 mm creep):

| Precursor | Detection collapses at | Usable revisit |
|-----------|------------------------|----------------|
| 5 days    | >= 2 days              | <= 1.7 days    |
| 10 days   | >= 4 days              | <= 3.3 days    |
| 20 days   | >= 8 days              | <= 6.7 days    |
| 40 days   | >= 16 days             | <= 13.3 days   |

**Revisit must be shorter than about one third of the precursor duration.**
The collapse lands at precursor/2.5 in every case tested - a sharp cliff, not
a gentle decline, because the detector needs at least three velocity estimates
inside the accelerating phase to fit a trend.

Applied to a Blatten-class 10-day precursor:

| Configuration | Revisit | Detection |
|---------------|---------|-----------|
| Single Sentinel-1 track | 12 d | **never** |
| NISAR + Sentinel-1, all 5 tracks combined | 4 d | 22% - marginal |
| Commercial tasking (ICEYE / Capella / Umbra) | 1 d | 100%, ~6 days warning |

This is the argument for the pipeline being sensor-agnostic: **revisit is a
procurement decision, not a scientific limit.** A national agency that tasks
daily commercial SAR gets operational early warning from the same code.

Counter-intuitive detail worth keeping: velocity noise scales as
`sigma*sqrt(2)/dt`, so a *longer* revisit gives a *quieter* velocity estimate.
Sampling density still wins, but the trade-off is why the cliff is sharp rather
than gradual.

**These are model results**, and the sweep above assumes measurement is always
possible. It is not - see the Blatten calibration below, where the phase
unwrapping ceiling removes the signal exactly when the slope accelerates.

---

## Impoundment susceptibility

`src/impoundment.py` maps where a landslide could dam a river and how much
water it would hold. Terrain only - no SAR, no credentials.

The 26 August 2026 cascade was not a single failure: ice/rock detachment ->
channel blockage -> breach -> surge. The blockage is the multiplier that turns
a local slope failure into a downstream flood, and blockage potential is a
property of the terrain, so it can be mapped in advance.

Priority-flood depression filling -> D8 flow routing -> flow accumulation ->
channel extraction -> flood the upstream contributing area of each channel cell
to a range of dam heights -> rank, with non-maximum suppression so one valley
reach cannot fill the table.

```bash
python src/impoundment.py --api --aoi lhende --heights 25 50 100 150     --rank-by efficiency --geojson outputs/dam_sites_lhende.geojson
```

Default ranking is **volume per metre of blockage**, not maximum volume. Ranking
by maximum volume just returns the largest dam height every time; what matters
operationally is which sites impound a lot for a *small* blockage, because a
25 m dam is a common event and a 150 m dam is not.

Preliminary result, SRTM 30 m sampled to ~105 m cells:

| AOI | Best site | Mm3/m | 25 m dam | 150 m dam | Sites >2 Mm3 at 25 m |
|-----|-----------|-------|----------|-----------|----------------------|
| Langtang | 28.2657 N, 85.5649 E | 0.10 | 2.4 Mm3 | 53.6 Mm3 | 3 of 10 |
| Lhende Khola | 28.4180 N, 85.5577 E | 0.21 | 5.1 Mm3 | 78.9 Mm3 | 8 of 8 |

**The Lhende Khola catchment is roughly twice as dammable as Langtang**, and
every site tested there impounds Thame-scale water (>2 Mm3) from only a 25 m
blockage. That is consistent with a cascade being possible there, though it is
terrain susceptibility, not a reconstruction of what happened.

Cross this layer with InSAR deformation on the flanking slopes and you get a
two-factor alert: a site that is both **dammable** and **moving**. That is the
actual thesis of SARGuardian, and it is defensible in a way a single AI risk
score is not.

**Limits, and state them in the paper.** Pools are not tested for spilling over
cols, so volumes are upper bounds near divides. Dam height is imposed, not
predicted - couple to a runout model before calling any number a forecast. SRTM
dates from 2000, so terrain reshaped since (the 2015 Gorkha avalanche through
Langtang) is out of date.

---

## Multi-geometry merge

`src/geometry_merge.py` combines ascending and descending LOS series into one
denser record of downslope motion.

**You can never interfere ascending with descending** - different look
geometry, zero coherence. But once each geometry has its own inverted LOS
series you are combining *measurements*, not phase, and those merge fine.

InSAR sees only the component along the line of sight. On a hillslope the
standard fix is to assume motion is downslope, take that direction from the
DEM, and solve for its magnitude:

    d_downslope = d_los / (slope_hat . los_hat)

That denominator is the **sensitivity**. Near zero, the track is blind and
dividing by it amplifies noise without limit. The tool computes it per track
and rejects geometries below a threshold.

Heading comes from orbital inclination and latitude
(`sin(heading_asc) = cos(i)/cos(lat)`), so no per-scene metadata is needed to
plan.

```bash
# Which tracks can even see motion here? Needs no data.
python src/geometry_merge.py --sensitivity --lat 28.29 --lon 85.51

# Merge the per-geometry series from timeseries.py
python src/geometry_merge.py --merge --ts outputs/ts.csv     --lat 28.29 --lon 85.51 --csv outputs/merged.csv --plot outputs/merged.png
```

Example, a north-facing 17-degree slope at 28.29 N 85.51 E:

| Track | Heading | Sensitivity | Noise x | Verdict |
|-------|---------|-------------|---------|---------|
| S1 ASC 85 | 350.7 | -0.437 | 2.3 | usable |
| NISAR ASC 98 | 350.5 | -0.437 | 2.3 | usable |
| NISAR DESC 48 | 189.5 | -0.225 | 4.5 | **blind** |
| S1 DESC 19 | 189.3 | -0.215 | 4.7 | **blind** |
| S1 DESC 121 | 189.3 | -0.215 | 4.7 | **blind** |

Two of five tracks usable, giving ~6-day combined sampling instead of 12.

### Sensitivity is not the same thing as layover

Earlier analysis found descending is the better geometry here for **visibility**
- 1.9% of the AOI in layover/shadow against 9.3% for ascending. This table says
ascending is better for **sensitivity** on this particular slope. Both are true
and they are different questions: layover asks whether the sensor can image the
ground at all, sensitivity asks whether the motion projects into the line of
sight once it can. Check both before choosing a track.

### The cancellation trap

A slope facing the satellite at close to the incidence angle has **near-zero**
sensitivity, because horizontal approach and vertical drop cancel in range. At
39 degrees incidence on a descending pass, an ESE-facing 35-degree slope gives
sensitivity +0.07 - effectively invisible - while a WNW-facing slope of the same
steepness gives -0.96, near ideal. A fast-moving slope can therefore show
nothing at all. This is the quantitative version of the geometry warning in the
methods brief.

When ascending and descending have **opposite** sensitivity sign at a site, that
is a genuine cross-check: real downslope motion must appear with opposite LOS
sign in the two geometries, and atmospheric artefacts will not do that.

---

## Calibration against Blatten (28 May 2025)

Published velocities for the Birch Glacier before it destroyed Blatten:

| Date | Days to failure | Velocity |
|------|-----------------|----------|
| 14 May 2025 | 14 | instability first observed on Kleines Nesthorn |
| 19 May 2025 | 9 | ~300 residents evacuated |
| 21 May 2025 | 7 | flow speed begins rising |
| 22 May 2025 | 6 | 0.5-0.8 m/day |
| 24 May 2025 | 4 | 4-4.5 m/day |
| 27 May 2025 | 1 | **10 m/day** |
| 28 May 2025 15:24 | 0 | collapse, ~9 Mm3 |

A Voight fit to those four velocity points gives alpha ~ 1.5, but with only
four published, rounded values that is indicative rather than authoritative.
The number that matters is the **scale**: metres per day, not millimetres.

### The finding: phase-based InSAR cannot measure a failing slope

Interferometric phase is ambiguous once displacement between passes exceeds
lambda/4. That sets a hard velocity ceiling:

| Band | lambda | lambda/4 | ceiling at 12 d | at 4 d | at 1 d |
|------|--------|----------|-----------------|--------|--------|
| L (NISAR) | 23.8 cm | 6.0 cm | 5.0 mm/day | 14.9 mm/day | 59.6 mm/day |
| C (Sentinel-1) | 5.5 cm | 1.4 cm | 1.2 mm/day | 3.5 mm/day | 13.9 mm/day |

Blatten was at **650 mm/day six days before failure** - already 10x above even
the 1-day L-band ceiling. Simulating the calibrated event both ways:

| Measurement | Revisit | Saturation | Detection | Warning |
|-------------|---------|------------|-----------|---------|
| Phase, L-band | 1 d | **100%** | 10% | unreliable (pred. error 21 d) |
| Phase, L-band | 4 d | **100%** | never | - |
| Offset tracking | 1 d | 0% | **100%** | **4.9 days** |
| Offset tracking | 2 d | 0% | **100%** | 2.6 days |
| Offset tracking | 4 d | 0% | 9% | - |

**For a Blatten-class failure, interferometry is the wrong instrument.** Phase
saturates the moment the slope starts running away. Offset tracking - GOFF, not
GUNW - is the only satellite measurement that survives to failure, and it needs
1-2 day revisit to give useful warning.

The handful of "detections" in the phase rows fire during the pre-saturation
noise and predict failure ~21 days out for a 7-day precursor. They are false
alarms, not warnings; the prediction error column is what exposes them.

### Two regimes, two instruments

The ESA analysis of ALOS-2 and SAOCOM L-band found the Kleines Nesthorn flank
creeping years ahead of failure: ~50 cm/yr by 2023, >150 cm/yr by August 2024.
That is 1.4-4.1 mm/day - comfortably **inside** the 12-day L-band ceiling.

| | Timescale | Velocity | Instrument | Answers |
|---|-----------|----------|------------|---------|
| **Site identification** | years | mm/day | InSAR phase, 12-day | *where* to watch |
| **Failure timing** | days | m/day | offset tracking, 1-2 day | *when* to evacuate |

Blatten was visible to phase InSAR for **years** before it failed, and invisible
to it during the fortnight that actually mattered. Both statements are true and
a credible warning system needs both instruments. This is the strongest argument
in the project for treating revisit and product type as procurement decisions:
a national agency buying daily commercial SAR with offset tracking gets the
second regime, which no free mission currently provides.

```bash
# Blatten calibrated, phase-limited
python src/detectability.py --sweep --precursor 7 --creep 27000     --noise 5 --wavelength 0.2384 --revisit 1 2 4 6 12

# same event, offset tracking
python src/detectability.py --sweep --precursor 7 --creep 27000     --noise 300 --revisit 1 2 4 6 12
```

Sources: published velocity timeline via Swiss cantonal reporting; long-term
creep rates from the ESA EO4Society analysis of ALOS-2 / SAOCOM. Neither the
raw monitoring series nor the InSAR time series is public, so the calibration
rests on reported summary figures.

---

## GOFF: offset tracking

`src/goff_reader.py` reads NISAR L2 pixel offsets. Use it wherever the slope is
moving too fast for phase.

```bash
python src/goff_reader.py --inspect FILE.h5
python src/goff_reader.py --read FILE.h5 --quicklook outputs/goff.png
python src/goff_reader.py --noise-floor data/nisar_l2/GOFF/2025-12_winter
```

Products carry `slantRangeOffset` and `alongTrackOffset` in **metres**, plus
per-pixel variance, `correlationSurfacePeak` and `snr`, across **three** layers
at different correlation window sizes. NASA labels them *raw, unculled,
unfiltered*, so gating, deramping and outlier rejection are the caller's job.

Processing chain: correlation and SNR gates, AOI clip, auto-reference on the
best fully-valid block outside the AOI, **planar deramp** fitted on stable
ground outside the AOI, then MAD-based culling.

Deramping matters. A stable winter pair showed a -646 mm median before it and
-146 mm after: residual orbit and coregistration error leaves a tilt across the
field that a single constant cannot remove.

### Measured noise floor

Winter pairs over a static glacier, so the scatter in the referenced field *is*
the detection threshold:

| Layer | Range sigma | 3-sigma floor |
|-------|-------------|---------------|
| layer1 (fine window) | 171 mm | 43 mm/day |
| layer2 | 114 mm | 29 mm/day |
| **layer3 (coarse window)** | **72 mm** | **18 mm/day** |

Best-case pairs reach 18 mm/day; the median across all five winter pairs is
75 mm/day, so quality varies by a factor of four between acquisitions. Quote
the median, not the best.

**The blind band is therefore about 5 to 18 mm/day at best, 5 to 75 mm/day
typically** - above the GUNW phase ceiling, below GOFF sensitivity. That is
measured, not assumed, and it is much narrower than the 5-125 mm/day I had
estimated before the data arrived.

### Polarisation

Quad-pol (`QD`) products carry both HH and VV. In these scenes **VV has no
valid pixels at all** while HH is fine, so layers are keyed by polarisation and
layer together - keying on the layer alone lets VV silently overwrite HH and
the whole read collapses to nothing.

---

## Inverse-velocity forecasting on measured data

`src/inverse_velocity.py` is the step that turns the pipeline into an alarm.
Fukuzono's construction: for accelerating creep, `1/v` falls linearly and
crosses zero at the failure time, so fitting a trailing window and
extrapolating the x-intercept gives a predicted failure date.

```bash
python src/inverse_velocity.py --ts outputs/ts_goff_summer.csv     --noise-floor 75 --event-date 2026-08-26 --plot outputs/inverse_velocity.png
```

It differs from the simulation in `detectability.py` in one decisive way: a
**significance gate**. Every velocity must exceed a multiple of the *measured*
noise floor before it may enter a fit. Without that an inverse-velocity fit
will happily forecast failure from three noise samples. Measure the floor
first with `goff_reader.py --noise-floor`.

When no alarm fires it reports which gate stopped it and what velocity would
have been required. A null with that bound attached is publishable; a null
without it is silence.

### First result on real data

GOFF layer3, both geometries, Lhende Khola box, July-August 2026:

| Geometry | Intervals | Fastest measured | Required | Verdict |
|----------|-----------|------------------|----------|---------|
| ASC 98 | 3 | +1.28 mm/day | 75 mm/day | no alarm |
| DESC 48 | 2 | +1.55 mm/day | 75 mm/day | no alarm |

**Shortfall: 48x.** No motion above 75 mm/day in the six weeks before the
collapse. That is a bounded non-detection, and the bound is the result.

### The finding that follows from it

Apparent motion is ~1.5 mm/day. That is:

- **below** the GOFF detection floor of 75 mm/day, so offsets cannot see it
- **below** the GUNW phase ceiling of 5.1 mm/day, so interferometry *can*

So for this site and this period, **GUNW is the correct product and GOFF is
the wrong one.** The two-regime argument is no longer theoretical: the ground
here sits squarely in the phase regime, and the offset product was never going
to resolve it. Prioritise the GUNW summer pairs.

---

## Working with a teammate who holds the data

Products are 1-2.4 GB each. Derived results are not. `--export` writes only the
AOI clip, so a 2.4 GB GUNW becomes a **0.1-0.2 MB** GeoTIFF carrying LOS
displacement and coherence, georeferenced and tagged with its pair dates.

Measured on the three winter GUNW: **488 KB total from 6.6 GB of source.**
Thirteen pairs come to roughly 2 MB - small enough to email.

### What the teammate runs

```bash
git pull
pip install -r requirements.txt          # pyproj and rasterio are required

python src/organise.py --src <wherever their files are>
python src/organise.py --src <wherever their files are> --apply

python src/gunw_reader.py --inspect data/nisar_l2/GUNW/*/NISAR_L2_*GUNW*.h5

python src/gunw_reader.py --batch data/nisar_l2/GUNW     --aoi lhende --auto-ref     --export outputs/export --csv outputs/gunw_stats.csv
```

Then zip `outputs/` and send it back. That is the whole handoff: **the .h5 files
never move.**

Run it twice, once per AOI, if both boxes are still in play - change `--aoi` and
the export directory.

---

## The noise floor is seasonal, and not in the direction you would guess

Same AOI (Lhende), same product, same layer3:

| Season | Range sigma | 3-sigma floor |
|--------|-------------|---------------|
| Winter (Nov-Dec) | 189 mm | **47 mm/day** |
| Monsoon (Jul-Aug) | 38 mm | **7 mm/day** |

**Monsoon is six times better than winter.** Snow, not rain, is what destroys
offset-tracking correlation in high mountains: fresh cover changes the surface
texture between passes, while in summer you are tracking bare rock and ice that
holds its speckle.

This matters twice. It reverses the usual assumption that monsoon acquisitions
are the poor ones. And it nearly closes the gap between products:

| | Velocity |
|---|---|
| GUNW phase ceiling | 5.1 mm/day |
| GOFF monsoon floor | 7.3 mm/day |
| **Blind band** | **5.1 - 7.3 mm/day, a factor of 1.4** |

Against the 5-125 mm/day I assumed before any data existed. In monsoon
conditions the two products very nearly meet, so almost nothing is invisible to
both - and the earlier winter-based floor of 75 mm/day was pessimistic by an
order of magnitude because it was measured in the wrong season and on the wrong
AOI.

Revised bound on the Nepal non-detection: **no motion above 7.3 mm/day** in the
six weeks before the collapse, a shortfall of 4.7x rather than 48x.

---

## Method notes

Two hazard classes, different physics, different tooling:

- **Class A — precursory creep.** Bare rock, moraine, ice, above treeline.
  InSAR time series works; a real forecast with lead time is possible.
- **Class B — rainfall-triggered.** Vegetated tropical slopes. InSAR early
  warning is *not* possible; use susceptibility (terrain + soil moisture +
  rainfall) plus post-event coherence-change mapping.

The Himalayan AOIs above are Class A, which is why they were chosen.

Data spine is NASA: NISAR (L-band), NASADEM (topography), GPM IMERG (rainfall),
SMAP (antecedent soil moisture), Landsat and FIRMS (optical, thermal).
Sentinel-1/2 remain as the partner-agency historical archive.

---

## Layout

```
src/     nisar_acquisition.py   catalogue search + download
         gunw_reader.py         GUNW → LOS displacement
data/    nisar_l2/  dem/        products (gitignored)
outputs/                        GeoTIFFs, quicklooks, CSVs (gitignored)
notebooks/  docs/
```
