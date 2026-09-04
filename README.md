# SARGuardian

Multi-hazard Earth intelligence from SAR. Built on NASA NISAR L-band
interferometry, tested against the **26 August 2026 Langtang massif collapse**.

---

## The finding

We set out to detect a Himalayan slope failure with NISAR. We could not, and
the reason is the result.

**Spaceborne L-band SAR has four regimes over this terrain, and the failure
falls in the one where both products stop working.** Interferometric phase
measures slow creep and only in winter, because in the monsoon it retains 0-3%
of its pixels. Offset tracking survives the monsoon but its floor is 19 mm/day,
a thousand times slower than a failing slope. When the slope actually goes, the
surface is destroyed, correlation collapses, and offset tracking loses the
ground it was tracking - which is itself the only signal either product records
of the event.

So the answer to "can NISAR give warning of a collapse like this one" is: not
at 12-day repeat, not with L2 products, and we can say precisely which limit
stopped each one. **No motion above 18.6 mm/day in the seven weeks before
failure**, on the one geometry with a stable floor, with the last observation
seven days out. A Blatten-class precursor would have exceeded that floor by
16-26x. It was not there.

That is a bounded null with a measured floor behind it, paired with a measured
positive - and it is an argument about instruments and revisit, not about this
one mountain.

---

## Areas of interest

| Name | Extent | Role |
|------|--------|------|
| **Source zone** | 28.2453-28.3529 N, 85.4645-85.5562 E | **The analysis AOI.** Confirmed source zone of the 26 Aug 2026 collapse; contains the failure point at 28.28771 N, 85.52809 E. |
| Langtang | 28.2447-28.3297 N, 85.4591-85.5649 E | The wider massif box. Contains the source zone; 78.4% overlap with it. Glacier monitoring, good validation literature. |
| Lhende Khola | 28.3400-28.4700 N, 85.4400-85.6200 E | **Control region only.** See the correction below. |

Select at runtime with `--aoi source`, `--aoi langtang` or `--aoi lhende`.
Both boxes are covered by NISAR paths **48 (descending)** and **98
(ascending)**.

Sentinel-1 tracks over the same ground: ASC 85 (frame 88), DESC 19 (frame 497),
DESC 121 (frames 498-499), all 12-day repeat.

### A correction, because it invalidated four months of results

An earlier version of this README carried a warning in bold: *"the Langtang box
does not contain the 26 Aug failure zone - that sits ~9 km north, in the Lhende
Khola catchment."* **That was wrong.** The Lhende extent came from a report
published immediately after the collapse and was never checked against
anything.

The confirmed source zone is **5.8 km south of the southern edge of the Lhende
box** and overlaps it by 12%. Every result computed with `--aoi lhende` - the
bounded non-detection, the detection floors, the product-choice conclusion, the
co-event analysis - described ground that does not contain the failure. Those
results were not wrong so much as vacuous: they correctly characterise a piece
of the Himalaya where nothing happened.

We keep Lhende as a labelled control region, because a well-characterised patch
of nearby ground where nothing occurred is genuinely useful. It is no longer
the analysis AOI.

---

## The four regimes

Measured over the source zone, on real products, except where marked.

| Regime | Velocity | GUNW - phase | GOFF - offsets | Status |
|--------|----------|--------------|----------------|--------|
| Winter, slow creep | mm/day | **Works.** 3% coverage, below the phase ceiling | Blind. Floor 9-118 mm/day | measured |
| Monsoon, slow creep | mm/day | **Fails.** 0% coverage - decorrelation | Marginal. ASC floor 9-25 mm/day | measured |
| Pre-failure acceleration | m/day | Fails. Above the ceiling, and decorrelates | Works at 1-2 day revisit | simulated |
| **The failure itself** | metres | **Fails.** 0-1% coverage | **Loses the surface it tracks** | **measured** |

The two products have **opposite** seasonal behaviour, which is the physically
sensible answer and not the one we expected. Offsets track amplitude speckle
and improve modestly in the monsoon; phase needs coherence and collapses in it.

---

## Why NISAR L2

NASA publishes NISAR **Level-2** products that are already coregistered,
unwrapped and geocoded:

| Product | What it gives you |
|---------|-------------------|
| `GUNW`  | Geocoded unwrapped interferogram -> LOS displacement |
| `GOFF`  | Geocoded pixel offsets -> large/fast motion, needs no coherence |
| `GCOV`  | Geocoded covariance -> backscatter change detection |

So the first displacement time series needs **no SNAP, no ISCE2, no DEM, no
orbit files, no burst handling**.

L-band (lambda = 23.84 cm) penetrates vegetation far better than Sentinel-1's
C-band (5.55 cm) and raises the unwrapping ceiling roughly fourfold.

**Two constraints worth knowing before you plan anything.** The NISAR archive
starts mid-2025; for anything earlier, use Sentinel-1 C-band. And NISAR L2 ships
**one interferogram per consecutive acquisition pair** - the catalogue over this
AOI returns 17 GUNW pairs and every one is between consecutive acquisitions.
There are no loop-closing pairs, so **every L2-only time series is a chain with
zero redundancy and no internal error estimate.** That is not specific to this
site. Drop to `RSLC` + ISCE2 when you need a network that can check itself.

---

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate     # Linux/macOS;  .venv\Scripts\activate on Windows
pip install -r requirements.txt
```

On most Linux distributions the command is `python3`, not `python`.

A `.venv` inside the repository is fine - the tools skip virtualenvs when
searching for data. Without that they would report h5py's own bundled test
files (`vlen_string_dset.h5` and friends) as NISAR products.

Credentials - either a `.env` from `.env.example`, or `~/.netrc`
(`~/_netrc` on Windows), which most geospatial tools already expect:

```
machine urs.earthdata.nasa.gov login YOUR_USER password YOUR_PASS
```

Register at <https://urs.earthdata.nasa.gov>. **Never commit either file.**

---

## Usage

Every script resolves paths against the **repository root**, not your current
directory, so the same command works from anywhere. Run `python src/paths.py`
to see where the code thinks everything is and how many products it can find.

### 1. Find out what exists (no credentials, no downloads)

```bash
python src/nisar_acquisition.py --recon --sentinel1
```

### 2. Download, then organise

```bash
python src/nisar_acquisition.py --download GUNW GOFF
python src/organise.py            # dry run: finds files anywhere in the repo
python src/organise.py --apply    # then move them
```

Products land in `data/nisar_l2/_incoming/`. `organise.py` searches the whole
repository recursively, leaves files already in place alone, and writes
`data/nisar_l2/MANIFEST.csv`. Re-running is always safe.

### 3. Confirm the layout of any new product type

```bash
python src/gunw_reader.py --inspect data/nisar_l2/GUNW/2025-11_winter/NISAR_L2_PR_GUNW_*.h5
python src/goff_reader.py --inspect data/nisar_l2/GOFF/co_event/NISAR_L2_UR_GOFF_*.h5
```

Do this once per product type. Every real-data bug so far was found this way.

### 4. Read

```bash
# one pair
python src/gunw_reader.py --read FILE.h5 --aoi source --auto-ref \
    --quicklook outputs/d.png

# every GUNW, with AOI-clipped export small enough to email
python src/gunw_reader.py --batch --aoi source --auto-ref \
    --export outputs/export_src --csv outputs/gunw_stats_source.csv

# offsets, and the measured detection floor
python src/goff_reader.py --batch --aoi source --layer layer2 \
    --csv outputs/goff_stats_source.csv
python src/goff_reader.py --noise-floor data/nisar_l2/GOFF/2026-07_summer --aoi source
```

`--batch` defaults to `data/nisar_l2/GUNW` (or `GOFF`) and searches recursively.
**Always use `--auto-ref`** - see the reference note below.

### 5. Coherence-change mapping

Every export lands on a **fixed AOI grid**: cell edges anchored to absolute
multiples of the pixel size in the projected CRS, and the AOI bounds snapped
outward onto them. Over the source zone at 80 m posting that is 116 x 155,
origin (349280, 3137440), EPSG:32645 - for **every** product, GUNW and GOFF,
ascending and descending, pre-event and co-event. Products fill what they
cover; the rest is nodata.

So two rasters of the same ground subtract directly, with no cropping and no
half-pixel guess:

```python
import rasterio
a = rasterio.open("outputs/export_src/GUNW_20251128_20251210_PR.tif")
b = rasterio.open("outputs/export_src/GUNW_20251210_20251222_PR.tif")
change = b.read(2) - a.read(2)        # coherence band, 5,843 common pixels
```

**Then the caveat that matters more than the fix.** Mechanically this now works
for any pair. Scientifically it only works in winter here: the best winter
pairing shares **6,467** valid pixels, while every summer pairing shares
**between 0 and 6**. Coherence-change mapping of the 26 August failure is
therefore not available from GUNW at this site - not because the pipeline
cannot do it, but because there is no monsoon coherence to change. The
equivalent measurement on GOFF correlation does work, and is what the co-event
detection above rests on.

That is the standard way to map a failure after the fact, and it is the method
these notes recommend for Class B hazards. It was unreachable from this pipeline
until the grid was fixed, because each export was clipped to the bounding box of
its own valid pixels and came out a different size every time - 179x220 for one
pair, 181x142 for the next. Nothing failed; the files were correct and simply
could not be compared.

A product whose own grid does not share the lattice is **refused, not
resampled**: it falls back to a per-product clip and warns that the file will
not align. Resampling would invent values and hide the problem.

> Exports written before this change are on the old per-product grids. Re-run
> `--export` on any directory you intend to difference.

### 6. Time series, then forecast

```bash
python src/timeseries.py --dir --product GOFF --network          # structure only
python src/timeseries.py --dir data/nisar_l2/GOFF --product GOFF \
    --goff-layer layer2 --aoi source --invert --auto-ref --jackknife \
    --csv outputs/ts_goff_source.csv

python src/inverse_velocity.py --ts outputs/ts_goff_source.csv \
    --floors outputs/goff_stats_source.csv --floors-layer layer2 \
    --noise-floor 18.6 --event-date 2026-08-26
```

**Pass `--floors`.** It gates each interval against the floor of the pair that
produced it instead of one scalar for the whole stack, and per-pair floors here
run from 9.4 to 117.9 mm/day. Without it, the descending interval
2026-06-29 -> 2026-07-11 reads +19.60 mm/day and clears a global 18.6 mm/day
gate at 1.05x; against its own 114.3 mm/day floor it is **0.17x**, plainly
inside the noise. `--noise-floor` is then only the fallback for pairs the
stats file does not cover.

`--invert` is required to do more than print the network. **Always pass
`--jackknife`** - see the redundancy note below. The noise floor passed to the
detector must be the one you **measured** for that product, not a guess.

Note: `--from-stats` only understands the GUNW column schema. With
`--product GOFF` it raises `KeyError: 'median'`. Use `--dir`.

### 7. Geometry and terrain (no data needed)

```bash
python src/geometry_merge.py --sensitivity --lat 28.2877 --lon 85.5281
python src/detectability.py --sweep --plot outputs/detectability.png
python src/impoundment.py --api --aoi source --rank-by efficiency \
    --geojson outputs/dam_sites.geojson
```

### Watching for a new product

```bash
0 */6 * * * cd /path/to/SARGuardian && python src/nisar_acquisition.py \
    --watch --new-since 2026-08-28 || notify-send "New NISAR product"
```

Exit code **10** means something new appeared.

---

## Six things that will bite you

**Unwrapped phase is relative.** Every connected component carries an arbitrary
constant, so an absolute displacement is meaningless until it is referenced.
**Use `--auto-ref`.** Hand-picking a reference does not work: on a real winter
scene all three "obvious" choices had *zero* usable pixels under snow, while a
block 8 km away sat at 0.95 coherence.

**The GUNW `mask` layer is not a layover/shadow flag.** It is a three-digit
code: hundreds = water, tens = subswath in the reference image, units =
subswath in the secondary. A usable pixel is dry land inside a real subswath in
**both** acquisitions. Keeping `mask == 0` keeps exactly the pixels that were
invalid in both - which took valid coverage from 39.5% down to 0.5% before this
was found.

**Sign convention.** `d_los = -(lambda/4pi) * phi`, positive = away from the
satellite. Check it against a signal whose direction you already know before
interpreting anything; `--flip-sign` inverts it.

**Every L2 time series has zero redundancy.** The network is a chain, so the
inversion fits every observation exactly and reports error bars that are not
error bars. On the ascending summer block, four of five "pairs" are load-bearing
and the fifth is the same acquisition pair processed twice - the residual RMS
of 0.58 mm measures how well two processings of identical data agree, not
measurement scatter. **Run `--jackknife` and believe what it says.** On the
descending winter block it says the trend reverses sign when one interferogram
is removed.

**Sensitivity is not determined on gentle terrain.** It is a projection onto the
downslope direction, and on nearly flat ground that direction is whichever way
the DEM noise tilts. At one target here, widening the DEM stencil from 60 m to
300 m moves the aspect 74 degrees and takes the ascending sensitivity from
+0.020 to -0.642, sign included. A 62-degree slope 3 km away holds its aspect to
4 degrees over the same range. **Do not quote sensitivity to three decimals
below about 10 degrees of slope**, and treat the multi-geometry table below as
indicative on shallow ground.

**Routine and urgent processing of the same acquisitions disagree.** For the
co-event GUNW pair, PR and UR give AOI medians of -85.93 mm and +251.00 mm - a
336.93 mm gap, 2.76 fringes. `gunw_reader.py` catches this itself and refuses
to average them. Keep the product that agrees with the rest of the stack.

---

## Evidence

### The pre-event bound

Ascending path 98, GOFF layer2, over the source zone. Summer block, 2 July to
19 August 2026 - the last observation seven days before failure.

| Interval | Days | Velocity | Against the 18.6 mm/day gate |
|----------|------|----------|------------------------------|
| 2026-07-02 -> 2026-07-14 | 12 | -3.39 mm/day | below floor, 0.18x |
| 2026-07-14 -> 2026-07-26 | 12 | +1.72 mm/day | below floor, 0.09x |
| 2026-07-26 -> 2026-08-19 | 24 | -0.41 mm/day | below floor, 0.02x |

Fitted linear velocity **-0.711 mm/day**, sign stable across every removable
subset (-0.723 to -0.455) but **untested** - the block has zero redundancy, so
we quote no interval. Nothing approaches the floor.

**The bound is the result, and the bound does not come from the time series.**
It comes from the MAD scatter of eight independent pairs, none of which depends
on the network:

| Track | Pairs | 3-sigma floor, mm/day | Usable? |
|-------|-------|----------------------|---------|
| **ASC 098** | 8 | median **19.8**, range **8.9-24.5** | Yes - stable to a factor of 2.8 across nine months |
| DESC 048 | 8 | median 98.3, range 9.4-117.9 | No - varies twelvefold |

Descending shows one interval at +19.60 mm/day that clears a global gate, but
that pair's own scatter is 169.7 mm. **Gate each interval against the floor of
the pair that produced it**, not against one scalar, and it disappears.

**Two coverage caveats, stated because they are the most attackable points.**
The final seven days before failure are unobserved. The interval that covers
late August is a 24-day average, which dilutes a 7-day precursor about
threefold.

### The co-event detection

The measurement that carries the event is not displacement - it is loss of
correlation. Largest connected region of new decorrelation between consecutive
pairs, GOFF layer2:

| Track | Event-free baseline (5 transitions each) | Co-event pair |
|-------|------------------------------------------|---------------|
| DESC 048 | largest blob 22-81 px | **2,801 px = 17.93 km2**, 54.9% coverage lost |
| ASC 098 | largest blob 25-50 px | **712 px = 4.56 km2**, 23.0% lost |

Both geometries record their largest decorrelation event of the entire archive
in the one pair that spans 26 August - **35x and 14x the event-free baseline**.
The reported failure point falls inside the descending footprint, whose centroid
is 0.94 km away.

Two controls:

- **Random-loss null.** Losing the same number of pixels at random from the same
  pre-valid mask, 200 trials: largest blob median 201 px, 95th percentile 365,
  maximum 523. The observed 2,801 px sits at percentile 100.
- **Processing chain.** The routine product for the identical pair loses 45.3%
  with a 2,266 px / 14.50 km2 blob. Same phenomenon in both chains, so it is in
  the data and not the processor.

**What this is not.** The footprint fills only 38% of its 6.81 x 6.90 km
bounding box - sprawling, not a compact scar - and the terrain under it
(3,700-5,819 m, slope 36.8 deg) is statistically indistinguishable from the AOI
background (3,553-6,522 m, 37.1 deg). It maps co-event surface disturbance, not
the failure scar. Late August is also peak monsoon, which decorrelates broadly
on its own; the preceding monsoon transition lost only 12.2%, which argues
against a purely seasonal explanation but does not eliminate one.

### Product coverage

Valid-pixel fraction over the source zone, all 15 GUNW pairs:

| Season | Median | Range |
|--------|--------|-------|
| Winter | 3% | 1-4% |
| Monsoon | **0%** | 0-1% |

`gunw_reader.py` prints the warning itself on every summer pair: *"under 10% of
the scene survived gating... Consider GOFF instead."* **There is no usable
summer phase measurement at the source zone.** The reasoning that once pointed
the other way asked whether the motion was below the phase ceiling - it is - and
never asked whether there was any phase to measure.

### The stratified troposphere

LOS displacement regressed against DEM elevation, per pair, on 698 samples
inside the source polygon spanning 3,075-7,166 m.

| Product | Pairs | \|r\| median | Variance explained | p<0.001 |
|---------|-------|-------------|--------------------|---------|
| GUNW - phase | 7 | **0.40** | **16.0%** | 5 / 7 |
| GOFF - offsets | 14 | **0.10** | **1.0%** | 3 / 14 |

The four best-sampled winter GUNW pairs give r = -0.49, +0.47, +0.40, -0.43 -
**alternating sign on consecutive 12-day pairs.** Ground does not reverse
direction every twelve days; a water-vapour field does. Over the 4,091 m of
relief the GUNW term reaches **+/-78 to +/-116 mm**, larger than the
quarter-wavelength ceiling of 59.5 mm per 12-day pair.

The sixteenfold difference between products is what the physics predicts. Phase
measures optical path length and is perturbed directly by stratified water
vapour; offset tracking measures a geometric pixel shift and is affected only at
second order. **The bounded null rests on GOFF, which is the product that
carries 1% of this.**

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

The number that matters is the **scale**: metres per day, not millimetres.

**On the Voight exponent** [1]**.** The published record does not constrain it. Three
numeric velocities at rounded dates give alpha ~ 1.77, but shifting the assumed
failure time by half a day moves the estimate across **1.62 to 2.07**. That is
why `detectability.py` assumes alpha = 2 rather than fitting it - the assumption
is honest and the fit would not be.

### Phase cannot measure a failing slope

Two independent limits, and the second is the decisive one.

**Ambiguity.** Interferometric phase is ambiguous once the displacement
*difference between adjacent pixels* exceeds lambda/4 - the Itoh condition [3]. That
sets a ceiling on the resolvable phase gradient:

| Band | lambda | lambda/4 | at 12 d | at 4 d | at 1 d |
|------|--------|----------|---------|--------|--------|
| L (NISAR) | 23.84 cm | 5.96 cm | 5.0 mm/day | 14.9 mm/day | 59.6 mm/day |
| C (Sentinel-1) | 5.55 cm | 1.39 cm | 1.2 mm/day | 3.5 mm/day | 13.9 mm/day |

**Decorrelation, which binds first and harder.** At 0.65 m/day over a 12-day
pair the surface moves roughly 7.8 m. That rearranges the scatterers inside
every resolution cell completely: coherence goes to zero and there is no phase
to unwrap at all. Blatten was at 650 mm/day six days before failure - eleven
times the 1-day L-band ceiling, and far past any coherence.

**For a Blatten-class failure, interferometry is the wrong instrument**, and it
is the wrong instrument twice over.

### Two regimes, two instruments

The ESA analysis of ALOS-2 and SAOCOM L-band found the Kleines Nesthorn flank
creeping years ahead of failure: ~50 cm/yr by 2023, >150 cm/yr by August 2024 -
1.4-4.1 mm/day, comfortably **inside** the 12-day L-band ceiling.

| | Timescale | Velocity | Instrument | Answers |
|---|-----------|----------|------------|---------|
| **Site identification** | years | mm/day | InSAR phase, 12-day | *where* to watch |
| **Failure timing** | days | m/day | offset tracking, 1-2 day | *when* to evacuate |

Blatten was visible to phase InSAR for **years** before it failed, and invisible
to it during the fortnight that actually mattered. Both statements are true and
a credible warning system needs both instruments.

---

## Detectability: how much revisit do you need?

`src/detectability.py` simulates accelerating creep, runs an inverse-velocity
detector over samples at a fixed revisit, and measures the warning that
survives. It needs no data.

```bash
python src/detectability.py --sweep --precursor 5 10 20 40
```

**Revisit must be shorter than about one third of the precursor duration.** The
collapse lands at precursor/2.5 in every case tested - a sharp cliff, because
the detector needs at least three velocity estimates inside the accelerating
phase.

| Precursor | Detection collapses at | Usable revisit |
|-----------|------------------------|----------------|
| 5 days | >= 2 days | <= 1.7 days |
| 10 days | >= 4 days | <= 3.3 days |
| 20 days | >= 8 days | <= 6.7 days |
| 40 days | >= 16 days | <= 13.3 days |

### Daily revisit is worse, not better

This is the counter-intuitive result and it survives every parameterisation we
have run. Velocity noise scales as `sigma*sqrt(2)/dt`, so a shorter revisit
gives a **noisier** velocity estimate. Reading the false-alarm column:

| Precursor | Revisit | Detection | False alarm | Prediction error |
|-----------|---------|-----------|-------------|------------------|
| 10 days | 1 d | 91% | **11.3%** | 4.7 d |
| 10 days | **3 d** | **95%** | **0.4%** | **2.0 d** |
| 20 days | 1 d | 90% | 14.6% | 8.2 d |
| 20 days | **4 d** | **98%** | **0.0%** | 5.0 d |
| 40 days | 1 d | 92% | **20.5%** | 11.8 d |
| 40 days | **6 d** | **96%** | **0.0%** | 8.0 d |

**Buy revisit near precursor/3. Below that, false alarms rise without improving
detection.** An agency tasking daily SAR on the strength of a detection rate
alone would evacuate on one alarm in five. This is the procurement conclusion,
and it is the opposite of the obvious one.

### The Blatten case, both instruments

```bash
# phase-limited
python src/detectability.py --sweep --precursor 7 --creep 27000 \
    --noise 5 --wavelength 0.2384 --revisit 1 2 4 6 12

# same event, offset tracking
python src/detectability.py --sweep --precursor 7 --creep 27000 \
    --noise 300 --revisit 1 2 4 6 12
```

| Measurement | Revisit | Saturation | Detection | False alarm | Warning |
|-------------|---------|------------|-----------|-------------|---------|
| Phase, L-band | 1 d | **100%** | never | 10.5% | - |
| Phase, L-band | 4 d | **100%** | never | 0.0% | - |
| Offset tracking | 1 d | 0% | 89% | 14.1% | 4.8 d (err 3.2 d) |
| Offset tracking | 2 d | 0% | 94% | 7.2% | 2.6 d (err 1.2 d) |
| Offset tracking | 4 d | 0% | 5% | 2.6% | - |

Phase saturates the moment the slope runs away. Offset tracking survives to
failure and needs 1-2 day revisit to give useful warning - but note that even
at 1 day the predicted failure date carries a 3.2-day error against a 7-day
precursor.

**These are model results.** Absolute warning times depend on the assumed
precursor duration, creep amplitude and measurement noise.

---

## GOFF: offset tracking

`src/goff_reader.py` reads NISAR L2 pixel offsets. Use it wherever the slope is
moving too fast for phase, or wherever phase has no coherence to work with.

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

### Track dominates everything

The strongest single control on data quality here is not season, product or
layer - it is **geometry**.

All measured over the source zone, so the comparison is like for like.

| Measurement | ASC 098 | DESC 048 | Ratio |
|-------------|---------|----------|-------|
| GOFF range sigma, layer2 (median of 8) | 79.2 mm | 393.2 mm | **5.0x** |
| GOFF 3-sigma floor (median of 8) | 19.8 mm/day | 98.3 mm/day | **5.0x** |
| GUNW valid pixels, winter (median) | 3.5% | 1.0% | **3.5x** |
| Co-event coverage lost | 23.0% | 54.9% | **2.4x** |

Two products, two different physical measurements, the same geometric control.
**Ascending path 98 is the usable geometry over this terrain**, and that is an
acquisition-planning result, not a footnote.

Season is real but secondary, and it only shows its true size once you stratify
by track. Over the source zone, GOFF layer2:

| | Winter sigma | Summer sigma | Seasonal ratio |
|---|---|---|---|
| ASC 098 | 90.4 mm | 66.7 mm | **1.36x** |
| DESC 048 | 434.9 mm | 343.6 mm | **1.27x** |
| Both pooled | 235.2 mm | 74.7 mm | 3.15x |

The pooled row is the trap. It reads as a threefold seasonal improvement, and
it is not: it is a **1.3x** seasonal effect riding on a **5.0x** geometric one,
produced because the two seasonal subsets are not balanced across tracks.
Control for the AOI, the product and the layer and you still get the wrong
answer if you do not control for geometry.

### Polarisation

Quad-pol (`QD`) products carry both HH and VV. In these scenes **VV has no
valid pixels at all** while HH is fine, so layers are keyed by polarisation and
layer together - keying on the layer alone lets VV silently overwrite HH and
the whole read collapses to nothing.

---

## Inverse-velocity forecasting on measured data

`src/inverse_velocity.py` turns the pipeline into an alarm. Fukuzono's
construction [2]: for accelerating creep, `1/v` falls linearly and crosses zero at
the failure time, so fitting a trailing window and extrapolating the x-intercept
gives a predicted failure date.

It differs from the simulation in one decisive way: a **significance gate**.
Every velocity must exceed a multiple of the *measured* noise floor before it
may enter a fit. Without that an inverse-velocity fit will happily forecast
failure from three noise samples.

When no alarm fires it reports which gate stopped it and what velocity would
have been required. **A null with that bound attached is publishable; a null
without it is silence.**

Two limits to state when you use it. The default window of three points leaves
one degree of freedom, so the reported interval on the failure date comes from a
delta-method propagation that is unreliable near a zero denominator - which is
exactly the regime the detector works in. And the gate takes a single scalar
floor, which is wrong when per-pair floors vary twelvefold; gate per pair.

---

## Multi-geometry merge

`src/geometry_merge.py` combines ascending and descending LOS series into one
denser record of downslope motion. You can never interfere ascending with
descending, but once each geometry has its own inverted LOS series you are
combining *measurements*, not phase.

    d_downslope = d_los / (slope_hat . los_hat)

That denominator is the **sensitivity**. Near zero the track is blind and
dividing by it amplifies noise without limit.

### At the failure point, four of five tracks are blind

```bash
python src/geometry_merge.py --sensitivity --lat 28.2877 --lon 85.5281
```

The failure point sits at 5,166 m on a **27.4-degree west-facing slope** -
steep enough that the sensitivity is well determined there, unlike the gentle
ground discussed above.

| Track | Heading | Sensitivity | Noise x | Verdict |
|-------|---------|-------------|---------|---------|
| S1 ASC 85 | 350.7 | **-0.904** | 1.1 | usable |
| **NISAR ASC 98** | 350.5 | **-0.890** | 1.1 | **usable** |
| S1 DESC 19 | 189.3 | +0.197 | 5.1 | **blind** |
| S1 DESC 121 | 189.3 | +0.197 | 5.1 | **blind** |
| NISAR DESC 48 | 189.5 | +0.163 | 6.1 | **blind** |

**Only the ascending tracks can see downslope motion at the place the slope
actually failed**, and they see it almost perfectly - a sensitivity of -0.890
amplifies noise by only 1.1x. Every descending track, NISAR and Sentinel-1
alike, is geometrically blind there: a west-facing slope viewed from the east
puts horizontal approach and vertical drop into near-cancellation in range.

This is why ascending is the usable geometry in every table above, and it is not
a data-quality accident - it is the terrain. It also means the descending
track's larger co-event coverage loss reflects its worse viewing geometry over
this slope, not that it saw more of the event.

**For a warning system, that is the finding.** Doubling your revisit by
combining ascending and descending is not available here; on this aspect you
have one usable look direction and a 12-day repeat, and no amount of processing
recovers the other one.

### The cancellation trap

A slope facing the satellite at close to the incidence angle has **near-zero**
sensitivity, because horizontal approach and vertical drop cancel in range. At
39 degrees incidence on a descending pass, an ESE-facing 35-degree slope gives
sensitivity +0.07 - effectively invisible - while a WNW-facing slope of the same
steepness gives -0.96, near ideal. A fast-moving slope can therefore show
nothing at all.

### Telling motion from atmosphere

Downslope motion has one magnitude projected onto two lines of sight, so

    v_desc / v_asc = sens_desc / sens_asc

is fixed by geometry alone. Where the sensitivities have opposite signs, real
motion **must** appear with opposite LOS signs. A path delay is not a vector -
atmosphere and snow add the same extra path length whichever direction the radar
looks from, so they arrive with the **same** sign and a ratio near +1.

`candidate_check.py` compares the measured ratio against both predictions and
reports which it is closer to. It refuses to decide in four cases: no
coordinates, a blind track, a slope below 10 degrees, or a geometry where the
two predictions coincide.

---

## Rejecting our own candidate

The strongest-looking signal in the dataset is at **28.27484 N, 85.47405 E**,
inside the source polygon, 5.47 km west-south-west of the failure point. It
accumulates across four consecutive ascending dry-season pairs (-54.3, -27.6,
-24.3, -61.5 mm on full 49-pixel windows), is spatially coherent over
0.56 km2, and its coherence is 0.68 against a scene mean of 0.615 - better than
its surroundings, not worse.

**It is not ground motion.** Four tests, and what each one says:

| Test | Verdict |
|------|---------|
| Seasonality | Suggestive. Dry-season signal, monsoon quiescence - but the monsoon windows carry 0-5 valid pixels, below the 20-pixel threshold. That is measurement failure, not measured quiescence. |
| Geometry ratio | **Undetermined.** The slope is 5.0 degrees at native DEM spacing; the sensitivity is not determined and neither is any ratio built on it. |
| Elevation plane | **Holds.** Removing the fitted stratified term changes the signal by -2%. Stratified troposphere is excluded as the explanation. |
| Terrain | **Holds.** 4-7 degrees at every stencil width. Nothing fails at 5 degrees. |

The honest description is **a geometry-independent path delay, most plausibly a
seasonal snowpack**. Dry snow is nearly transparent at L-band, so it preserves
the very coherence that makes the signal look trustworthy while adding a path
delay that grows with accumulated water equivalent. Turbulent water vapour with
a persistent local pattern is not excluded; separating it needs a weather model
or snow-depth data, not more SAR.

**Record the location and the reason. Do not carry it forward as a detection.**

---

## Impoundment susceptibility

`src/impoundment.py` maps where a landslide could dam a river and how much water
it would hold. Terrain only - no SAR, no credentials.

A collapse of this kind is rarely a single failure: detachment -> channel
blockage -> breach -> surge. The blockage is the multiplier that turns a local
slope failure into a downstream flood, and blockage potential is a property of
the terrain, so it can be mapped in advance.

Priority-flood depression filling -> D8 flow routing -> flow accumulation ->
channel extraction -> flood the upstream contributing area of each channel cell
to a range of dam heights -> rank, with non-maximum suppression so one valley
reach cannot fill the table.

Default ranking is **volume per metre of blockage**, not maximum volume. Ranking
by maximum volume just returns the largest dam height every time; what matters
operationally is which sites impound a lot for a *small* blockage, because a
25 m dam is a common event and a 150 m dam is not.

Preliminary result, SRTM 30 m sampled to ~105 m cells:

| AOI | Best site | Mm3/m | 25 m dam | 150 m dam | Sites >2 Mm3 at 25 m |
|-----|-----------|-------|----------|-----------|----------------------|
| Langtang | 28.2657 N, 85.5649 E | 0.10 | 2.4 Mm3 | 53.6 Mm3 | 3 of 10 |
| Lhende Khola | 28.4180 N, 85.5577 E | 0.21 | 5.1 Mm3 | 78.9 Mm3 | 8 of 8 |

**Not yet run on the source zone.** Both rows above predate the AOI correction.
`--aoi source` is registered and the run is one command, but until it is done
this table describes the wider massif box and the control region, not the ground
that failed. Do not quote it as if it did.

Two further AOIs are registered and unrun: `source`, and `blatten` - the
Loetschental, where the Birch Glacier deposit obstructed the Lonza and impounded
a lake on 28 May 2025. That second one is the only available test of whether any
of this travels outside the Himalaya, and it has a documented outcome to check
against.

Cross this layer with InSAR deformation on the flanking slopes and you get a
two-factor alert: a site that is both **dammable** and **moving**. That is the
actual thesis of SARGuardian, and it is defensible in a way a single risk score
is not.

**Limits, and state them in the paper.** Pools are not tested for spilling over
cols, so volumes are upper bounds near divides. Dam height is imposed, not
predicted - couple to a runout model before calling any number a forecast. SRTM
dates from 2000, so terrain reshaped since is out of date. And note that the box
containing the confirmed source zone is the *less* dammable of the two tested;
this is terrain susceptibility, not a reconstruction of what happened.

---

## Working with a teammate who holds the data

Products are 1-2.4 GB each. Derived results are not. `--export` writes only the
AOI clip, so a 2.4 GB GUNW becomes a **0.1-0.2 MB** GeoTIFF carrying LOS
displacement and coherence, georeferenced and tagged with its pair dates.
Measured on the three winter GUNW: **488 KB total from 6.6 GB of source.**

```bash
git pull
pip install -r requirements.txt          # pyproj and rasterio are required

python src/organise.py --src <wherever their files are>
python src/organise.py --src <wherever their files are> --apply
python src/gunw_reader.py --inspect data/nisar_l2/GUNW/*/NISAR_L2_*GUNW*.h5
python src/gunw_reader.py --batch --aoi source --auto-ref \
    --export outputs/export_src --csv outputs/gunw_stats_source.csv
```

Then zip `outputs/` and send it back. The .h5 files never move.

**Exports land on a fixed AOI grid**, so anything they send back can be
differenced directly - see below.

---

## Testing

```bash
python -m pytest tests/ -q      # 37 tests
python tests/mutate.py          # 15 mutations, all must be caught
```

A test that passes proves nothing on its own - it may assert something that was
never in danger. `mutate.py` puts each historical bug back and checks the suite
notices. It refuses to run on a dirty working tree, because a crash
mid-mutation would leave a sabotaged reader on disk.

Fixtures are built **backwards from a known displacement field**: the phase is
computed from the answer rather than the answer from the phase, so every step is
checkable against a number we already hold. Every serious fault in this project
produced plausible numbers instead of an error, and none of them could be caught
by looking at the output.

---

## Method notes

Two hazard classes, different physics, different tooling:

- **Class A - precursory creep.** Bare rock, moraine, ice, above treeline.
  InSAR time series works; a real forecast with lead time is possible.
- **Class B - rainfall-triggered.** Vegetated tropical slopes. InSAR early
  warning is *not* possible; use susceptibility plus post-event
  coherence-change mapping.

The source zone is Class A, which is why it was chosen. Note that Class A
failures are **not** necessarily monsoon-driven - Blatten failed on 28 May,
before the monsoon - so a seasonality argument that assumes monsoon-driven creep
is Class B reasoning and does not transfer.

Data spine is NASA: NISAR (L-band), NASADEM and SRTM (topography), GPM IMERG
(rainfall), SMAP (antecedent soil moisture), Landsat and FIRMS (optical,
thermal). Sentinel-1/2 remain as the partner-agency historical archive.

---

## Known limitations

Stated here rather than left for a reader to find.

1. **The final seven days before the failure are unobserved**, and the interval
   covering late August is a 24-day average.
2. **Every L2 time series has zero redundancy.** No internal error estimate is
   possible without a custom pair network from RSLC.
3. **Sensitivity is undetermined below ~10 degrees of slope**, which affects the
   multi-geometry table and every downslope magnitude derived by division.
4. **The stratified troposphere is measured but not removed** from the GUNW
   series - 16% of variance, +/-78 to +/-116 mm over the relief.
5. **The inverse-velocity detector uses one scalar floor** where per-pair floors
   vary twelvefold.
6. **The co-event footprint is not a scar map.** It is broader and less
   terrain-selective than the failure, and peak monsoon is a confound we can
   argue against but not eliminate.
7. **GOFF displacement values are not reproducible between processing chains**
   at this site (median difference +218 mm, scatter 267 mm). Coverage is; use it.
8. **No independent ground validation.** No GNSS, no field survey, no optical
   confirmation of the deformation field.

---

## References

Methods this project uses, rather than a survey.

1. **Voight, B.** (1989). A relation to describe rate-dependent material failure.
   *Science* **243**(4888), 200-203.
   [doi:10.1126/science.243.4888.200](https://doi.org/10.1126/science.243.4888.200)
   - the `dv/dt = A v^alpha` law behind `detectability.py` and the alpha
     discussion in the Blatten section.

2. **Fukuzono, T.** (1985). A method to predict the time of slope failure caused
   by rainfall using the inverse number of velocity of surface displacement.
   *Journal of the Japan Landslide Society* **22**(2), 8-13.
   [J-STAGE](https://www.jstage.jst.go.jp/article/jls1964/22/2/22_2_8/_article)
   - the inverse-velocity construction implemented in `inverse_velocity.py`.

3. **Itoh, K.** (1982). Analysis of the phase unwrapping algorithm.
   *Applied Optics* **21**(14), 2470.
   [doi:10.1364/AO.21.002470](https://doi.org/10.1364/AO.21.002470)
   - unwrapping is unique only where adjacent samples differ by less than pi,
     which is where the lambda/4 ceiling comes from and why it constrains the
     phase *gradient* rather than absolute displacement.

4. **Berardino, P., Fornaro, G., Lanari, R., Sansosti, E.** (2002). A new
   algorithm for surface deformation monitoring based on small baseline
   differential SAR interferograms. *IEEE Transactions on Geoscience and Remote
   Sensing* **40**(11), 2375-2383.
   [doi:10.1109/TGRS.2002.803792](https://doi.org/10.1109/TGRS.2002.803792)
   - the small-baseline temporal inversion `timeseries.py` reduces to, given
     products that are already unwrapped and geocoded.

### Still to cite

Three sources are used as evidence in this README and are not yet properly
referenced. They are load-bearing, so they need real citations before
submission anywhere:

- **The Blatten velocity record** in the calibration table (0.5-0.8 m/day at
  six days, 10 m/day at one day). Currently unsourced. The whole detectability
  argument is calibrated against it.
- **The ESA ALOS-2 / SAOCOM analysis** of Kleines Nesthorn creep rates
  (~50 cm/yr by 2023, >150 cm/yr by August 2024), quoted in "Two regimes, two
  instruments".
- **The NISAR L2 product specification**, for the GUNW and GOFF layer
  definitions, the three-digit mask encoding, and the "raw, unculled,
  unfiltered" characterisation of the offset layers.

---

## Layout

```
src/     nisar_acquisition.py   catalogue search + download
         organise.py            sort products into dated buckets
         gunw_reader.py         GUNW -> LOS displacement
         goff_reader.py         GOFF -> pixel offsets, measured noise floor
         pixel_stack.py         per-pixel velocity, correlation-length null
         timeseries.py          SBAS inversion, network, jackknife
         inverse_velocity.py    Fukuzono forecasting with a significance gate
         geometry_merge.py      LOS -> downslope, sensitivity
         candidate_check.py     four tests against one location
         detectability.py       revisit vs warning time simulation
         impoundment.py         landslide-dam susceptibility
data/    nisar_l2/  dem/        products (gitignored)
outputs/                        GeoTIFFs, quicklooks, CSVs (derived stats kept)
tests/   synth.py mutate.py     fixtures built from a known answer
```
