# SARGuardian

Multi-hazard Earth intelligence from SAR. Current focus: **glacier and landslide
deformation in the Nepal Himalaya**, built on NASA NISAR L-band interferometry.

> Folder is still named `Weatherly` from an earlier project. Contents are
> SARGuardian only. The previous project is archived at
> `../Weatherly_backup_20260829.zip`.

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

**These are model results.** Absolute warning times depend on assumed precursor
duration, creep amplitude and noise. Calibrate against a documented failure
before quoting any number.

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
