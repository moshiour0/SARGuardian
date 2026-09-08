# Reproducing the headline results

Every number in the README comes from one of the commands below. This page is
the path from public data identifiers to those numbers, so a reader can check
them without asking us anything.

Two of the four results need **no data at all** - start there if you only have
ten minutes.

---

## 0. What you need

```bash
git clone https://github.com/moshiour0/SARGuardian.git
cd SARGuardian
python3 -m venv .venv && source .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python -m pytest tests/ -q
python tests/mutate.py
```

`mutate.py` refuses to run on a dirty working tree - that is deliberate, it
edits source in place and must be able to restore it exactly.

Earthdata credentials are needed **only** for the two results that read
products. Register free at <https://urs.earthdata.nasa.gov> and put them in
`~/.netrc` (`~/_netrc` on Windows):

```
machine urs.earthdata.nasa.gov login YOUR_USER password YOUR_PASS
```

---

## Result 1 - revisit near precursor/3, and daily revisit is *not* worse

**Needs: nothing. Runs in about a minute.**

```bash
python src/detectability.py --sweep --precursor 5 10 20 40
```

**Expect** the detection-collapse threshold at `precursor / 2.5` in all four
cells, and - with the default noise-scaled significance gate - a **0.0%
false-alarm rate at every revisit tested**:

| Precursor | Revisit | Detection | False alarm | Prediction error |
|-----------|---------|-----------|-------------|------------------|
| 10 days | 1 day | ~98% | ~0.0% | ~0.7 d |
| 10 days | 3 days | ~79% | ~0.0% | ~0.9 d |
| 40 days | 1 day | ~66% | ~0.0% | ~0.3 d |
| 40 days | 12 days | ~53% | ~0.0% | ~2.7 d |

Monte Carlo, so exact percentages move by a point or two. The **ordering** is
the result: shorter revisit gives a smaller prediction error and a shorter
warning; longer revisit gives more warning, a fuzzier date, and eventually no
detection at all.

### The correction, and how to reproduce the old numbers

This page used to promise the opposite - ~11% false alarms at daily revisit
rising to ~21% - and that came from a detector whose threshold was a hardcoded
1.0 mm/day, unreachable from the command line and independent of noise. At 5 mm
displacement noise the velocity noise is 7.1 mm/day at daily revisit, so the
gate sat an order of magnitude below it.

`--gate fixed` reproduces the old behaviour so the two are comparable:

```bash
python src/detectability.py --sweep --precursor 10 40 --gate fixed
python src/detectability.py --sweep --precursor 10 40 --gate noise   # default
```

**Expect** 12.0% and 22.0% false alarms at 1-day revisit under `fixed`, and
0.0% under `noise`. Both sweeps are committed:
`outputs/detectability.csv` and `outputs/detectability_fixed_gate.csv`, each
carrying its `gate` and `gate_mm_day` columns so a table cannot be quoted
without them.

The Blatten-calibrated pair, now that `--blatten` actually does something:

```bash
python src/detectability.py --sweep --blatten --noise 5 --wavelength 0.2384  # phase
python src/detectability.py --sweep --blatten --noise 300                    # offsets
```

Phase saturates at 100% and never detects. Offset tracking reaches ~100% at 1
and 2 days with no false alarms, and nothing at 4 days or longer.

---

## Result 2 - three of five tracks can see the failure point

**Needs: nothing but a network connection (SRTM via OpenTopoData).**

```bash
python src/geometry_merge.py --sensitivity --lat 28.2877 --lon 85.5281 --stencil-sweep
```

**Expect** slope ~27 deg, aspect ~273 deg (west-facing), elevation ~5166 m, and
**three usable tracks, not two**:

| Track | Look side | Sensitivity | Verdict |
|-------|-----------|-------------|---------|
| S1 DESC 19 / 121 | right | -0.913 | usable |
| NISAR ASC 98 | left | -0.889 | usable |
| S1 ASC 85 | right | +0.188 | blind |
| NISAR DESC 48 | left | +0.163 | blind |

If you get the mirror image of this - both ascending usable, all descending
blind - you are running a version that models Sentinel-1 as left-looking.
NISAR looks left; Sentinel-1 looks right. `pytest tests/test_geometry_merge.py`
asserts it against the look vectors the products themselves carry.

The `--stencil-sweep` output should show the sign **stable** across DEM stencil
widths here. Run it at `--lat 28.27484 --lon 85.47405` instead and it is not -
that is the gentle-terrain limit, and the tool refuses to give a verdict there.

---

## Result 3 - the impoundment grid floor, found at Blatten

**Needs: nothing but a network connection.**

Run **all four at the same heights** — a ratio between two areas is meaningless
otherwise. Blatten against Langtang is 1.69x at 150 m and 1.15x at 10 m.

```bash
for A in source langtang lhende blatten; do
  python src/impoundment.py --api --aoi $A --heights 10 25 50 100 150 \
      --rank-by efficiency --geojson outputs/dam_sites_$A.geojson
done
```

PowerShell:

```powershell
foreach ($A in "source","langtang","lhende","blatten") {
  python src/impoundment.py --api --aoi $A --heights 10 25 50 100 150 `
      --rank-by efficiency --geojson outputs/dam_sites_$A.geojson
}
```

**Expect** sites responding at a 10 m blockage: **6/12** for Blatten, 4/12
source zone, 2/12 Langtang, 1/12 Lhende. That first column is what separates the
Loetschental, not the volumes.

**And expect** a Blatten site at **46.4179 N, 7.8141 E**, about 500 m from the
village, that impounds **nothing at 10 m** and 0.73 Mm3 at 25 m. The observed
lake after the 28 May 2025 Birch Glacier collapse was about 10 m deep, so the
tool identifies the valley and ranks it correctly, then cannot resolve the
blockage at the reach where it happened. That is the grid floor.

---

## Result 4 - no motion above 40.4 mm/day at the failure point

**Needs: 16 NISAR L2 GOFF products, about 17 GB.**

### The products

Search ASF for `NISAR` `GOFF` over `85.4645,28.2453,85.5562,28.3529`, or let the
tool do it:

```bash
python src/nisar_acquisition.py --recon              # what exists, no credentials
python src/nisar_acquisition.py --download GOFF      # needs credentials
python src/organise.py                               # dry run
python src/organise.py --apply
```

The 16 routine (`PR`) granules the result uses, ascending path 098 and
descending path 048:

```
NISAR_L2_PR_GOFF_006_048_D_074_007_2000_SH_20251125T125813...
NISAR_L2_PR_GOFF_006_098_A_016_007_4000_SH_20251128T233919...
NISAR_L2_PR_GOFF_007_048_D_074_008_2000_SH_20251207T125814...
NISAR_L2_PR_GOFF_007_098_A_016_008_4000_SH_20251210T233920...
NISAR_L2_PR_GOFF_008_048_D_074_009_2000_QD_20251219T125815...
NISAR_L2_PR_GOFF_008_098_A_016_009_4000_SH_20251222T233921...
NISAR_L2_PR_GOFF_009_048_D_074_010_2000_QD_20251231T125815...
NISAR_L2_PR_GOFF_009_098_A_016_010_4000_SH_20260103T233921...
NISAR_L2_PR_GOFF_024_048_D_074_025_2000_SH_20260629T125814...
NISAR_L2_PR_GOFF_024_098_A_016_025_4000_SH_20260702T233920...
NISAR_L2_PR_GOFF_025_048_D_074_026_2000_SH_20260711T125814...
NISAR_L2_PR_GOFF_025_098_A_016_026_4000_SH_20260714T233920...
NISAR_L2_PR_GOFF_026_048_D_074_028_4000_SH_20260723T125813...
NISAR_L2_PR_GOFF_026_098_A_016_028_4000_SH_20260726T233919...
NISAR_L2_PR_GOFF_028_048_D_074_029_2000_SH_20260816T125812...
NISAR_L2_PR_GOFF_028_098_A_016_029_4000_SH_20260819T233918...
```

The last two span the 26 August 2026 collapse. Urgent-response (`UR`) copies of
those two also exist and are used only for the processing-chain control.

### The commands

```bash
python src/goff_reader.py --batch --aoi source --layer layer2 \
    --csv outputs/goff_stats_source.csv --export outputs/export_goff_src

python src/timeseries.py --dir data/nisar_l2/GOFF --product GOFF \
    --goff-layer layer2 --aoi source --invert --auto-ref --jackknife \
    --csv outputs/ts_goff_source.csv

python src/inverse_velocity.py --ts outputs/ts_goff_source.csv \
    --floors outputs/goff_stats_source.csv --floors-layer layer2 \
    --noise-floor 18.6 --event-date 2026-08-26
```

`--event-date` is also the **forecast cutoff**. The series it reads contains 28
and 31 August, both after the collapse, and the run must print:

```
FORECAST CUTOFF 2026-08-26: dropped 1 interval(s)
  2026-08-19 -> 2026-08-31  -2.78 mm/day   ends on or after the cutoff
```

with a second drop of `2026-08-16 -> 2026-08-28  +8.74 mm/day` on descending.
If those lines are missing you are running a version that lets the event into
its own forecast. Both intervals sit below the floor here, so the bound is
unchanged either way - but on data where they do not, the detector will
announce a lead time built entirely on hindsight.

### Expected numbers

`outputs/goff_stats_source.csv`, HH/layer2, routine products only:

| Track | n | 3-sigma floor, mm/day |
|-------|---|----------------------|
| ASC 098 | 8 | median 19.8, range 8.9-24.5 |
| DESC 048 | 8 | median 98.3, range 9.4-117.9 |

Ascending summer block, 2 July to 19 August 2026:

| Interval | Days | Velocity |
|----------|------|----------|
| 2026-07-02 -> 2026-07-14 | 12 | -3.39 mm/day |
| 2026-07-14 -> 2026-07-26 | 12 | +1.72 mm/day |
| 2026-07-26 -> 2026-08-19 | 24 | -0.41 mm/day |

Fitted linear velocity **-0.711 mm/day**, and `--jackknife` must report the
block as a chain with redundancy 0 - **untested, not confirmed**. That warning
is part of the result.

### The floor at the failure point, which is the one to quote

Everything above is measured over the whole 82 km2 polygon. Add the target flag
and the answer changes by 1.7x:

```bash
python src/local_floor.py --dir outputs/export_goff_src --match layer2 \
    --lat 28.28771 --lon 85.52809 --sweep 3 6 12 --exclude 20260828 _UR_
```

**Expect** an ASC 098 median floor of **19.8 mm/day over the AOI against 33.4 at
the point**, and the last pre-event pair `20260726_20260819` to go from **8.9 to
34.3 mm/day on 26 of 169 valid pixels** - the pair with the best AOI floor in
the archive is the worst one at the failure point. At radius 3 that pair holds
2 valid pixels and is refused; at radius 12 it recovers, because the window has
pulled in terrain that did not fail.

**The bound to quote is 40.4 mm/day, not 33.4.** 33.4 is the median across all
eight ascending pairs, five of them winter pairs outside the window being
bounded. Over the seven weeks before failure the three covering intervals give
40.4, 19.2 and 34.3 at the point, and a bound that holds across a window is set
by its weakest interval.

Radii are `(2r+1)`-cell windows at 80 m posting: radius 3 is 560 m, radius 6 is
**1.04 km**, radius 12 is 2.00 km.

Then confirm the conclusion is unchanged at the point:

```bash
python src/timeseries.py --dir data/nisar_l2/GOFF --product GOFF \
    --goff-layer layer2 --aoi source --invert --auto-ref \
    --target-lat 28.28771 --target-lon 85.52809 --target-radius 6
```

**Expect** a summer ascending velocity of **-0.14 mm/day** at radius 6 (1 km),
+2.40 at radius 3 (560 m) and -0.06 at radius 12 (2 km) - each far below the
local floor at that window.

The `+/-` figures this page used to quote alongside them are not usable: they
are an OLS slope error on a *cumulative* series, whose residuals are correlated
by construction, and the summer blocks' only degrees of freedom come from a
routine/urgent duplicate of the same two acquisitions. Compare against the
measured floor instead.

### Shortcut, no products needed

The derived statistics are committed, so the floors and the time series can be
checked without downloading anything:

```bash
python src/timeseries.py --from-stats outputs/goff_stats_source.csv \
    --product GOFF --goff-layer layer2 --aoi source
```

---

## Result 5 - the troposphere is a phase problem, not an offset problem

**Needs: the exported rasters from Result 4, plus the GUNW exports.**

```bash
python src/troposphere.py --dir outputs/export_src --aoi source --report-only
python src/troposphere.py --dir outputs/export_goff_src --aoi source     --match layer2 --report-only
```

The first fetch of the DEM takes a few minutes - OpenTopoData allows one
request per second - and is cached under `outputs/dem_cache/` afterwards.
`--match layer2` is required on the GOFF directory: it holds three layers per
pair and mixing them in one summary is meaningless.

**Expect**, over the source zone:

| Product | Pairs | \|r\| median | Variance explained | Sign reversals (within track) |
|---------|-------|-------------|--------------------|-------------------------------|
| GUNW - phase | 15 | 0.46 | 21.1% | 9 of 12, 6 expected by chance |
| GOFF - offsets | 18 | 0.12 | 1.5% | consistent with chance |

Reversals are counted **within one track, in date order, one observation per
acquisition pair**. Counted across the alphabetical file list - which
interleaves ascending and descending, and double-counts routine/urgent pairs -
the phase figure was 6 of 14, which is *below* the 7 that random signs give.
The variance split, not the reversal count, is the evidence here.

Fourteenfold in explained variance. That is the reason the bound in Result 4
survives: it rests on GOFF.

**Do not read the slope column as the answer.** GOFF slopes are the larger of
the two - median 19.1 mm/km against 8.8 for phase - because offset fields are
one to two orders of magnitude noisier, so a steeper line explains less of what
is there.

**And expect ten of fifteen phase pairs flagged `!` for leverage.** The valid
pixels span about 1.1 km of a 4.0 km elevation range, so the fit is
extrapolated across ~3.4x its own IQR. That is why the tool is run with
`--report-only` and the trend is quoted as an error bar rather than subtracted.

### Shortcut, no products needed

Both per-pair fits are committed:

```bash
python -c "import csv,statistics as st; r=[x for x in csv.DictReader(open('outputs/troposphere_goff_source.csv')) if x['usable']=='True']; print(len(r), round(st.median([abs(float(x['r'])) for x in r]),2))"
```

Expect `18 0.12`.

---

## What you cannot reproduce, and why

- **Independent ground truth.** There is none. No GNSS, no field survey, no
  optical confirmation of the deformation field. The bound is what the satellite
  can say, not what the ground did.
- **A tighter bound than the floor.** Any precursor slower than the measured
  floor is invisible to this product. That is the point of quoting the floor.
- **Anything from GUNW over the source zone in the monsoon.** Coverage is 0-1%.
  The reader is not doing it wrong; there is no usable phase there.

---

## If a number does not match

Open an issue with the command, the output, and the product filenames. A result
that does not reproduce is a defect in this repository, not in your run.
