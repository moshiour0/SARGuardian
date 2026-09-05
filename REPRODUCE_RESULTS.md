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

## Result 1 - revisit near precursor/3, not daily

**Needs: nothing. Runs in about a minute.**

```bash
python src/detectability.py --sweep --precursor 5 10 20 40
```

**Expect** the detection-collapse threshold at `precursor / 2.5` in all four
cells, and the false-alarm column to rise as revisit shortens:

| Precursor | Revisit | Detection | False alarm |
|-----------|---------|-----------|-------------|
| 10 days | 1 day | ~91% | ~11% |
| 10 days | 3 days | ~95% | ~0.4% |
| 40 days | 1 day | ~92% | ~21% |
| 40 days | 6 days | ~96% | ~0% |

Monte Carlo, so exact percentages move by a point or two between runs. The
**ordering** is the result and it is stable: daily revisit has the highest
false-alarm rate and the worst prediction error in every case tested.

The Blatten-calibrated pair:

```bash
python src/detectability.py --sweep --precursor 7 --creep 27000 \
    --noise 5 --wavelength 0.2384 --revisit 1 2 4 6 12     # phase
python src/detectability.py --sweep --precursor 7 --creep 27000 \
    --noise 300 --revisit 1 2 4 6 12                       # offset tracking
```

Phase saturates at 100% and never detects. Offsets reach ~89% at 1 day and ~94%
at 2 days, with false-alarm rates of ~14% and ~7%.

---

## Result 2 - four of five tracks are blind at the failure point

**Needs: nothing but a network connection (SRTM via OpenTopoData).**

```bash
python src/geometry_merge.py --sensitivity --lat 28.2877 --lon 85.5281 --stencil-sweep
```

**Expect** slope ~27 deg, aspect ~273 deg (west-facing), elevation ~5166 m, and:

| Track | Sensitivity | Verdict |
|-------|-------------|---------|
| S1 ASC 85 | -0.904 | usable |
| NISAR ASC 98 | -0.890 | usable |
| S1 DESC 19 / 121 | +0.197 | blind |
| NISAR DESC 48 | +0.163 | blind |

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

## Result 4 - no motion above 18.6 mm/day before the failure

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

### Shortcut, no products needed

The derived statistics are committed, so the floors and the time series can be
checked without downloading anything:

```bash
python src/timeseries.py --from-stats outputs/goff_stats_source.csv \
    --product GOFF --goff-layer layer2 --aoi source
```

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
