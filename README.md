# Defaunation From Space — Synthesis Analysis

**Working title**: *Defaunation leaves detectable structural and functional footprints in satellite data*
**Target journal**: Nature Ecology & Evolution

## Scientific Overview

Two independent satellite signals of tropical forest defaunation, Congo and Amazon basins:

| Signal | Metric | Dataset | Causal chain |
|---|---|---|---|
| **Structural** | Understory Openness Index (UOI) | GEDI L2B PAVD | Short: fauna → physical structure |
| **Functional** | Flooding Role in Productivity (FRIP) | MODIS NPP × JRC GLOFAS | Longer: fauna → nutrient pump → productivity pattern |

### Principled Hypotheses
*Based on the theoretical premise that the Amazon is more heavily depleted of large megaherbivores than the Congo (which retains forest elephants).*

- **H1 (Structural)**: Forest understories are more open where megafauna are intact. *Prediction: Congo > Amazon for UOI.*
- **H2 (Functional)**: Flooding more strongly predicts productivity where megafauna are depleted (nutrient pump broken). *Prediction: Amazon > Congo for FRIP.*
- **H3 (Convergence)**: Both signals spatially converge. Protected areas and pixels with high UOI (intact structure) should have low FRIP (intact function).
- **H4 (Temporal)**: The functional signal (FRIP) has strengthened over time where defaunation is increasing.

---

## Temporal Alignment Strategy

The two signals operate on different temporal spans:
- **FRIP**: 2001–2023 (requires >20 years of annual data to compute a robust Spearman rank correlation).
- **GEDI**: 2019–2023 (instrument lifecycle).

**Justification**: Defaunation is a long-term "press" disturbance. The structural state measured by GEDI in 2020 is the *cumulative result* of the faunal regime that operated over the preceding decades. The FRIP signal measured over 2001–2023 captures the functional consequence of that same regime. By using the JRC TMF **Class 10 (Undisturbed since ~1982)** mask, we guarantee that no major pulse disturbances (logging, fire) have occurred during either window. Therefore, the asynchronous windows are ecologically aligned: GEDI provides the current structural endpoint, and FRIP provides the long-term functional regime that produced it.

---

## Pipeline Architecture

```
01_BaseStack_GEE.ipynb          GEE — 34-band base stack at ~463m → GEE Assets
         │
         ▼  (load Assets)
02_Signals_GEE.ipynb            GEE — FRIP at 20 scales + masked GEDI → GEE Assets
         │
         ▼  (load Assets)
03_Build_Analysis_Dataset.ipynb  GEE/local — join signals + labels → GeoTIFFs
         │
         ▼  analysis_stack_{scale}.tif
04_Analysis/  (R scripts)        local R — H1–H4 tests + figures → covariate_models.json
         │
         ▼
05_Denoising_Maps.ipynb          GEE — apply model coefficients → adjusted rasters + maps
```

**Run order**: 1 → 2 → 3 → 4 → 5

### Design Principles

**NB1 + NB2 resolve all "Reprojection output too large" errors via three key strategies:**

1. **Independent `reduceResolution` chains**: Each fine-resolution dataset (25m GEDI, 30m JRC, 30m SRTM, 90m MERIT) is aggregated to MODIS resolution independently. No cross-dataset dependencies — the computation that caused errors (e.g., applying a 30m mask to 25m GEDI within one `reduceResolution`) is eliminated.

2. **No `reproject()` in NB1**: Following [GEE best practices](https://developers.google.com/earth-engine/guides/best_practices), `reproject` is avoided entirely in the base stack export. The export's `scale`/`crs` parameters define the output grid. `reproject` is only used in NB2 where the input is already a ~463m asset and the output is 5–100km (trivially small).

3. **Per-basin exports**: Congo and Amazon exported separately to avoid bounding boxes spanning the Atlantic.

4. **Masking deferred**: `forest_fraction`, `elevation`, `slope` exported as continuous bands. Threshold masks applied only when loading the pre-computed asset in NB2/NB3 — never during a `reduceResolution` chain.

---

## Shared Spatial Parameters

```python
# Bounding boxes — forest mask handles ecological precision
CONGO_BBOX  = ee.Geometry.Rectangle([8,  -12, 35,  8])
AMAZON_BBOX = ee.Geometry.Rectangle([-73, -18, -44, 8])

FOREST_MASK  = 'projects/JRC/TMF/v1_2024/TransitionMap_MainClasses'
FOREST_CLASS = 10   # continuously undisturbed since ~1982
FOREST_COVER_THRESHOLD = 0.95
SCALES       = list(range(5000, 105000, 5000))
GEE_PROJECT  = 'quantum-bonus-434714-t2'
ASSET_ROOT   = 'projects/quantum-bonus-434714-t2/assets/DefaunationFromSpace'
```

---

## NB1: Base Stack Export (`01_BaseStack_GEE.ipynb`)

**Compute**: GEE Python Colab | **Exports**: GEE Assets (2 tasks — one per basin)

Exports a single 34-band raster per basin at MODIS sinusoidal resolution (~463m). Contains ALL variables needed for downstream FRIP and GEDI analysis.

### Base Stack Bands (34 total)

| # | Band | Source | Native res → 463m |
|---|---|---|---|
| 1 | `Npp_median` | MODIS MOD17A3HGF | native |
| 2–24 | `NPP_2001` – `NPP_2023` | MODIS MOD17A3HGF | native |
| 25 | `flood_freq` | JRC GLOFAS v1 (binary ≥0, summed across 7 return periods) + MERIT HND mask | ~4km/90m → 463m |
| 26 | `forest_fraction` | JRC TMF v1_2024 class 10 | 30m → 463m |
| 27 | `elevation` | SRTM | 30m → 463m |
| 28 | `slope` | SRTM-derived | 30m → 463m |
| 29 | `hnd` | MERIT Hydro | 90m → 463m |
| 30 | `GEDI_UOI` | GEDI L2B: 1 − (pavd_z0/pai) | 25m → 463m |
| 31 | `GEDI_N` | GEDI L2B footprint count | 25m → 463m |
| 32 | `GEDI_rh98` | GEDI L2B canopy height | 25m → 463m |
| 33 | `precip` | CHIRPS daily → annual mean | ~5km → 463m |
| 34 | `clay` | OpenLandMap SoilGrids clay fraction | 250m → 463m |

**Assets exported**:
```
projects/.../DefaunationFromSpace/BaseStack_Congo
projects/.../DefaunationFromSpace/BaseStack_Amazon
```

---

## NB2: Signal Exports (`02_Signals_GEE.ipynb`)

**Compute**: GEE Python Colab | **Loads**: Base stack assets from NB1 | **Exports**: GEE Assets

### FRIP Exports (80 tasks)

Cross-sectional + annual FRIP at 20 scales (5–100km) × 2 basins. Applies `forest_fraction >= 0.95` mask before computing Spearman correlation between `flood_freq` and `Npp_median`.

```
projects/.../DefaunationFromSpace/FRIP/FRIP_{scale}_{basin}
projects/.../DefaunationFromSpace/FRIP/FRIP_Annual_{scale}_{basin}
```

### GEDI Exports (2 tasks)

Masked GEDI + covariates at MODIS resolution. Applies `forest_fraction >= 0.95`, `elevation < 1000m`, `slope < 10°`.

Bands: `GEDI_UOI`, `GEDI_N`, `GEDI_rh98`, `elevation`, `slope`, `hnd`, `precip`, `clay`, `forest_fraction`.

```
projects/.../DefaunationFromSpace/GEDI/GEDI_masked_{basin}
```

---

## NB3: Analysis Dataset (`03_Build_Analysis_Dataset.ipynb`)

**Compute**: GEE Python Colab | **Loads**: GEE Assets from NB1 + NB2 | **Exports**: Drive GeoTIFFs

Joins FRIP + GEDI signals with spatial labels (basin, country, protection status). Exports multi-band GeoTIFFs for R.

| Band / Label | Source | Use |
|---|---|---|
| FRIP correlation | NB2 FRIP assets | H2 signal |
| GEDI_UOI | NB2 GEDI assets | H1 signal |
| Covariates | NB1 base stack | Adjustment models |
| protection | `WCMC/WDPA/current/polygons` IUCN I–IV | Protected / unprotected flag |
| basin, sub_basin | `WWF/HydroSHEDS/v1/Basins/hybas_2` | ANOVA blocking |
| country | GAUL boundaries | Country-level grouping |
| mk_tau | Annual FRIP | Mann-Kendall τ |

---

## NB4: Hypothesis Testing (`04_Analysis/`)

**Compute**: Local R scripts | **Input**: `analysis_stack_{scale}.tif` + local DI rasters

| Script | Purpose |
|---|---|
| `Functions.r` | Shared helpers |
| `01_Load_and_Join.r` | Load stacks, join DI rasters locally |
| `02_H1_GEDI_Structural.r` | Regional t-test + protection ANOVA on UOI |
| `03_H2_FRIP_Functional.r` | Same structure on FRIP; DI validation |
| `04_H3_Convergence.r` | PA-scale + pixel-scale UOI–FRIP correlation |
| `05_H4_Temporal_Trends.r` | mk_tau maps, trend distributions |
| `06_Figures.r` | Manuscript figures |

---

## NB5: Covariate Adjustment Maps (`05_Denoising_Maps.ipynb`)

Applies OLS adjustment coefficients from NB4 to produce "denoised" signal maps. Only run if raw signals are validated.

---

## Directory Structure

```
DefaunationSynthesis/
├── README.md
├── 01_BaseStack_GEE.ipynb         ← NEW: 34-band base stack export
├── 02_Signals_GEE.ipynb           ← NEW: FRIP + GEDI signal exports
├── 03_Build_Analysis_Dataset.ipynb
├── 04_Analysis/
│   ├── Functions.r
│   ├── 01_Load_and_Join.r
│   ├── 02_H1_GEDI_Structural.r
│   ├── 03_H2_FRIP_Functional.r
│   ├── 04_H3_Convergence.r
│   ├── 05_H4_Temporal_Trends.r
│   └── 06_Figures.r
├── 05_Denoising_Maps.ipynb
├── covariate_models.json
├── data/
│   ├── processed/
│   └── local/
│       └── DefInd/
├── figures/
└── legacy/
    ├── DefaunationFromSpace_Paper/
    └── GEDI_openness/
```

---

## Requirements

**NB1, NB2, NB3, NB5** (GEE Colabs): `earthengine-api`, `geemap`, `numpy`, `pandas`

**NB4** (R): `terra`, `tidyterra`, `tidyverse`, `ggplot2`, `lme4`, `MuMIn`, `spdep`, `multcompView`

---

## Status

| Component | Status |
|---|---|
| `01_BaseStack_GEE.ipynb` | ✅ Built — unit tests configured |
| `02_Signals_GEE.ipynb` | ✅ Built — unit tests configured |
| `03_Build_Analysis_Dataset.ipynb` | 🔲 To build |
| `04_Analysis/` (R scripts) | 🔲 To build |
| `05_Denoising_Maps.ipynb` | 🔲 To build |
| Legacy FRIP R pipeline | ✅ Reference |
| Legacy GEDI Python pipeline | ✅ Reference |
