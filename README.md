# Defaunation From Space — Synthesis Analysis

**Working title**: *Defaunation leaves detectable structural and functional footprints in satellite data*
**Target journal**: Nature Ecology & Evolution

## Scientific Overview

Two independent satellite signals of tropical forest defaunation, Congo and Amazon basins:

| Signal | Metric | Dataset | Causal chain |
|---|---|---|---|
| **Structural** | Understory Openness Index (UOI) | GEDI L2B PAVD | Short: fauna → physical structure |
| **Functional** | Flood-Referenced Integrated Productivity (FRIP) | MODIS NPP × JRC GLOFAS | Longer: fauna → nutrient pump → productivity pattern |

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

## Five-Notebook Architecture

```
01_FRIP_GEE.ipynb              GEE — raw FRIP signal → GEE Assets
02_GEDI_GEE.ipynb              GEE — raw GEDI UOI signal → GEE Assets
         │            │
         └─────┬───────┘
               ▼  (load Assets)
03_Build_Analysis_Dataset.ipynb  GEE — predictor stack + join → Drive (GeoTIFFs)
               │
               ▼  analysis_stack_{scale}.tif
04_Analysis/  (R scripts)        local R — H1–H4 tests + figures → covariate_models.json
               │
               ▼
05_Denoising_Maps.ipynb          GEE — apply model coefficients → adjusted rasters + maps
```

**Run order**: 1 → 2 (parallel) → 3 → 4 → 5

### Design Principles
- **NB1 + NB2**: Raw signal production. To apply the 95% forest cover filter correctly, the binary 30m TMF mask is aggregated via `reduceResolution` to the target scale (5–100km), and *then* filtered to `fraction >= 0.95`.
- **NB3**: Data assembly. Computes predictor stack, exports stacked GeoTIFFs. We accept the heavy compute load to generate Mann-Kendall τ at all 20 scales.
- **NB4 (R)**: All hypothesis testing and visualisation, including spatial maps of the raw signals. Computes statistical covariate models.
- **NB5**: Downstream utility. Only if hypothesis tests deem it useful, NB5 applies the covariate coefficients to produce a "denoised" product and its corresponding map.

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
WORKING_SCALE = 25000
PROTECTED_IUCN = ['Ia', 'Ib', 'II', 'III', 'IV']
GEE_PROJECT  = 'quantum-bonus-434714-t2'
ASSET_ROOT   = 'users/JakeWilliams844/DefaunationFromSpace'
```

> **Bounding boxes**: Simpler and transparent. TMF mask excludes non-forest. Sub-basin HYBAS labels added in NB3 via local shapefiles for ANOVA blocking.
> **TMF class 10**: Continuously undisturbed across the full archive — ensures FRIP (2001–2023) and GEDI (2020–2023) operate on stable intact forest.
> **WDPA I–IV**: Maximises protected vs unprotected contrast for hypothesis testing.

---

## NB1: FRIP Signal (`01_FRIP_GEE.ipynb`)

**Compute**: GEE Python Colab | **Exports**: GEE Assets only

**FRIP** = pixel-level Spearman correlation between JRC GLOFAS flood depth (7 return periods summed) and MODIS annual NPP (2001–2023). Forest pixels only (TMF class 10, ≥95% per pixel).

| Dataset | GEE ID |
|---|---|
| MODIS NPP | `MODIS/061/MOD17A3HGF` |
| JRC GLOFAS | `JRC/CEMS_GLOFAS/FloodHazard/v1` |
| JRC TMF | `projects/JRC/TMF/v1_2024/TransitionMap_MainClasses` |
| MERIT Hydro | `MERIT/Hydro/v1_0_1` (HND > 0 mask) |

**Assets exported**:
```
users/JakeWilliams844/DefaunationFromSpace/FRIP_{scale}         (20 — Spearman r)
users/JakeWilliams844/DefaunationFromSpace/FRIP_Annual_{scale}  (20 — 23-band annual)
```

---

## NB2: GEDI Signal (`02_GEDI_GEE.ipynb`)

**Compute**: GEE Python Colab | **Exports**: GEE Assets only

**UOI** = `1 − (pavd_z0 / pai)`. Quality filters: TMF class 10, elevation < 1000m, slope < 10°. No covariate values exported — all labelling in NB3.

| Dataset | GEE ID | Role |
|---|---|---|
| GEDI L2B | `LARSE/GEDI/GEDI02_B_002_MONTHLY` | Signal |
| JRC TMF | `projects/JRC/TMF/v1_2024/TransitionMap_MainClasses` | Forest filter |
| SRTM | `USGS/SRTMGL1_003` | Quality filter only |

**Assets exported** (two-band rasters: UOI + footprint count N):
```
users/JakeWilliams844/DefaunationFromSpace/GEDI_{scale}   (20 assets, bands: UOI_mean, N)
```

---

## NB3: Analysis Dataset (`03_Build_Analysis_Dataset.ipynb`)

**Compute**: GEE Python Colab | **Loads**: GEE Assets from NB1 + NB2 | **Exports**: Drive GeoTIFFs

Computes predictor stack at matching resolutions, joins all signals pixel-by-pixel, pre-computes Mann-Kendall τ per pixel. Exports stacked multi-band GeoTIFFs — R/terra loads these directly.

**Predictor stack + labels (all computed in GEE)**:

| Band / Label | Source | Use |
|---|---|---|
| ndvi_mean, ndvi_var | MODIS MOD09A1 | FRIP covariate adjustment |
| hand_mean, hand_var | MERIT Hydro | FRIP covariate adjustment |
| elev_mean, slope_mean | SRTM | GEDI covariate adjustment |
| flood_prob | JRC GLOFAS | Additional covariate |
| forest_cover | JRC TMF | Coverage QC |
| protection | `WCMC/WDPA/current/polygons` IUCN I–IV | Protected / unprotected flag |
| basin, sub_basin | `WWF/HydroSHEDS/v1/Basins/hybas_2` | Congo / Amazon + PFAF_ID for ANOVA |
| country | GAUL boundaries | Country-level grouping |
| PA_ID | WDPA (focal PAs only) | Numeric ID (with separate CSV lookup) |
| mk_tau | Annual FRIP Asset | Mann-Kendall τ (computed for all 20 scales) |

**Exported to Drive** (`DefaunationSynthesis/AnalysisStack/`):
```text
analysis_stack_5000.tif … analysis_stack_100000.tif
  Multi-band raster for spatial mapping and statistical modelling. 
  All bands above + frip + uoi.

pa_id_lookup.csv
  Mapping of numeric PA_ID to string PA_name.
```

> **Why GeoTIFF?** R/terra handles multi-band rasters natively. Spatial structure is preserved for map figures (H3/H4) and spatial autocorrelation tests, while tabular stats (H1/H2) can be run quickly by casting the rasters to data.frames (`terra::as.data.frame()`).

---

## NB4: Hypothesis Testing (`04_Analysis/`)

**Compute**: Local R scripts | **Input**: `analysis_stack_{scale}.tif` + local DI rasters

R is chosen for its statistical depth (lme4, spdep, MuMIn) and publication-quality spatial visualisation (terra, tidyterra, ggplot2) — matching the existing FRIP R pipeline.

### Scripts

| Script | Purpose |
|---|---|
| `Functions.r` | Shared helpers (label_df, CI extraction, Tukey labelling) |
| `01_Load_and_Join.r` | Load stacks, join DI rasters (Benítez-López + Bogoni) locally |
| `02_H1_GEDI_Structural.r` | Regional t-test + protection ANOVA on raw UOI; multi-scale CI |
| `03_H2_FRIP_Functional.r` | Same structure on raw FRIP; DI validation models |
| `04_H3_Convergence.r` | PA-scale + pixel-scale Spearman r between UOI and FRIP |
| `05_H4_Temporal_Trends.r` | mk_tau maps, PA time series, trend distribution |
| `06_Figures.r` | All manuscript figures (Fig 2–5) |

### Model output
Scripts export `covariate_models.json` — OLS coefficients from the covariate adjustment models (fitted as part of H1/H2 testing). This JSON contains dictionaries for **all 20 scales**. These coefficients are consumed by NB5 **only if** the metrics are confirmed useful.

> Models in NB4 are **statistical testing models**, not denoising models. They test whether covariates (NDVI, HAND, elevation, slope) explain variation in the signals. The coefficients happen to be reusable for downstream raster-level adjustment in NB5.

---

## NB5: Covariate Adjustment Maps (`05_Denoising_Maps.ipynb`)

**Compute**: GEE Python Colab | **Inputs**: Raw signal Assets + `covariate_models.json`

If the raw signals are validated in NB4, NB5 applies the multi-scale OLS adjustment coefficients as raster band math. Exports adjusted rasters and a final "denoised" spatial map product.

**Outputs → Drive** (`DefaunationSynthesis/AdjustedSignals/`):
```
frip_adjusted_25000.tif
uoi_adjusted_25000.tif
figures/fig_denoised_product.png
```

---

## Directory Structure

```
DefaunationSynthesis/
├── README.md
├── 01_FRIP_GEE.ipynb
├── 02_GEDI_GEE.ipynb
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
├── covariate_models.json          ← produced by NB4, consumed by NB5
├── data/
│   ├── processed/                 ← Drive outputs (analysis_stack_*.tif, adjusted_*.tif)
│   └── local/
│       └── DefInd/                ← Benítéz-López + Bogoni DI rasters (not in GEE)
├── figures/
│   ├── statistical/               ← from NB4 R scripts
│   └── spatial/                   ← from NB5
└── legacy/
    ├── DefaunationFromSpace_Paper/ ← original FRIP R pipeline (reference)
    └── GEDI_openness/              ← original GEDI Python pipeline (reference)
```

---

## Requirements

**NB1, NB2, NB3, NB5** (GEE Colabs): `earthengine-api`, `geemap`, `numpy`, `pandas`

**NB4** (R): `terra`, `tidyterra`, `tidyverse`, `ggplot2`, `lme4`, `MuMIn`, `spdep`, `multcompView`

---

## Status

| Component | Status |
|---|---|
| `01_FRIP_GEE.ipynb` | 🔲 To build |
| `02_GEDI_GEE.ipynb` | 🔲 To build |
| `03_Build_Analysis_Dataset.ipynb` | 🔲 To build |
| `04_Analysis/` (R scripts) | 🔲 To build |
| `05_Denoising_Maps.ipynb` | 🔲 To build |
| Legacy FRIP R pipeline | ✅ Reference |
| Legacy GEDI Python pipeline | ✅ Reference |
