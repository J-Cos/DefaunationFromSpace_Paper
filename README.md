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

## Two-Notebook Pipeline Architecture

The entire GEE analysis data production is decoupled into a robust, two-stage architecture that resolves all computational bottlenecks:

```
01_BaseStack_GEE.ipynb          GEE — 34-band base stack at ~463m → GEE Assets
         │
         ▼  (load Assets)
02_Signals_And_Exports_GEE.ipynb GEE — computes FRIP & GEDI, stacks 33 bands → Drive (GeoTIFFs)
         │
         ▼  analysis_stack_{scale}_{basin}.tif
04_Analysis/  (R scripts)        local R — H1–H4 tests + figures → covariate_models.json
         │
         ▼
05_Denoising_Maps.ipynb          GEE — apply model coefficients → adjusted rasters + maps
```

**Run order**: 1 → 2 → R Analysis → 5

### Design Principles

**This architecture resolves all GEE "Reprojection output too large" and projection limits via three key strategies:**

1. **Independent `reduceResolution` chains (NB1)**: Each fine-resolution dataset (25m GEDI, 30m JRC, 30m SRTM, 90m MERIT) is aggregated to MODIS resolution independently. No cross-dataset dependencies are evaluated during the Stage 1 reduction.
2. **Standard Geographic Projection (`EPSG:4326`)**: Standardizing the Stage 1 base stack exports to standard geographic WGS84 coordinates avoids sinusoidal projection boundary limits (`Can't transform` coordinate error) at the edges of the Amazon basin.
3. **Decoupled Drive Exports (NB2)**: Rather than exporting temporary GEE assets, NB2 computes FRIP (cross-sectional + annual) in memory from the materialized Stage 1 asset, aggregates covariates, and writes the final **33-band stack** directly to Google Drive as a GeoTIFF.
4. **Deferred Vector Masking**: All vector layers (HydroSHEDS basins, WDPA protected areas, countries) are rasterized and masked *locally in R* onto the exported GeoTIFF grids. This keeps GEE fully raster-based and extremely fast.

---

## Shared Spatial Parameters

```python
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

Exports a single 34-band raster per basin at MODIS WGS84 resolution (~463m equivalent). Contains ALL variables needed for downstream FRIP and GEDI analysis.

### Base Stack Bands (34 total)

| # | Band | Source | Native res → 463m |
|---|---|---|---|
| 1 | `Npp_median` | MODIS MOD17A3HGF | native |
| 2–24 | `NPP_2001` – `NPP_2023` | MODIS MOD17A3HGF | native |
| 25 | `flood_freq` | JRC GLOFAS v2_1 (binary depth ≥0, summed across 7 return periods) + MERIT HND mask | ~4km/90m → 463m |
| 26 | `forest_fraction` | JRC TMF v1_2024 class 10 | 30m → 463m |
| 27 | `elevation` | SRTM | 30m → 463m |
| 28 | `slope` | SRTM-derived | 30m → 463m |
| 29 | `hnd` | MERIT Hydro | 90m → 463m |
| 30 | `GEDI_UOI` | GEDI L2B: 1 − (pavd_z0/pai) | 25m → 463m |
| 31 | `GEDI_N` | GEDI L2B footprint count | 25m → 463m |
| 32 | `GEDI_rh98` | GEDI L2A canopy height | 25m → 463m |
| 33 | `precip` | CHIRPS daily → annual mean | ~5km → 463m |
| 34 | `clay` | OpenLandMap SoilGrids clay fraction | 250m → 463m |

**Assets exported**:
```
projects/.../DefaunationFromSpace/BaseStack_Congo
projects/.../DefaunationFromSpace/BaseStack_Amazon
```

---

## NB2: Signals and Drive Exports (`02_Signals_And_Exports_GEE.ipynb`)

**Compute**: GEE Python Colab | **Loads**: Base stack assets from NB1 | **Exports**: Drive GeoTIFFs (40 tasks)

Computes FRIP and GEDI structural indicators in memory, aggregates environmental covariates, and compiles them directly into a unified 33-band GeoTIFF per scale and basin.

### Stacked Bands (33 total)

- **`frip`** (1 band): Cross-sectional Spearman correlation
- **`FRIP_2001` ... `FRIP_2023`** (23 bands): Annual Spearman correlations
- **`uoi`**, **`rh98`**, **`gedi_n`** (3 GEDI bands): Openness, height, footprint count
- **`elevation`**, **`slope`**, **`hnd`**, **`precip`**, **`clay`**, **`forest_fraction`** (6 covariate bands)

**Exported to Drive** (`DefaunationSynthesis/AnalysisStack/`):
```text
analysis_stack_5000_Congo.tif … analysis_stack_100000_Congo.tif
analysis_stack_5000_Amazon.tif … analysis_stack_100000_Amazon.tif
```

---

## NB4: Hypothesis Testing (`04_Analysis/`)

**Compute**: Local R scripts | **Input**: `analysis_stack_{scale}_{basin}.tif` + local DI rasters

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

## Requirements

**NB1, NB2, NB5** (GEE Colabs): `earthengine-api`, `geemap`, `numpy`, `pandas`

**NB4** (R): `terra`, `tidyterra`, `tidyverse`, `ggplot2`, `lme4`, `MuMIn`, `spdep`, `multcompView`

---

## Status

| Component | Status |
|---|---|
| `01_BaseStack_GEE.ipynb` | ✅ Built — running server-side (WGS84) |
| `02_Signals_And_Exports_GEE.ipynb` | ✅ Built — verified, ready to run |
| `04_Analysis/` (R scripts) | 🔲 To build |
| `05_Denoising_Maps.ipynb` | 🔲 To build |
| Legacy FRIP R pipeline | ✅ Reference |
| Legacy GEDI Python pipeline | ✅ Reference |
