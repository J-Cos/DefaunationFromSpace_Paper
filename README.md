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

## Three-Notebook Pipeline Architecture

The GEE analysis data production is decoupled into a robust, sequential three-notebook architecture that resolves all computational bottlenecks:

```
01_BaseStack_GEE.ipynb (Notebook 1)
  │  (Exports NPP [24b], CovStack [7b] at MODIS scale; 9x GediStack grids [3b each] at native 25m per basin)
  ▼
02_GEDI_ForestStructure_MODIS_Aggregation_GEE.ipynb (Notebook 2)
  │  (Loads 25m gridded GEDI stacks, applies JRC undisturbed forest mask at 25m, aggregates to MODIS scale, exports GEE Asset)
  ▼
03_FRIP_Signals_And_Drive_Exports_GEE.ipynb (Notebook 3)
     (Loads NPP & CovStack from NB1 + MODIS GEDI Asset from NB2, computes FRIP & MK-trend, exports 42 GeoTIFFs to Drive)
```

**Run order**: 1 ➔ 2 ➔ 3 ➔ R Analysis ➔ 5

### Design Principles

**This architecture resolves all GEE "User memory limit exceeded" and "Reprojection output too large" errors via four key strategies:**

1. **Modular Parallel Base Stacks (NB1)**: To bypass GEE's memory ceiling, we split the 34-band composite into **three lightweight modular assets** per basin exported in parallel.
   * *NppStack (24 bands)*: Already at MODIS scale; zero `reduceResolution` memory overhead.
   * *GediStack (3 bands, 9 spatial grids)*: Stacks the GEDI L2B Understory Openness Index (UOI), footprint count (N), and L2A height (rh98) into a unified 3-band asset exported at native 25m resolution. Stacking these bands and gridding them into a 3x3 layout bypasses Earth Engine's spatial boundaries and memory limits, completely eliminating server-side "User memory limit exceeded" errors.
   * *CovStack (7 bands)*: Only 7 active `reduceResolution` chains from static datasets.
2. **High-Fidelity GEDI MODIS Aggregation (NB2)**: Mosaics the gridded 25m GediStack assets, applies JRC undisturbed forest mask at native 25m resolution, aggregates to MODIS scale, and exports to a static GEE asset. This isolates the expensive high-resolution masking and aggregation, avoiding GEE memory limits.
3. **Decoupled Drive Exports (NB3)**: NB3 loads the NPP and Covariate stacks from NB1, and the pre-aggregated GEDI MODIS-scale asset from NB2, concatenates them instantly (`ee.Image.cat`), and writes the final **42 GeoTIFFs** directly to Google Drive. Since all inputs are pre-materialized, there is zero memory pressure in Notebook 3.
4. **Basin-Specific Geometry Clipping**: Every high-resolution dataset (SRTM, GLOFAS, JRC TMF, MERIT Hydro, SoilGrids, GEDI) is explicitly clipped to the target basin bounding box (`basin_geom`) *before* executing operations. This bounds the reprojection grid and keeps the pixel grid well below the ~300M pixel limit.
5. **Standard Geographic Projection (`EPSG:4326`)**: Standardizing all GEE exports to standard geographic WGS84 coordinates avoids sinusoidal projection boundary limits (`Can't transform` coordinate error) at the edges of the Amazon basin.
6. **GEE-based Mann-Kendall Trend (mk_tau)**: We compute the Mann-Kendall trend τ of annual FRIP across 2001–2023 directly in Earth Engine. This compresses 23 annual bands into a single highly informative `frip_mk_tau` band, saving 70% in export file size.
7. **Deferred Vector Masking**: All vector layers (HydroSHEDS basins, WDPA protected areas, countries) are rasterized and masked *locally in R* onto the exported GeoTIFF grids. This keeps GEE fully raster-based and extremely fast.

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

## NB1: Base Stack Exports (`01_BaseStack_GEE.ipynb`)

**Compute**: GEE Python Colab | **Exports**: GEE Assets (22 parallel tasks — 11 per basin)

Exports three modular assets per basin (with the high-resolution GEDI stack spatially gridded) to bypass Earth Engine's memory ceiling:

### 1. `NppStack_{basin}` (24 bands)
*   **Bands**: `Npp_median`, `NPP_2001` – `NPP_2023` (MODIS MOD17A3HGF)
*   **Resolution**: Native MODIS scale (0 aggregation overhead).

### 2. `GediStack_{basin}_{gi}` (3 bands, 9 grids per basin)
*   **Bands**: `GEDI_UOI` (Understory Openness Index), `GEDI_N` (footprint count), `GEDI_rh98` (height)
*   **Resolution**: **Native 25m resolution**
*   **Grid Partition**: Split into a 3x3 spatial grid (`gi` from 0 to 8) to circumvent Earth Engine size limits for exporting high-density native 25m layers.
*   **Projection**: Default projection set to `EPSG:4326` at 25m to avoid cross-zone UTM reprojection overhead.

### 3. `CovStack_{basin}` (7 bands)
*   **Bands**: `flood_freq` (GLOFAS + HND mask), `forest_fraction` (JRC TMF), `elevation` (SRTM), `slope` (SRTM slope), `hnd` (MERIT Hydro), `precip` (CHIRPS), `clay` (SoilGrids)
*   **Resolution**: Aggregated from native high-res datasets (only 7 `reduceResolution` memory chains).

---

## NB2: GEDI Forest Structure MODIS Aggregation (`02_GEDI_ForestStructure_MODIS_Aggregation_GEE.ipynb`)

**Compute**: GEE Python Colab | **Loads**: `GediStack` grids from NB1 | **Exports**: GEE Assets (2 tasks)

Loads the 9 gridded native 25m GediStack assets, mosaics them seamlessly, applies high-fidelity JRC TMF undisturbed forest cover masking directly at 25m, and aggregates the masked structural indicators to a single 3-band MODIS-scale image (`uoi`, `rh98`, `gedi_n`) exported as a static asset to GEE Assets.

### Output Structure (3 bands):
*   **`uoi`** (1 band): Understory Openness Index aggregated via spatial `mean`
*   **`rh98`** (1 band): 98th canopy height percentile aggregated via spatial `mean`
*   **`gedi_n`** (1 band): footprint shot count aggregated via spatial `sum`

---

## NB3: FRIP Signals and Drive Exports (`03_FRIP_Signals_And_Drive_Exports_GEE.ipynb`)

**Compute**: GEE Python Colab | **Loads**: Modular assets from NB1 & GEDI MODIS from NB2 | **Exports**: Drive GeoTIFFs (42 tasks)

Loads `NppStack`, `CovStack`, and `GediUndisturbedModis` assets, concatenates them instantly (`ee.Image.cat([npp, covs, gedi])`), computes FRIP and GEDI structural indicators in memory, aggregates environmental covariates, and compiles them directly into 42 unified GeoTIFFs per scale and basin.

### 1. Multi-scale Stacks (40 total — 20 scales × 2 basins)
Contains **11 bands** for multi-scale analysis:
*   **`frip`** (1 band): Cross-sectional Spearman correlation
*   **`frip_mk_tau`** (1 band): Mann-Kendall trend τ of annual FRIP across 2001–2023
*   **`uoi`**, **`rh98`**, **`gedi_n`** (3 GEDI bands): Openness, height, footprint count (pre-masked with undisturbed forest at 25m scale in Notebook 2)
*   **`elevation`**, **`slope`**, **`hnd`**, **`precip`**, **`clay`**, **`forest_fraction`** (6 covariate bands)

### 2. Native-scale Stacks (2 total — 1 per basin)
Contains **10 bands** at native MODIS resolution (~463m) for high-resolution spatial modeling *without* FRIP:
*   **`uoi`**, **`rh98`**, **`gedi_n`** (3 GEDI bands)
*   **`elevation`**, **`slope`**, **`hnd`**, **`precip`**, **`clay`**, **`forest_fraction`** (6 covariate bands)
*   **`Npp_median`** (1 NPP productivity band)

**Exported to Drive** (`DefaunationSynthesis/AnalysisStack/`):
```text
analysis_stack_5000_Congo.tif … analysis_stack_100000_Congo.tif
analysis_stack_native_Congo.tif
analysis_stack_5000_Amazon.tif … analysis_stack_100000_Amazon.tif
analysis_stack_native_Amazon.tif
```

---

## NB4: Hypothesis Testing (`04_Analysis/`)

**Compute**: Local R scripts | **Input**: `analysis_stack_*.tif` + local vector layers

| Script | Purpose |
|---|---|
| `Functions.r` | Shared helpers (including `rasteriseAndMask` cover >= 0.99) |
| `01_Load_and_Join.r` | Load multi-scale and native stacks, join vector layers locally |
| `02_H1_GEDI_Structural.r` | Regional t-test + protection ANOVA on UOI (both native & multi-scale) |
| `03_H2_FRIP_Functional.r` | Same structure on FRIP; temporal DI validation |
| `04_H3_Convergence.r` | PA-scale + pixel-scale UOI–FRIP correlation |
| `05_H4_Temporal_Trends.r` | Analysis of GEE pre-computed `frip_mk_tau` trends |
| `06_Figures.r` | Manuscript figures |

---

## NB5: Covariate Adjustment Maps (`05_Denoising_Maps.ipynb`)

Applies OLS adjustment coefficients from NB4 to produce "denoised" signal maps. Only run if raw signals are validated.

---

## Requirements

**NB1, NB2, NB3, NB5** (GEE Colabs): `earthengine-api`, `geemap`, `numpy`, `pandas`

**NB4** (R): `terra`, `tidyterra`, `tidyverse`, `ggplot2`, `lme4`, `MuMIn`, `spdep`, `multcompView`

---

## Status

| Component | Status |
|---|---|
| `01_BaseStack_GEE.ipynb` | ✅ Built — modular stacks ready to run |
| `02_GEDI_ForestStructure_MODIS_Aggregation_GEE.ipynb` | ✅ Built — standalone GEDI aggregator |
| `03_FRIP_Signals_And_Drive_Exports_GEE.ipynb` | ✅ Built — FRIP computing and Drive exporter |
| `04_Analysis/` (R scripts) | 🔲 To build |
| `05_Denoising_Maps.ipynb` | 🔲 To build |
| Legacy FRIP R pipeline | ✅ Reference |
| Legacy GEDI Python pipeline | ✅ Reference |
