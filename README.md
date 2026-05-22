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
- **GEDI**: 2020–2023 (instrument lifecycle).

**Justification**: Defaunation is a long-term "press" disturbance. The structural state measured by GEDI in 2020–2023 is the *cumulative result* of the faunal regime that operated over the preceding decades. The FRIP signal measured over 2001–2023 captures the functional consequence of that same regime. By using the JRC TMF **Class 10 (Undisturbed since ~1982)** mask, we guarantee that no major pulse disturbances (logging, fire) have occurred during either window. Therefore, the asynchronous windows are ecologically aligned: GEDI provides the current structural endpoint, and FRIP provides the long-term functional regime that produced it.

---

## Pipeline Architecture

```
01_BaseStack_GEE.ipynb (NB1)
  │  Block 1–4: Exports NPP [24b], CovStack [7b] at MODIS scale;
  │              9× GediStack grids [3b each] at native 25m per basin
  │  Block 5:   Exports GediQuality [2b] masks at 25m per basin
  ▼
02_GEDI_ForestStructure_MODIS_Aggregation_GEE.ipynb (NB2)
  │  Loads 25m GediStack tiles + GediQuality masks
  │  Applies 3-layer quality masking at 25m:
  │    1. GEDI shot quality (l2b_quality_flag == 1, degrade_flag == 0)
  │    2. JRC TMF Class 10 (undisturbed forest)
  │    3. SRTM slope < 10°
  │  Aggregates to MODIS scale, exports GEE Asset + Drive GeoTIFF
  ▼
03_FRIP_Signals_And_Drive_Exports_GEE.ipynb (NB3)
  │  Loads NPP & CovStack (NB1) + GEDI MODIS asset (NB2)
  │  Computes FRIP & Mann-Kendall trend
  │  Exports 42 GeoTIFFs to Drive
  ▼
code/ (R Analysis Pipeline)
  │  00: Synthetic data generation for pipeline testing
  │  01: Load GeoTIFFs + vector masking
  │  02–05: Hypothesis testing (H1–H4)
  │  06: PNAS-style manuscript figures
  │  07: Unit tests for all R functions
```

**Run order**: NB1 → NB2 → NB3 → R Analysis

---

## Shared Spatial Parameters

All GEE notebooks use harmonized constant names:

```python
# Project
GEE_PROJECT  = 'quantum-bonus-434714-t2'
ASSET_ROOT   = f'projects/{GEE_PROJECT}/assets/DefaunationFromSpace'

# Study regions
CONGO_BBOX   = ee.Geometry.Rectangle([8, -12, 35, 8])
AMAZON_BBOX  = ee.Geometry.Rectangle([-73, -18, -44, 8])

# Masking thresholds
FOREST_CLASS_UNDISTURBED = 10  # JRC TMF: undisturbed since ~1982
SLOPE_MAX    = 10              # Degrees: GEDI waveform reliability threshold
FOREST_COVER_THRESHOLD = 0.95  # NB3 only: fraction for FRIP computation

# Analysis scales
SCALES       = list(range(5000, 105000, 5000))  # 5km to 100km
MODIS_SCALE  = 463.3127165279165                 # MODIS equatorial pixel size (m)
```

---

## GEE Design Principles

**This architecture resolves all GEE "User memory limit exceeded" and "Reprojection output too large" errors via:**

1. **Modular Parallel Stacks (NB1)**: The 34-band composite is split into three lightweight, independently exportable asset types per basin.
2. **High-Fidelity GEDI Masking (NB2)**: All quality masks are applied at 25m native resolution *before* spatial aggregation, preventing dilution of signal by contaminated or low-quality pixels.
3. **Idempotent Exports (`safe_start`)**: All asset exports use `safe_start()` which deletes any existing asset before starting, ensuring re-runs never fail with "Asset already exists".
4. **Decoupled Drive Exports (NB3)**: All inputs are pre-materialized static assets, so NB3 has zero memory pressure.
5. **Basin-Specific Geometry Clipping**: Every high-resolution dataset is clipped to the target basin bounding box before operations, bounding the reprojection grid.
6. **GEE-based Mann-Kendall Trend**: Compresses 23 annual FRIP bands into a single `frip_mk_tau` band in-engine, saving 70% export size.
7. **Deferred Vector Masking**: HydroSHEDS basins, WDPA protected areas, and country boundaries are rasterized and masked in R, keeping GEE fully raster-based.

---

## NB1: Base Stack Exports (`01_BaseStack_GEE.ipynb`)

**Compute**: GEE Python Colab | **Exports**: GEE Assets (22 + 2 tasks per run)

Five code blocks: Setup → Functions → Unit Tests [7/7] → Modular Exports → Quality Mask Exports

### Exported Assets

#### 1. `NppStack_{basin}` (24 bands)
- **Bands**: `Npp_median`, `NPP_2001` – `NPP_2023` (MODIS MOD17A3HGF)
- **Resolution**: Native MODIS scale (zero aggregation overhead)

#### 2. `GediStack_{basin}_{gi}` (3 bands × 9 grids per basin)
- **Bands**: `GEDI_UOI`, `GEDI_N`, `GEDI_rh98`
- **Resolution**: Native 25m
- **Grid Partition**: 3×3 spatial grid to circumvent GEE memory limits

#### 3. `CovStack_{basin}` (7 bands)
- **Bands**: `flood_freq`, `forest_fraction`, `elevation`, `slope`, `hnd`, `precip`, `clay`
- **Resolution**: Aggregated from native high-res datasets (7 `reduceResolution` chains)

#### 4. `GediQuality_{basin}` (2 bands)
- **Bands**: `quality_min` (min of `l2b_quality_flag`), `degrade_max` (max of `degrade_flag`)
- **Resolution**: Native 25m
- **Purpose**: Shot-level quality mask applied in NB2 before spatial aggregation. Since 97% of GEDI pixels have exactly 1 shot (N=1), post-average quality filtering is functionally equivalent to pre-average filtering.

---

## NB2: GEDI Forest Structure MODIS Aggregation (`02_GEDI_ForestStructure_MODIS_Aggregation_GEE.ipynb`)

**Compute**: GEE Python Colab | **Loads**: GediStack + GediQuality from NB1 | **Exports**: GEE Assets + Drive (4 tasks)

Four code blocks: Setup → Functions → Unit Tests [5/5] → Export

### Quality Masking (applied at 25m before aggregation)

Three masks are applied in `apply_quality_masks()`, all at native 25m resolution:

1. **GEDI shot quality** — `quality_min == 1 AND degrade_max == 0` (from pre-exported `GediQuality` asset). Removes pixels where any contributing shot had low waveform fidelity or degraded pointing/positioning.
2. **JRC TMF Class 10** — Continuously undisturbed moist forest since ~1982.
3. **SRTM slope < 10°** — Excludes steep terrain where GEDI waveform processing degrades due to footprint-scale elevation spread.

A polarity unit test [5/5] verifies the mask logic using synthetic constant images (keep quality=1/degrade=0, remove quality=0 or degrade>0).

### Output Structure (3 bands at MODIS scale)
- **`uoi`**: Understory Openness Index — spatial `mean`
- **`rh98`**: 98th canopy height percentile — spatial `mean`
- **`gedi_n`**: Footprint shot count — spatial `sum`

---

## NB3: FRIP Signals and Drive Exports (`03_FRIP_Signals_And_Drive_Exports_GEE.ipynb`)

**Compute**: GEE Python Colab | **Loads**: Assets from NB1 + NB2 | **Exports**: Drive GeoTIFFs (42 tasks)

### 1. Multi-scale Stacks (40 total — 20 scales × 2 basins, 11 bands each)
- **`frip`**: Cross-sectional Spearman correlation
- **`frip_mk_tau`**: Mann-Kendall trend τ of annual FRIP (2001–2023)
- **`uoi`**, **`rh98`**, **`gedi_n`**: GEDI structural signals (pre-masked at 25m in NB2)
- **`elevation`**, **`slope`**, **`hnd`**, **`precip`**, **`clay`**, **`forest_fraction`**: Environmental covariates

### 2. Native-scale Stacks (2 total — 1 per basin, 10 bands each)
- GEDI signals + covariates + `Npp_median` at native MODIS resolution (~463m)

**Export path**: `DefaunationSynthesis/AnalysisStack/`

---

## R Analysis Pipeline (`code/`)

**Compute**: Local R | **Input**: `analysis_stack_*.tif` + vector layers

### Orchestration Scripts

| Script | Purpose |
|---|---|
| `00_Generate_Synthetic_Data.R` | Generates spatially realistic synthetic data at full MODIS resolution for pipeline testing |
| `01_Load_And_Join.R` | Load multi-scale and native GeoTIFF stacks, rasterize and join WDPA/country/basin vector layers |
| `02_H1_GEDI_Structural.R` | Regional comparison + protection effect on UOI |
| `03_H2_FRIP_Functional.R` | Same structure for FRIP; temporal defaunation index validation |
| `04_H3_Convergence.R` | PA-scale + pixel-scale UOI–FRIP correlation |
| `05_H4_Temporal_Trends.R` | Analysis of `frip_mk_tau` trends |
| `06_Figures.R` | PNAS-formatted manuscript figures |
| `07_Test_Functions.R` | Unit test suite (60 tests across all function modules) |

### Function Modules (`code/functions/`)

| Module | Responsibility |
|---|---|
| `load_data.R` | GeoTIFF and vector layer loading, rasterization, masking |
| `gedi_analysis.R` | UOI regional comparison, protection effects |
| `frip_analysis.R` | FRIP analysis, defaunation proxy models |
| `convergence_analysis.R` | Structural–functional convergence at PA and pixel scales |
| `temporal_analysis.R` | Mann-Kendall trend analysis |
| `covariate_adjustment.R` | OLS covariate denoising |
| `matched_pairs.R` | Topographic matched-pairs design |
| `pa_pairs.R` | Protected area paired comparisons |
| `plotting.R` | Visualization helpers with downsampling for large datasets |
| `theme_pnas.R` | PNAS journal formatting theme |

---

## Joint Camera Trapping & Biophysical Scaling Analysis

To empirically ground the remote-sensing structural (GEDI) and functional (FRIP) hypotheses, the repository includes a Python-based camera trap community and biophysical scaling pipeline (`code/process_camera_traps.py` and `code/visualise_camera_traps.py`). This pipeline ingests raw Wildlife Insights data packages from the Amazon and Congo basins, collapses image series using a 30-minute independence window, matches species to EltonTraits database masses (`MamFuncDat.txt` & `BirdFuncDat.txt`), corrects relative abundance indices (RAI) for allometric day-range scaling, and estimates landscape cluster-level biophysical metrics (richness $S$, standing Biomass Index $B_H$, and Metabolism Index $M_H$).

### ⚠️ Sampling Effort Bias & Robustness Thresholds
As shown by the diagnostic bias checks (**Figure 4**), camera trap community metrics are highly sensitive to cumulative sampling effort at the cluster scale:
* Taxon richness ($S$) is extremely sensitive to effort ($r_S = 0.733, P < 0.0001$), reflecting standard species-accumulation behaviors.
* Standing biomass ($B_H$) and megafaunal biomass indices ($B_{H, >50}, B_{H, >100}$) show moderate to high positive correlations with cumulative effort ($r_S = 0.40 \text{ to } 0.52$).

> [!IMPORTANT]
> **Robustness Filtering Rule**: For future downstream linkages between camera trap biophysical indices and remote sensing covariates (GEDI/MODIS), **clusters with $\le 100$ cumulative trap-days must be excluded**. Clusters below this threshold are severely under-sampled, resulting in highly biased richness and biomass estimates that confound ecological scaling relationships.

---

## Repository Structure

```
DefaunationSynthesis/
├── README.md
├── 01_BaseStack_GEE.ipynb           # NB1: Base raster exports (GEE → Assets)
├── 02_GEDI_...Aggregation_GEE.ipynb # NB2: GEDI quality masking + MODIS aggregation
├── 03_FRIP_...Exports_GEE.ipynb     # NB3: FRIP computation + Drive exports
│
├── code/
│   ├── 00–07_*.R                    # Analysis orchestration scripts
│   └── functions/                   # Modular R function library (10 modules)
│
├── data/                            # Vector layers, WDPA shapefiles, etc.
├── figures/                         # Generated manuscript figures
├── outputs/                         # Analysis outputs (tables, model summaries)
├── legacy/                          # Reference implementations
│   ├── DefaunationFromSpace_Paper/  # Original FRIP pipeline
│   └── GEDI_openness/              # Original GEDI L2B processing
└── scratch/                         # Notebook modification scripts, one-off utilities
```

---

## Requirements

**NB1, NB2, NB3** (GEE Colabs): `earthengine-api`, `geemap`, `numpy`, `pandas`

**R Analysis**: `terra`, `sf`, `tidyverse`, `ggplot2`, `lme4`, `MuMIn`, `spdep`, `DHARMa`, `patchwork`

---

## Status

| Component | Status |
|---|---|
| `01_BaseStack_GEE.ipynb` | ✅ Complete — all base stacks + quality masks exported |
| `02_GEDI_...Aggregation_GEE.ipynb` | ✅ Complete — 3-layer quality masking + MODIS aggregation |
| `03_FRIP_...Exports_GEE.ipynb` | ✅ Built — awaiting NB2 re-export with quality masks |
| R analysis pipeline (`code/`) | ✅ Built — tested on full-scale synthetic data (60/60 tests pass) |
| Real data analysis | 🔲 Awaiting NB2 + NB3 re-export with quality-filtered GEDI |
| Manuscript figures | 🔲 To produce from real data |
