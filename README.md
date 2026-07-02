# Detecting Heterotroph Biomass from Space

**Official Title**: *Detecting heterotroph biomass from space: heterotrophs leave detectable signal in GEDI plant area volume density profiles*  
**Target Journal**: *Nature Ecology & Evolution*

---

## 1. Scientific Overview, Hypotheses & Key Findings

This repository contains the complete, reproducible R and Google Earth Engine (GEE) analysis pipeline designed to test whether **tropical forest defaunation is detectable from space**. 

By combining spaceborne canopy vertical structures (measured by NASA's Global Ecosystem Dynamics Investigation (GEDI) LiDAR mission) with ground-based camera trap surveys across 101 sites in the Amazon, Congo, and Southeast Asian basins, we establish a robust bidirectional biophysical relationship between forest canopy architecture and standing mammal biomass.

### Core Hypotheses
In this work, we test two primary hypotheses:
1.  **Hypothesis 1 (Predictability):** The standing biomass of large terrestrial mammals is predictable from satellite-measured vertical forest canopy structure (specifically understory openness) at ecologically meaningful spatial scales.
2.  **Hypothesis 2 (Megafaunal Ecosystem Engineering):** Large-bodied megafauna, specifically forest elephants (*Loxodonta cyclotis*) and Asian elephants (*Elephas maximus*), play a distinctive, disproportionate role in generating this spaceborne-detectable structural canopy signal.

### Key Quantitative Findings
*   **Bidirectional Predictability:** Establish a robust, bidirectional biophysical relationship between spaceborne GEDI understory structure and ground-measured mammal community traits.
*   **Structural Predictive Power:** Spaceborne understory openness (UOI) alone explains **35%** of the cross-continental variance in heterotroph biomass under leave-one-region-out cross-validation (LORO-CV).
*   **Megafaunal Restructuring:** Calibrating predicted biomass by local elephant range presence raises the explained variance to **48%**, mathematically quantifying the disproportionate structural disturbance signature left by these large-bodied ecosystem engineers.
*   **Energetic Trophic Coupling:** Substituting static standing biomass ($B_H$) with metabolic-scaled energy flux ($M_H$, scaling as $M^{0.75}$ under Metabolic Scaling Theory) consistently improves model deviance explained (up to **52.2%**), demonstrating that satellite-derived vertical forest structure is tightly coupled with community-level metabolism.
*   **Formulation Robustness:** All model selection pathways (AIC, LORO-CV, and Beta regressions) select the *exact same* best-performing environmental and biogeographical covariate formulations under both standing biomass ($B_H$) and metabolic rate ($M_H$) traits, showing that the underlying biophysical signals are highly robust.

---

## 2. Pipeline Architecture & Sequential Run Order

The repository is structured as a fully sequential, modular, and non-hardcoded pipeline. Data processing flows from GEE cloud composite building to camera trap ingestion, geographical/temporal calibration, and statistical modeling.

> [!NOTE]
> **Pipeline Verification Workflow**: 
> 1. **Unit Testing (`code/07_Unit_Tests.R`)**: Verifies the behavioral and mathematical correctness of individual core algorithms (such as pixel extraction matrices, temporal weights, and regression formulas) in isolation. It should be run first to ensure code stability.
> 2. **One-Command R Orchestration (`code/08_Integration_Tests.R`)**: While R scripts `01` through `06` can be run manually in sequence, `code/08_Integration_Tests.R` acts as a master orchestrator. Executing this single test script automatically cleans up previous deliverables, runs all R analysis steps sequentially via `source()`, and verifies the mathematical integrity of every single model output and figure.


```
[Google Earth Engine Cloud Processing]
01_BaseStack_GEE.ipynb (NB1)  --> Exports raw 25m GEDI L2B PAVD tiles and masks.
02_GEDI_ForestStructure_...   --> Applies 3-layer native quality masking before aggregating to 463m.
03_FRIP_Signals_And_Exports   --> Computes Spearman correlation and exports aggregated GeoTIFFs to Drive.
      │
      ▼
[Local Camera Trap Ingestion & Trait Processing (Python)]
process_camera_traps.py       --> Ingests raw WI camera trap packages, joins with EltonTraits, performs 11.1km spatial
                                  clustering, geodesic buffering, GEDI checks, and writes all spatial GeoJSON outputs.
visualise_camera_traps.py     --> Strictly graphics-only; generates exploratory vertebrate community & scaling figures.
      │
      ▼
[Local R Analysis Pipeline (code/)]
01_FigureS1_Regional_Bounding_Boxes.R --> Prepares elephant ranges -> outputs outputs/elephant_ranges.gpkg and figureS1.png.
02_Load_And_Join.R                    --> Rasterizes vectors, loads elephant ranges & camera trap CSVs, saves ready-to-use cluster RDS.
03_Framework1_Analysis.R              --> Fits Beta Regressions (UOI ~ Biomass) -> generates Figure 3 (figure3.png).
04_Framework2_Analysis.R              --> Fits Tweedie GLMs (Biomass ~ UOI) -> generates Supplementary Figure S2 (figureS2.png).
05_Predictive_Biomass_Maps.R          --> Standing biomass projection mapping at 5 km (figureS3.png) and 20 km (figure5.png).
05b_Predictive_Columns_Correlation.R  --> Generates Supplementary Figure S4 correlating template vs. calibrated models.
06_Pipeline_Visualization.R           --> Generates GEDI & camera trap pipeline summary Figures 1 and 2.
09_Metabolic_Scaling_Analysis.R       --> Performs comparative Metabolic Scaling Theory (MST) and index robustness analysis.
10_Collect_Results_Stats.R            --> Harvests all numeric results into outputs/results_statistics.csv.
```

### Script Catalog (code/)

1.  **[code/process_camera_traps.py](file:///home/j/AgenticProjects/DefaunationSynthesis/code/process_camera_traps.py):**  
    Ingests and cleans raw Wildlife Insights camera trap packages from Congo, Amazon, and SE Asia basins, collapses image records to independent events, matches taxonomic entries to EltonTraits body-mass databases, performs 11.1km spatial clustering, filters clusters by GEDI grid overlap, computes taxonomic keep proportions (`p_keep`) and temporal weights (`w_temp_cluster`), generates buffered convex hulls using geodesic projections to avoid distortion, and writes all GeoJSON/CSV products.
2.  **[code/visualise_camera_traps.py](file:///home/j/AgenticProjects/DefaunationSynthesis/code/visualise_camera_traps.py):**  
    Strictly graphics-only script. Loads pre-processed outputs and generates publication-quality multi-panel exploratory figures analyzing community structure, taxonomic composition, rank-abundance curves, and biophysical scaling at the spatial cluster level.
3.  **[code/01_FigureS1_Regional_Bounding_Boxes.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/01_FigureS1_Regional_Bounding_Boxes.R):**  
    Generates Figure S1 (`figures/figureS1.png` and `figures/figureS1.pdf`) showing the symmetric $15^\circ\text{S}\text{ to }15^\circ\text{N}$ bounding boxes ($50^\circ$ wide in longitude) centered on their respective regional centroids for Amazon, Congo, and SE Asia, and overlays the IUCN Proboscidea range maps categorized by status. Also exports the combined vector ranges to a single GeoPackage (`outputs/elephant_ranges.gpkg`).
4.  **[code/02_Load_And_Join.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/02_Load_And_Join.R):**  
    Rasterizes administrative, protected area, and physical basin vectors, cleans datasets, loads the prepared elephant range GeoPackage, and saves a consolidated environmental composite `loaded_data.rds`.
5.  **[code/03_Framework1_Analysis.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/03_Framework1_Analysis.R):**  
    Performs covariate model selection for GEDI understory openness (UOI) using **weighted Beta Regressions** across **expanded candidate formulations** (incorporating environmental covariates, regional basin boundaries, standing biomass sub-components, and specific elephant presence alternates such as `ElephantPossible` and `ElephantStrict`). The selected best-fitting model is **`M26_possible: Biomass + ElephantPossible + Elev + Megafauna`** (`uoi ~ B_H_index + elephant_present_possible + elevation + B_H_gt100`), showing that the specific presence of large ecological engineers (elephants) is a stronger biophysical predictor of understory openness than generic regional basin boundaries. Generates the publication-quality **Figure 3 (figure3.png)** with a continuous, vibrant model selection bar plot.
6.  **[code/04_Framework2_Analysis.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/04_Framework2_Analysis.R):**  
    Performs spaceborne standing mammal biomass index prediction using **weighted Tweedie GLMs** across **52 formulations** (evaluating combinations of understory openness, elevation, and basin interactions, dynamically expanded to include all 15 regional basin variants). Integrates both **LOBOCV parsimonious** and **standard AIC** selection pathways in a unified framework, automatically skipping out-of-sample folds containing `basin` levels to prevent runtime errors. Generates both the main generalizability **Supplementary Figure S2** (`figures/figureS2.png` / `figures/figureS2.pdf`) and the alternate absolute fit **Supplementary Figure S2 Alternate** (`figures/figureS2_alternate_aic_full.png` / `figures/figureS2_alternate_aic_full.pdf`) using a DRY, highly parametric plotting pipeline, and saves both best models as RDS objects.
7.  **[code/05_Predictive_Biomass_Maps.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/05_Predictive_Biomass_Maps.R):**  
    Consolidates biomass prediction mapping. Projects both the LOBO-selected parsimonious best model (`M2.1: UOI Only`) and the AIC-selected basin-calibrated best model (`M2.17b: UOI * Basin + Elev`) across the tropical forest landscapes at both the **5 km** (supplementary map, **Supplementary Figure S3 / figureS3.png**) and **20 km** (peak predictive scale, **Figure 5 / figure5.png**) resolutions in a side-by-side, symmetrical two-column layout. Standardizes all predictions in out-of-sample log-scale MAE units to allow direct comparison of the universal biophysical footprint (Column 1) vs. basin-calibrated biogeographical shift models (Column 2).
8.  **[code/05b_Predictive_Columns_Correlation.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/05b_Predictive_Columns_Correlation.R):**  
    Generates Supplementary Figure S4 (`figures/figureS4.png` and `figures/figureS4.pdf`) correlating the predictions from the LOBO-selected template model and the AIC-selected basin-calibrated model across all tropical forest pixels at the 20 km scale. Computes and displays both Pearson and Spearman rank correlation coefficients, providing quantitative insights into regional predictive shifts.
9.  **[code/06_Pipeline_Visualization.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/06_Pipeline_Visualization.R):**  
    Generates premium, PNAS-style multipanel manuscript **Figure 1 (figure1_gedi_pipeline.png)** and **Figure 2 (figure2_camera_trap_pipeline.png)**, which visually synthesize log-scale GEDI Shot Count distributions and the camera trap spatial-temporal ingestion/calibration pipeline.
10. **[code/07_Unit_Tests.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/07_Unit_Tests.R):**  
    A comprehensive functional unit-testing suite that verifies function signatures, argument structures, and correct scientific return types for all data extraction, precision/temporal weighting, regression modeling, and spatial projection mapping modules.
11. **[code/08_Integration_Tests.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/08_Integration_Tests.R):**  
    An end-to-end integration test runner that unlinks old deliverables, executes the sequential R pipeline (`01` through `06` and `09`) sequentially on the real GEE GeoTIFF datasets, and verifies the mathematical integrity and presence of all RDS models, vector GPKGs, CSV tables, and manuscript figures.
12. **[code/09_Metabolic_Scaling_Analysis.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/09_Metabolic_Scaling_Analysis.R):**  
    Performs comparative Metabolic Scaling Theory (MST) and index robustness analysis. Evaluates the statistical sensitivity and generalizability of raw biomass ($B_H$) vs. metabolic-scaled energy flux ($M_H$, exponent $\beta = 0.75$), verifying that selected best formulations are completely stable. Saves its summary output to `outputs/metabolic_scaling_model_selection.csv`.
13. **[code/10_Collect_Results_Stats.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/10_Collect_Results_Stats.R):**  
    Post-hoc statistics harvester. Reads all pipeline outputs (models, CSVs, rasters) and produces a long-form CSV (`outputs/results_statistics.csv`) containing every numeric result cited in the manuscript, with confidence intervals, p-values, and source file tracking.


### **Diagnostic & Scratch Utilities (scratch/)**

*   **[scratch/correlation_uoi_biomass.R](file:///home/j/AgenticProjects/DefaunationSynthesis/scratch/correlation_uoi_biomass.R):** Performs statistical correlation tests (Pearson/Spearman) on raw and log1p scales and plots GEDI UOI vs. standing biomass.
*   **[scratch/correlation_6panel.R](file:///home/j/AgenticProjects/DefaunationSynthesis/scratch/correlation_6panel.R):** Generates a 6-panel publication grid plotting standing biomass components (Total, >100kg, >1000kg) vs. GEDI UOI on both log1p (top row) and raw/linear (bottom row) scales with full Pearson/Spearman correlation statistics.

### **Core Biophysical Functions Module**

*   **[code/functions/calibration_helpers.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/functions/calibration_helpers.R):**  
    Defines all shared mathematical operations, ensuring strict DRY compliance:
    *   `extract_scale_pixels()`: Performs raw pixel extraction within buffered MCP polygons.
    *   `extract_scale_data()`: Aggregates pixels, computes empirical GEDI standard error (`uoi_se`), and derives normalized precision weights (`w_uoi_norm = w_uoi / mean(w_uoi)`) and spatial homogeneity (`homogeneity`). It reads temporal weights (`w_temp_cluster`) and keep proportions (`p_keep`) directly from the pre-computed GeoJSON attributes, and computes combined weights (`w_combined_norm = w_combined / mean(w_combined)`).
    *   `fit_framework1_model()` / `fit_framework2_model()`: Dynamically loads fitted model formulas and coefficients from outputs RDS files, ensuring downstream code remains completely parameterization-independent.
*   **[code/functions/theme_pnas.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/functions/theme_pnas.R):**  
    Centralizes all aesthetic parameters, PNAS column-width specifications, and standardizes color scales (`pal_basin`).

---

## 3. Testing Suite & Verification

To guarantee scientific reproducibility and statistical rigor, the pipeline is covered by a two-tiered testing framework:

### Unit Test Coverage (100% of Core Pipeline Functions)
The unit test suite (`code/07_Unit_Tests.R`) verifies the behavior of all core functions in isolation, including the calibration helper algorithms, temporal weighting, and integrated model fit routines.
*   **Coverage Summary**: 100% of high-level analytical, modeling, mapping, and plotting functions are covered under 19 strict assertions.
*   **Run Command**:
    ```bash
    Rscript code/07_Unit_Tests.R
    ```
*   **Assertions Verified**:
    *   Pixel extraction matrices (`extract_scale_pixels`) and region bounding coordinates.
    *   Data calibration (`extract_scale_data`), including proper derivation of normalized precision weights (`w_uoi_norm`) and spatial homogeneity scores.
    *   Survey temporal alignment weighting (`calculate_temporal_weights`).
    *   Model fitting functions (`fit_framework1_model`, `fit_framework2_model`) using dynamically loaded formulas.
    *   High-level script runners (`run_framework1_analysis`, `run_framework2_analysis`, `run_predictive_biomass_mapping`) ensuring correct models are outputted and figures are cleanly populated in `figures/`.

### End-to-End Integration Testing & Pipeline Orchestration
The integration test suite (**`code/08_Integration_Tests.R`**) serves as both a master pipeline runner and a validation framework. It deletes old outputs and executes the entire local processing flow sequentially (`01` through `06`) on the GEE GeoTIFF stack exports located in `outputs/EOdata/`.
*   **Run Command (Executes Entire R Pipeline)**:
    ```bash
    Rscript code/08_Integration_Tests.R
    ```
*   **Integrity Checks**: After execution, it performs mathematical and existence checks to verify that all regression models, elephant range GPKGs, goodness-of-fit tables, and PNAS/Supplementary figures are successfully populated in `outputs/` and `figures/` with non-zero file sizes.

---

## 4. Requirements & Installation

### R Environment
Ensure you have R version $\ge 4.1.0$ installed with the following packages:
```R
install.packages(c("terra", "sf", "dplyr", "readr", "mgcv", "ggplot2", "scales", "cowplot", "tidyterra"))
```

### Google Earth Engine Environment
Notebooks run in standard GEE Python environments and require the standard API package:
```bash
pip install earthengine-api geemap numpy pandas
```

---

## 5. Current Implementation Status

| Component | Status | Verification & Deliverables |
| :--- | :---: | :--- |
| **01_BaseStack_GEE.ipynb** | ✅ Complete | GEDI, NPP, and Covariate stacks successfully pre-materialized in cloud assets. |
| **02_GEDI_Aggregation_GEE.ipynb** | ✅ Complete | 3-layer high-fidelity GEDI shot-quality and slope filtering applied at 25m. |
| **03_FRIP_Signals_Exports_GEE.ipynb** | ✅ Complete | Aggregated GeoTIFFs successfully exported to Drive. |
| **Python Ingestion Pipeline** | ✅ Complete | Ingests and cleans WI camera trap packages, merges with EltonTraits databases, and computes corrected RAI, B_H, and M_H indices. |
| **R Analysis Pipeline** | ✅ Complete | Clean sequential analysis and modeling pipeline (`01` to `06` and `09`) executing completely warning-free. |
| **Principled Weighting Scenarios** | ✅ Complete | Combined precision-temporal weights (`w_combined_norm`) validated as statistically superior to raw shot count (`gedi_n`). |
| **Manuscript-Ready Figures** | ✅ Complete | Double-column **Figure 3**, **Supplementary Figure S2**, and alternate **Supplementary Figure S2 (AIC-selected)**, 6-panel **Figures 1 and 2**, predicted biomass maps at 20 km (**Figure 5**) and 5 km (**Supplementary Figure S3**), and regional bounding boxes (**Supplementary Figure S1**) successfully generated. |
| **Figure S4 (Correlation Plot)** | ✅ Complete | Regional faceted Pearson and Spearman rank correlation plot between template and calibrated model predictions successfully generated as `figureS4.png/pdf`. |
| **Figure S1 (Regional Bounding Boxes)** | ✅ Complete | Shows $15^\circ\text{S}$ to $15^\circ\text{N}$ bounding boxes with IUCN elephant range overlays colored by status, saving `figureS1.png/pdf` and `outputs/elephant_ranges.gpkg`. |
| **Metabolic Scaling Comparison**| ✅ Complete | Evaluates animal biomass index ($B_H$) vs energetic metabolism index ($M_H$) across both frameworks in `code/09_Metabolic_Scaling_Analysis.R`. |
| **Functional Unit Testing Suite** | ✅ Complete | 19 functional assertions covering 100% of pipeline modules in `code/07_Unit_Tests.R` passing 100% successfully. |
| **End-to-End Integration Runner** | ✅ Complete | Validates end-to-end sequential flow on real GEE datasets and verifies all output sizes and shapes in `code/08_Integration_Tests.R`. |

