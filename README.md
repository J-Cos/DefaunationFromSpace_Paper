# Defaunation From Space — Synthesis Analysis

**Working title**: *Defaunation leaves detectable structural and functional footprints in satellite data*  
**Target journal**: *PNAS / Nature Ecology & Evolution*

---

## 1. Scientific Overview & Principled Hypotheses

This repository contains the complete R and Google Earth Engine (GEE) analysis pipeline designed to test whether **tropical forest defaunation is detectable from space**. 

We evaluate two independent satellite signals of defaunation across the Congo and Amazon basins, grounded in the theoretical premise that the Amazon is more heavily depleted of large megaherbivores (such as forest elephants) than the Congo:

| Signal | Metric | Dataset | Biophysical Mechanism |
| :--- | :--- | :--- | :--- |
| **Structural** | Understory Openness Index (UOI) | GEDI L2B PAVD | Direct: physical disturbance, trampling, and understory browsing open the canopy understory (Congo > Amazon UOI sensitivity). |
| **Functional** | Flooding Role in Productivity (FRIP) | MODIS NPP × JRC GLOFAS | Indirect: animal nutrient pump redistribution breaks, making productivity highly dependent on seasonal flood dynamics. |

### Principled Core Hypotheses
*   **H1 (Structural Canopy Response):** Understories are significantly more open (higher UOI) where mammal standing biomass ($B_H$) is intact.
*   **H2 (Functional Ecosystem Response):** Flooding more strongly predicts net primary productivity (FRIP) where mammal populations are depleted (the animal nutrient pump is broken).
*   **H3 (Canopy Convergence):** Intact forest patches (high UOI) converge with intact ecosystem function (low FRIP).
*   **H4 (Decadal Trends):** The functional defaunation signal (FRIP) has systematically strengthened over time where defaunation has intensified.

---

## 2. The Multi-Scale GEDI Shot-Noise Averaging Law

A central finding of this synthesis is the scale dependency of spaceborne wildlife detection:
1.  **Native 500m Sparsity:** GEDI is a spaceborne orbital track lidar. At its raw $500\text{ m}$ native scale, track-level sampling sparsity introduces severe shot noise. While GEDI UOI and mammal biomass are biophysically coupled, standard errors at 500m are inflated, rendering local canopy relationships statistically non-significant.
2.  **Noise-Averaging Aggregation ($5-15\text{ km}$):** Aggregating GEDI UOI to intermediate grid cells averages out fine-scale GEDI orbital shot noise and camera trap positioning variance. This noise reduction reveals highly significant biophysical relationships, matching the biological home-range scales at which mammal engineering operates.
3.  **Spatial Dilution ($>20\text{ km}$):** Beyond $20\text{ km}$, the biotic footprint is diluted as regional environmental gradients (precipitation, soil texture, elevation) dominate canopy variance.

---

## 3. Pipeline Architecture & Sequential Run Order

The repository is structured as a fully sequential, modular, and non-hardcoded pipeline. Data processing flows from GEE cloud composite building to camera trap ingestion and statistical modeling:

```
[Google Earth Engine Cloud Processing]
01_BaseStack_GEE.ipynb (NB1)  --> Exports raw 25m GEDI L2B PAVD tiles and masks.
02_GEDI_ForestStructure_...   --> Applies 3-layer native quality masking before aggregating to 463m.
03_FRIP_Signals_And_Exports   --> Computes Spearman correlation and exports multi-scale GeoTIFFs to Drive.
      │
      ▼
[Local R Analysis Pipeline (code/)]
00_Generate_Synthetic_Data.R  --> Generates mock EO stacks & camera trap CSVs for local testing.
01_Load_And_Join.R            --> Rasterizes vectors, cleans data, and saves ready-to-use cluster RDS.
02_Framework1_Analysis.R      --> Fits weighted Beta Regressions (UOI ~ Biomass) & generates Figure 2.
03_Framework2_Analysis.R      --> Fits weighted Tweedie GLMs (Biomass ~ UOI) & generates Figure 3.
04_Predictive_Biomass_Maps.R  --> Functional, DRY projection mapping at 5 km and 20 km scales.
```

### **R Script Catalog (code/)**

1.  **[code/00_Generate_Synthetic_Data.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/00_Generate_Synthetic_Data.R):**  
    Generates synthetic multi-scale Earth Observation stacks (`analysis_stack_{scale}_{region}.tif`) and mock camera trap detections to enable robust local pipeline verification in the absence of GEE asset connections.
2.  **[code/01_Load_And_Join.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/01_Load_And_Join.R):**  
    Rasterizes administrative, protected area, and physical basin vectors, cleans datasets, and saves a consolidated multi-scale environmental composite `loaded_data.rds`.
3.  **[code/02_Framework1_Analysis.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/02_Framework1_Analysis.R):**  
    Performs covariate model selection for GEDI understory openness (UOI) using **weighted Beta Regressions** across 10 candidate models, runs residual temporal stability checks, and generates the publication-quality Figure 2.
4.  **[code/03_Framework2_Analysis.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/03_Framework2_Analysis.R):**  
    Performs spaceborne standing mammal biomass index prediction using **weighted Tweedie GLMs** across 10 formulations, runs residual checks, and generates the publication-quality Figure 3.
5.  **[code/04_Predictive_Biomass_Maps.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/04_Predictive_Biomass_Maps.R):**  
    Consolidates biomass prediction mapping. Projects the best-fitting Tweedie GLM from Framework 2 across Congo and Amazon landscapes at both the **$5\text{ km}$** (core analysis) and **$20\text{ km}$** (peak predictive) scales.

### **Core Biophysical Functions Module**

*   **[code/functions/calibration_helpers.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/functions/calibration_helpers.R):**  
    Defines all shared mathematical operations, ensuring strict DRY compliance:
    *   `extract_scale_pixels()`: Performs raw pixel extraction within buffered MCP polygons.
    *   `extract_scale_data()`: Aggregates pixels, computes empirical GEDI standard error (`uoi_se`), and derives normalized precision weights (`w_uoi_norm = w_uoi / mean(w_uoi)`) and spatial homogeneity (`homogeneity`).
    *   `calculate_temporal_weights()`: Runs geographical single-linkage clustering (11.1km threshold) on camera trap coordinates and computes deployment temporal alignment weights (`w_temp_cluster`) to account for GEDI temporal offset.
    *   `fit_framework1_model()` / `fit_framework2_model()`: Dynamically loads fitted model formulas and coefficients from outputs RDS files, ensuring downstream code remains completely parameterization-independent.
*   **[code/functions/theme_pnas.R](file:///home/j/AgenticProjects/DefaunationSynthesis/code/functions/theme_pnas.R):**  
    Centralizes all aesthetic parameters, PNAS column-width specifications, and standardizes color scales (`pal_basin`).

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
| **03_FRIP_Signals_Exports_GEE.ipynb** | ✅ Complete | Multi-scale GeoTIFFs successfully exported to Drive. |
| **R Sequential Analysis Pipeline** | ✅ Complete | Clean sequential pipeline (`00` to `04`) executing completely warning-free. |
| **Principled Weighting Scenarios** | ✅ Complete | Spatial homogeneity weighting (`uoi_se`) validated as statistically superior to raw shot count (`gedi_n`). |
| **Manuscript-Ready Figures** | ✅ Complete | Double-column Figure 2 and Figure 3, and single-column predicted maps (5km & 20km) successfully generated. |

---
*Defaunation synthesis modeling completed successfully. All outputs are fully reproducible and verified.*
