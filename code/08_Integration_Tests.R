# =============================================================================
# code/08_Integration_Tests.R
#
# Rigorous end-to-end integration testing suite. Executes the entire defaunation
# analysis pipeline sequentially (excluding GEE):
#   Step 1: Elephant Range Mapping (01_FigureS1_Regional_Bounding_Boxes.R)
#   Step 2: Vector loading and joining (02_Load_And_Join.R)
#   Step 3: Framework 1 analysis and model selection (03_Framework1_Analysis.R)
#   Step 4: Framework 2 analysis and model selection (04_Framework2_Analysis.R)
#   Step 5: Predictive biomass mapping (05_Predictive_Biomass_Maps.R)
#   Step 6: Pipeline Visualization (06_Pipeline_Visualization.R)
#
# Then checks the existence and integrity of all output RDS, CSV, and PNG
# figures in figures/ and outputs/ directories.
# =============================================================================

cat("============================================================\n")
cat("=== Starting End-to-End Pipeline Integration Test ===\n")
cat("============================================================\n\n")

# --- Helper function to verify file existence and non-zero size ---
verify_file <- function(file_path, desc) {
  if (file.exists(file_path)) {
    f_size <- file.info(file_path)$size
    if (f_size > 0) {
      cat(sprintf("  ✓ VERIFIED: %s (%s, %d bytes)\n", file_path, desc, f_size))
      return(TRUE)
    } else {
      cat(sprintf("  ✗ FAILURE: %s exists but is empty!\n", file_path))
      return(FALSE)
    }
  } else {
    cat(sprintf("  ✗ FAILURE: %s (%s) is missing!\n", file_path, desc))
    return(FALSE)
  }
}

# --- Cleanup existing key deliverables first to ensure a clean test run ---
cat("Cleaning up old files for clean integration testing...\n")
unlink("outputs/elephant_ranges.gpkg")
unlink("outputs/rds/loaded_data.rds")
unlink("outputs/framework1_best_model.RDS")
unlink("outputs/framework2_best_model.RDS")
unlink("figures/figureS1.png")
unlink("figures/figureS1.pdf")
unlink("figures/figure3.png")
unlink("figures/figure3.pdf")
unlink("figures/figure4.png")
unlink("figures/figure4.pdf")
unlink("figures/figure4_alternate_aic_full.png")
unlink("figures/figure4_alternate_aic_full.pdf")
unlink("figures/figureS5.png")
unlink("figures/figureS5.pdf")
unlink("figures/figure5.png")
unlink("figures/figure5.pdf")
unlink("figures/figureS6.png")
unlink("figures/figureS6.pdf")
unlink("figures/figure1_gedi_pipeline.png")
unlink("figures/figure2_camera_trap_pipeline.png")
cat("  Old files cleaned successfully.\n\n")

# =============================================================================
# RUN PIPELINE SEQUENTIALLY
# =============================================================================

# Step 1: Regional Bounding Boxes & Elephant Ranges
cat("--- Step 1: Running Bounding Boxes & Elephant Ranges (01_FigureS1_Regional_Bounding_Boxes.R) ---\n")
source("code/01_FigureS1_Regional_Bounding_Boxes.R")
cat("Step 1 Completed successfully.\n\n")

# Step 2: Vector Loading and Joining
cat("--- Step 2: Running Vector Loading & Joining (02_Load_And_Join.R) ---\n")
source("code/02_Load_And_Join.R")
cat("Step 2 Completed successfully.\n\n")

# Step 3: Framework 1 model selection & plotting
cat("--- Step 3: Running Framework 1 Beta regressions (03_Framework1_Analysis.R) ---\n")
source("code/03_Framework1_Analysis.R")
cat("Step 3 Completed successfully.\n\n")

# Step 4: Framework 2 Tweedie regressions & plotting
cat("--- Step 4: Running Framework 2 Tweedie GLMs (04_Framework2_Analysis.R) ---\n")
source("code/04_Framework2_Analysis.R")
cat("Step 4 Completed successfully.\n\n")

# Step 5: Predictive biomass mapping
cat("--- Step 5: Running Predictive Biomass Mapping (05_Predictive_Biomass_Maps.R) ---\n")
source("code/05_Predictive_Biomass_Maps.R")
cat("Step 5 Completed successfully.\n\n")

# Step 5b: Predictive columns correlation plotting (Supplementary Figure S6)
cat("--- Step 5b: Running Model Predictions Correlation Plotting (05b_Predictive_Columns_Correlation.R) ---\n")
source("code/05b_Predictive_Columns_Correlation.R")
cat("Step 5b Completed successfully.\n\n")

# Step 6: Pipeline Visualization
cat("--- Step 6: Running Pipeline Visualization (06_Pipeline_Visualization.R) ---\n")
source("code/06_Pipeline_Visualization.R")
cat("Step 6 Completed successfully.\n\n")

# =============================================================================
# DELIVERABLE VERIFICATION & INTEGRITY CHECKS
# =============================================================================
cat("============================================================\n")
cat("=== Verifying Generated Deliverables & Figure Integrity ===\n")
cat("============================================================\n\n")

checks <- c(
  # Data deliverables
  verify_file("outputs/elephant_ranges.gpkg", "Combined Elephant Range Vector GPKG"),
  verify_file("outputs/rds/loaded_data.rds", "Consolidated Multi-Scale RDS stack"),
  verify_file("outputs/framework1_covariate_model_selection.csv", "Framework 1 selection table"),
  verify_file("outputs/framework1_best_formula.RDS", "Framework 1 best formula RDS"),
  verify_file("outputs/framework1_best_model.RDS", "Framework 1 best model RDS"),
  verify_file("outputs/framework2_covariate_model_selection.csv", "Framework 2 selection table"),
  verify_file("outputs/framework2_best_formula.RDS", "Framework 2 best formula RDS"),
  verify_file("outputs/framework2_best_model.RDS", "Framework 2 best model RDS"),
  
  # Figure deliverables
  verify_file("figures/figureS1.png", "Manuscript Figure S1 (Bounding Boxes)"),
  verify_file("figures/figure3.png", "Manuscript Figure 3 (F1)"),
  verify_file("figures/figure4.png", "Manuscript Figure 4 (F2)"),
  verify_file("figures/figure4_alternate_aic_full.png", "Manuscript Figure 4 Alternate (AIC-selected)"),
  verify_file("figures/figureS5.png", "Predicted Biomass Map (5km, Fig S5)"),
  verify_file("figures/figure5.png", "Predicted Biomass Map (20km, Fig 5)"),
  verify_file("figures/figureS6.png", "Predictive Columns Correlation Plot (Fig S6)"),
  verify_file("figures/figure1_gedi_pipeline.png", "Manuscript Figure 1 (GEDI Pipeline)"),
  verify_file("figures/figure2_camera_trap_pipeline.png", "Manuscript Figure 2 (Camera Trap Pipeline)")
)

# =============================================================================
# INTEGRATION TEST FINAL STATUS REPORT
# =============================================================================
cat("\n", strrep("=", 60), "\n")
if (all(checks)) {
  cat("  ✓ INTEGRATION STATUS: SUCCESSFUL\n")
  cat("  ✓ The entire pipeline was executed end-to-end sequentially.\n")
  cat("  ✓ All model files, GPKGs, CSV tables, and figures were generated.\n")
  cat("  ✓ Figures populated directly in figures/ as requested.\n")
  cat(strrep("=", 60), "\n")
  quit(status = 0)
} else {
  cat("  ✗ INTEGRATION STATUS: FAILED\n")
  cat("  ✗ Some deliverables are missing or corrupted!\n")
  cat(strrep("=", 60), "\n")
  quit(status = 1)
}
