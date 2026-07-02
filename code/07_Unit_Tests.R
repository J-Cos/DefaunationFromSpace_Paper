# =============================================================================
# code/07_Unit_Tests.R
#
# Rigorous unit testing suite for all central functions in the defaunation pipeline:
#   1. calibration_helpers.R functions
#   2. Analytical scripts functions (run_framework1_analysis, run_framework2_analysis,
#      run_predictive_biomass_mapping)
#
# Returns standard exit codes (0 for pass, 1 for fail).
# =============================================================================

library(terra)
library(dplyr)
library(mgcv)

cat("=== Starting Comprehensive Pipeline Unit Tests ===\n\n")

# Counters
tests_passed <- 0L
tests_failed <- 0L
test_results <- list()

# --- Test helpers ------------------------------------------------------------
run_test <- function(test_name, expr) {
  result <- tryCatch(
    {
      val <- eval(expr)
      if (isTRUE(val)) {
        tests_passed <<- tests_passed + 1L
        cat(sprintf("  ✓ PASS: %s\n", test_name))
        test_results[[test_name]] <<- "PASS"
        TRUE
      } else {
        tests_failed <<- tests_failed + 1L
        cat(sprintf("  ✗ FAIL: %s (returned %s)\n", test_name, as.character(val)))
        test_results[[test_name]] <<- "FAIL"
        FALSE
      }
    },
    error = function(e) {
      tests_failed <<- tests_failed + 1L
      cat(sprintf("  ✗ FAIL: %s (error: %s)\n", test_name, conditionMessage(e)))
      test_results[[test_name]] <<- paste("ERROR:", conditionMessage(e))
      FALSE
    }
  )
  invisible(result)
}

check_signature <- function(fn_name, expected_args, env = globalenv()) {
  if (!exists(fn_name, envir = env, mode = "function")) {
    return(FALSE)
  }
  fn <- get(fn_name, envir = env, mode = "function")
  actual_args <- names(formals(fn))
  all(expected_args %in% actual_args)
}

# --- Create clean temp directories for sandbox testing ---
temp_out_dir <- "outputs/temp_test_out"
temp_fig_dir <- "outputs/temp_test_fig"
dir.create(temp_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(temp_fig_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# MODULE 1: code/functions/calibration_helpers.R
# =============================================================================
cat("\n--- Module: code/functions/calibration_helpers.R ---\n")
source("code/functions/calibration_helpers.R")

run_test("extract_scale_pixels exists", quote({
  exists("extract_scale_pixels", mode = "function")
}))

run_test("extract_scale_pixels signature check", quote({
  check_signature("extract_scale_pixels", c("scale_m", "mcps"))
}))

run_test("extract_scale_data exists", quote({
  exists("extract_scale_data", mode = "function")
}))

run_test("extract_scale_data signature check", quote({
  check_signature("extract_scale_data", c("scale_m", "mcps"))
}))


run_test("fit_framework1_model exists", quote({
  exists("fit_framework1_model", mode = "function")
}))

run_test("fit_framework1_model signature check", quote({
  check_signature("fit_framework1_model", c("data", "formula_path"))
}))

run_test("fit_framework2_model exists", quote({
  exists("fit_framework2_model", mode = "function")
}))

run_test("fit_framework2_model signature check", quote({
  check_signature("fit_framework2_model", c("data", "formula_path"))
}))

# Behavioral test for extract_scale_data on synthetic data
run_test("extract_scale_data returns valid structure (5km scale)", quote({
  dat <- extract_scale_data(5000)
  is.data.frame(dat) && 
    nrow(dat) > 0 && 
    all(c("cluster_id", "region", "basin", "trap_days", "uoi", "w_combined_norm", "homogeneity", "w_temp_cluster") %in% names(dat))
}))

# =============================================================================
# MODULE 2: code/03_Framework1_Analysis.R
# =============================================================================
cat("\n--- Module: code/03_Framework1_Analysis.R ---\n")
source("code/03_Framework1_Analysis.R")

run_test("run_framework1_analysis exists", quote({
  exists("run_framework1_analysis", mode = "function")
}))

run_test("run_framework1_analysis signature check", quote({
  check_signature("run_framework1_analysis", c("scale_m", "outputs_dir", "figures_dir"))
}))

run_test("run_framework1_analysis execution and returns list of data + fitted model", quote({
  res <- run_framework1_analysis(scale_m = 5000, outputs_dir = temp_out_dir, figures_dir = temp_fig_dir)
  
  is.list(res) &&
    is.data.frame(res$results_df) &&
    inherits(res$best_model, "gam") &&
    grepl("Beta regression", res$best_model$family$family) &&
    file.exists(file.path(temp_out_dir, "framework1_covariate_model_selection.csv")) &&
    file.exists(file.path(temp_fig_dir, "figure3.png"))
}))

# =============================================================================
# MODULE 3: code/04_Framework2_Analysis.R
# =============================================================================
cat("\n--- Module: code/04_Framework2_Analysis.R ---\n")
source("code/04_Framework2_Analysis.R")

run_test("run_framework2_analysis exists", quote({
  exists("run_framework2_analysis", mode = "function")
}))

run_test("run_framework2_analysis signature check", quote({
  check_signature("run_framework2_analysis", c("scale_m", "outputs_dir", "figures_dir"))
}))

run_test("run_framework2_analysis execution and returns list of data + fitted Tweedie GLM", quote({
  res <- run_framework2_analysis(scale_m = 5000, outputs_dir = temp_out_dir, figures_dir = temp_fig_dir)
  
  is.list(res) &&
    is.data.frame(res$results_df) &&
    inherits(res$best_model, "gam") &&
    grepl("Tweedie", res$best_model$family$family) &&
    file.exists(file.path(temp_out_dir, "framework2_covariate_model_selection.csv")) &&
    file.exists(file.path(temp_fig_dir, "figure4.png"))
}))

# =============================================================================
# MODULE 4: code/05_Predictive_Biomass_Maps.R
# =============================================================================
cat("\n--- Module: code/05_Predictive_Biomass_Maps.R ---\n")
source("code/05_Predictive_Biomass_Maps.R")

run_test("run_predictive_biomass_mapping exists", quote({
  exists("run_predictive_biomass_mapping", mode = "function")
}))

run_test("run_predictive_biomass_mapping signature check", quote({
  check_signature("run_predictive_biomass_mapping", c("scales", "outputs_dir", "figures_dir"))
}))

run_test("run_predictive_biomass_mapping successfully projects and outputs rasters and figures", quote({
  # Copy the best model RDS, AIC best model RDS, and selection table to the outputs folder first so mapping function can read it
  file.copy(file.path(temp_out_dir, "framework2_best_model.RDS"),
            file.path("outputs", "framework2_best_model.RDS"),
            overwrite = TRUE)
  file.copy(file.path(temp_out_dir, "framework2_best_model_aic.RDS"),
            file.path("outputs", "framework2_best_model_aic.RDS"),
            overwrite = TRUE)
  file.copy(file.path(temp_out_dir, "framework2_covariate_model_selection.csv"),
            file.path("outputs", "framework2_covariate_model_selection.csv"),
            overwrite = TRUE)
  
  rasts <- run_predictive_biomass_mapping(scales = c(5000, 20000), outputs_dir = "outputs", figures_dir = temp_fig_dir)
  
  is.list(rasts) &&
    length(rasts) == 2 &&
    inherits(rasts[["5000"]]$congo, "SpatRaster") &&
    inherits(rasts[["5000"]]$amazon, "SpatRaster") &&
    file.exists(file.path(temp_fig_dir, "figureS3.png")) &&
    file.exists(file.path(temp_fig_dir, "figure5.png"))
}))

# --- Cleanup temp test directories ---
unlink(temp_out_dir, recursive = TRUE)
unlink(temp_fig_dir, recursive = TRUE)

# =============================================================================
# SUMMARY REPORT
# =============================================================================
cat("\n", strrep("=", 60), "\n")
total_tests <- tests_passed + tests_failed
cat(sprintf("  TOTAL: %d tests | PASSED: %d | FAILED: %d\n",
            total_tests, tests_passed, tests_failed))

if (tests_failed == 0L) {
  cat("  ✓ All pipeline unit tests passed successfully!\n")
  cat(strrep("=", 60), "\n")
  quit(status = 0)
} else {
  cat("  ✗ Some unit tests failed. Review failure output above.\n")
  failed_list <- names(which(sapply(test_results, function(x) x != "PASS")))
  for (ft in failed_list) {
    cat(sprintf("    - %s: %s\n", ft, test_results[[ft]]))
  }
  cat(strrep("=", 60), "\n")
  quit(status = 1)
}
