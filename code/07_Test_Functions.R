# =============================================================================
# 07_Test_Functions.R
#
# Unit test runner for all analysis function files.
#
# Uses a tryCatch-based pattern with passed/failed counters. Tests are
# organised by function file (module) and test category. Each test
# verifies that the function exists, has the correct signature, and
# (where possible) returns the expected structure when run on synthetic data.
#
# Usage:
#   Rscript code/07_Test_Functions.R
#
# Dependencies:
#   terra, dplyr, tibble
# =============================================================================

# --- Setup -------------------------------------------------------------------

library(terra)
library(dplyr)
library(tibble)

cat("=== 07: Unit Tests ===\n\n")

# Counters
tests_passed <- 0L
tests_failed <- 0L
test_results <- list()


# --- Test helpers ------------------------------------------------------------

#' Run a single test with tryCatch
#'
#' Evaluates an expression in a tryCatch block, recording pass/fail.
#' On failure, captures and prints the error message.
#'
#' @param test_name Character. Descriptive name for the test.
#' @param expr An expression to evaluate. Should return TRUE on success.
#'
#' @return Invisible logical. TRUE if the test passed, FALSE otherwise.
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


#' Check that a function exists and has the expected formals
#'
#' @param fn_name Character. Function name.
#' @param expected_args Character vector. Expected argument names.
#' @param env Environment to search in (default: globalenv).
#'
#' @return Logical. TRUE if function exists with correct formals.
check_signature <- function(fn_name, expected_args, env = globalenv()) {
  if (!exists(fn_name, envir = env, mode = "function")) {
    return(FALSE)
  }
  fn <- get(fn_name, envir = env, mode = "function")
  actual_args <- names(formals(fn))
  all(expected_args %in% actual_args)
}


# --- Source all function files -----------------------------------------------

cat("Sourcing function files...\n")
source("code/functions/gedi_analysis.R")
source("code/functions/frip_analysis.R")
source("code/functions/convergence_analysis.R")
source("code/functions/temporal_analysis.R")
source("code/functions/plotting.R")
cat("  All function files sourced.\n\n")


# =============================================================================
# MODULE 1: gedi_analysis.R
# =============================================================================

cat("--- Module: gedi_analysis.R ---\n")

run_test("test_regional_uoi exists", quote(
  exists("test_regional_uoi", mode = "function")
))

run_test("test_regional_uoi signature", quote(
  check_signature("test_regional_uoi", c("native_stacks"))
))

run_test("test_regional_uoi returns error (not implemented)", quote(
  tryCatch(
    { test_regional_uoi(list()); FALSE },
    error = function(e) grepl("Not yet implemented", e$message)
  )
))

run_test("test_uoi_region_protection exists", quote(
  exists("test_uoi_region_protection", mode = "function")
))

run_test("test_uoi_region_protection signature", quote(
  check_signature("test_uoi_region_protection", c("native_stacks", "pa_rast"))
))

run_test("test_uoi_by_pa exists", quote(
  exists("test_uoi_by_pa", mode = "function")
))

run_test("test_uoi_by_pa signature", quote(
  check_signature("test_uoi_by_pa", c("native_stacks", "pa_rast"))
))

run_test("run_h1_multiscale exists", quote(
  exists("run_h1_multiscale", mode = "function")
))

run_test("run_h1_multiscale signature", quote(
  check_signature("run_h1_multiscale", c("multiscale_stacks", "pa_rast"))
))


# =============================================================================
# MODULE 2: frip_analysis.R
# =============================================================================

cat("\n--- Module: frip_analysis.R ---\n")

run_test("test_frip_by_basin_country_pa exists", quote(
  exists("test_frip_by_basin_country_pa", mode = "function")
))

run_test("test_frip_by_basin_country_pa signature", quote(
  check_signature("test_frip_by_basin_country_pa",
                   c("stack", "basins_r", "countries_r", "pa_r"))
))

run_test("test_frip_vs_di exists", quote(
  exists("test_frip_vs_di", mode = "function")
))

run_test("test_frip_vs_di signature", quote(
  check_signature("test_frip_vs_di", c("stack", "di_rast"))
))

run_test("compare_di_indices exists", quote(
  exists("compare_di_indices", mode = "function")
))

run_test("compare_di_indices signature", quote(
  check_signature("compare_di_indices", c("stack", "di_bl", "di_bogoni"))
))

run_test("run_denoising_cv exists", quote(
  exists("run_denoising_cv", mode = "function")
))

run_test("run_denoising_cv signature", quote(
  check_signature("run_denoising_cv", c("stack", "di_rast", "basins_v"))
))

run_test("run_basin_split_denoising exists", quote(
  exists("run_basin_split_denoising", mode = "function")
))

run_test("run_basin_split_denoising signature", quote(
  check_signature("run_basin_split_denoising", c("stack", "di_rast", "basins_v"))
))

run_test("run_h2_multiscale exists", quote(
  exists("run_h2_multiscale", mode = "function")
))

run_test("run_h2_multiscale signature", quote(
  check_signature("run_h2_multiscale", c("multiscale_stacks"))
))


# =============================================================================
# MODULE 3: convergence_analysis.R
# =============================================================================

cat("\n--- Module: convergence_analysis.R ---\n")

run_test("test_pa_convergence exists", quote(
  exists("test_pa_convergence", mode = "function")
))

run_test("test_pa_convergence signature", quote(
  check_signature("test_pa_convergence", c("native_stacks", "pa_rast"))
))

run_test("test_pixel_convergence exists", quote(
  exists("test_pixel_convergence", mode = "function")
))

run_test("test_pixel_convergence signature", quote(
  check_signature("test_pixel_convergence", c("native_stacks"))
))

run_test("classify_bivariate exists", quote(
  exists("classify_bivariate", mode = "function")
))

run_test("classify_bivariate signature", quote(
  check_signature("classify_bivariate", c("native_stacks"))
))

run_test("test_convergence_multiscale exists", quote(
  exists("test_convergence_multiscale", mode = "function")
))

run_test("test_convergence_multiscale signature", quote(
  check_signature("test_convergence_multiscale",
                   c("multiscale_stacks", "pa_rast"))
))

run_test("test_protection_dual_signal exists", quote(
  exists("test_protection_dual_signal", mode = "function")
))

run_test("test_protection_dual_signal signature", quote(
  check_signature("test_protection_dual_signal",
                   c("native_stacks", "pa_rast"))
))


# =============================================================================
# MODULE 4: temporal_analysis.R
# =============================================================================

cat("\n--- Module: temporal_analysis.R ---\n")

run_test("summarise_mk_tau_by_basin exists", quote(
  exists("summarise_mk_tau_by_basin", mode = "function")
))

run_test("summarise_mk_tau_by_basin signature", quote(
  check_signature("summarise_mk_tau_by_basin", c("stacks"))
))

run_test("test_mk_tau_by_pa exists", quote(
  exists("test_mk_tau_by_pa", mode = "function")
))

run_test("test_mk_tau_by_pa signature", quote(
  check_signature("test_mk_tau_by_pa", c("stacks", "pa_rast"))
))

run_test("test_mk_tau_protection exists", quote(
  exists("test_mk_tau_protection", mode = "function")
))

run_test("test_mk_tau_protection signature", quote(
  check_signature("test_mk_tau_protection", c("stacks", "pa_rast"))
))

run_test("run_h4_multiscale exists", quote(
  exists("run_h4_multiscale", mode = "function")
))

run_test("run_h4_multiscale signature", quote(
  check_signature("run_h4_multiscale", c("multiscale_stacks", "pa_rast"))
))


# =============================================================================
# MODULE 5: plotting.R
# =============================================================================

cat("\n--- Module: plotting.R ---\n")

run_test("make_basin_map exists", quote(
  exists("make_basin_map", mode = "function")
))

run_test("make_basin_map signature", quote(
  check_signature("make_basin_map", c("rast", "fill_col", "scale_fn"))
))

run_test("make_paired_maps exists", quote(
  exists("make_paired_maps", mode = "function")
))

run_test("make_paired_maps signature", quote(
  check_signature("make_paired_maps",
                   c("rast_congo", "rast_amazon", "fill_col", "scale_fn",
                     "countries"))
))

run_test("make_boxplot_with_letters exists", quote(
  exists("make_boxplot_with_letters", mode = "function")
))

run_test("make_boxplot_with_letters signature", quote(
  check_signature("make_boxplot_with_letters", c("df", "x", "y", "letters_df"))
))

run_test("make_pa_pairs_bar exists", quote(
  exists("make_pa_pairs_bar", mode = "function")
))

run_test("make_ranked_parks_boxplot exists", quote(
  exists("make_ranked_parks_boxplot", mode = "function")
))

run_test("make_openness_distribution exists", quote(
  exists("make_openness_distribution", mode = "function")
))

run_test("make_multiscale_ci_plot exists", quote(
  exists("make_multiscale_ci_plot", mode = "function")
))

run_test("make_multiscale_ci_plot signature", quote(
  check_signature("make_multiscale_ci_plot",
                   c("df", "x_col", "ymin_col", "ymax_col", "signif_col"))
))

run_test("make_denoising_ratio_plot exists", quote(
  exists("make_denoising_ratio_plot", mode = "function")
))

run_test("make_bivariate_map exists", quote(
  exists("make_bivariate_map", mode = "function")
))

run_test("make_bivariate_map signature", quote(
  check_signature("make_bivariate_map", c("rast", "countries"))
))

run_test("make_scatter_with_cor exists", quote(
  exists("make_scatter_with_cor", mode = "function")
))

run_test("make_scatter_with_cor signature", quote(
  check_signature("make_scatter_with_cor", c("df", "x", "y", "color_col"))
))

run_test("assemble_figure exists", quote(
  exists("assemble_figure", mode = "function")
))

run_test("save_pnas exists", quote(
  exists("save_pnas", mode = "function")
))

run_test("save_pnas signature", quote(
  check_signature("save_pnas", c("plot", "filename", "type"))
))

# --- Constant checks ---
run_test("BASIN_COLOURS defined", quote(
  exists("BASIN_COLOURS") && length(BASIN_COLOURS) == 2
))

run_test("BIVARIATE_COLOURS defined", quote(
  exists("BIVARIATE_COLOURS") && length(BIVARIATE_COLOURS) == 4
))


# =============================================================================
# SUMMARY
# =============================================================================

cat("\n", strrep("=", 60), "\n")
total <- tests_passed + tests_failed
cat(sprintf("  TOTAL: %d tests | PASSED: %d | FAILED: %d\n",
            total, tests_passed, tests_failed))

if (tests_failed == 0L) {
  cat("  ✓ All tests passed!\n")
} else {
  cat("  ✗ Some tests failed. Review output above.\n")
  failed_tests <- names(which(sapply(test_results, function(x) x != "PASS")))
  cat("  Failed tests:\n")
  for (ft in failed_tests) {
    cat(sprintf("    - %s: %s\n", ft, test_results[[ft]]))
  }
}
cat(strrep("=", 60), "\n")
