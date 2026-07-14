# =============================================================================
# code/03_Framework1_Analysis.R
#
# Performs covariate model selection for Framework 1 (GEDI UOI as Response) using
# Beta Regression (via mgcv::gam) across 37 candidate models, evaluates them via AIC.
# Figures are generated separately by 04b_Framework_Figures.R.
#
# All logic is encapsulated in a clean, unit-testable function.
# =============================================================================

library(terra)
library(dplyr)
library(readr)
library(mgcv)
library(ggplot2)

source("code/functions/calibration_helpers.R")
source("code/functions/model_convergence.R")
source("code/functions/framework_helpers.R")

#' Run Framework 1 Analysis
#'
#' Runs GEDI understory openness (UOI) beta regressions, performs model selection
#' on AIC, and saves model data.
#'
#' @param scale_m Numeric. Spatial resolution in meters (default: 5000)
#' @param outputs_dir Character. Directory to save model outputs (default: "outputs")
#'
#' @return A list containing the model selection results data frame and the best fitted model object.
#' @export
run_framework1_analysis <- function(scale_m = 5000, outputs_dir = "outputs") {
  cat(sprintf("=== Starting Framework 1 Model Selection (%d m scale) ===\n\n", scale_m))
  
  # --- 1. Ingest and Calibrate Scale-Specific Cluster Data ---------------------
  joined_data <- extract_scale_data(scale_m)
  
  cat("✓ Merged and calibrated data successfully. N =", nrow(joined_data), "clusters.\n")
  
  cat("\nRunning Beta Regression Covariate Model Selection (Excluding Elephant & MegaHx)...\n")
  base_formulas <- framework1_formulas

  expanded_base_formulas <- base_formulas

  # Fit each base model
  models_list <- list()
  for (name in names(expanded_base_formulas)) {
    f <- expanded_base_formulas[[name]]
    m <- gam(f, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "ML")
    check_model_convergence(m, name)
    models_list[[name]] <- m
  }
  
  # Compile results table (models fitted with ML for valid AICc comparison)
  n_obs <- nrow(joined_data)
  results_df <- data.frame(
    Model = names(models_list),
    AICc = sapply(models_list, function(m) {
      aic_val <- AIC(m)
      k <- attr(logLik(m), "df")
      aic_val + (2 * k * (k + 1)) / (n_obs - k - 1)
    }),
    LogLik = sapply(models_list, function(m) as.numeric(logLik(m))),
    edf = sapply(models_list, function(m) sum(m$edf)),
    R2 = sapply(models_list, function(m) {
      r2 <- summary(m)$r.sq
      if (is.null(r2) || is.na(r2)) return(summary(m)$dev.expl)
      return(r2)
    }),
    dev_expl = sapply(models_list, function(m) summary(m)$dev.expl),
    IsSignificant = sapply(models_list, function(m) {
      p_table <- summary(m)$p.table
      non_intercept_rows <- which(rownames(p_table) != "(Intercept)")
      if (length(non_intercept_rows) == 0) return(FALSE)
      p_vals <- p_table[non_intercept_rows, ncol(p_table), drop = TRUE]
      any(p_vals < 0.05)
    }),
    stringsAsFactors = FALSE
  )
  
  results_df <- results_df %>%
    mutate(delta_AICc = AICc - min(AICc)) %>%
    arrange(AICc)
  
  print(results_df)
  
  dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)
  
  write_csv(results_df, file.path(outputs_dir, "framework1_covariate_model_selection.csv"))
  saveRDS(models_list, file.path(outputs_dir, "framework1_full_models.RDS"))
  
  # Identify the best model (by AICc rank) and refit with REML for inference
  best_model_name <- results_df$Model[1]
  best_model_ml <- models_list[[best_model_name]]
  best_model <- gam(formula(best_model_ml), family = betar(link = "logit"),
                    weights = w_combined_norm, data = joined_data, method = "REML")
  check_model_convergence(best_model, paste(best_model_name, "(REML Refit)"))
  cat(sprintf("\n★ Selected Best-Fitting Model: %s (AICc: %.2f, d_AICc: 0.00)\n", best_model_name, results_df$AICc[1]))
  cat("  (Selection via ML/AICc; coefficients from REML refit)\n\n")
  
  # Save best model RDS objects to outputs
  saveRDS(formula(best_model), file.path(outputs_dir, "framework1_best_formula.RDS"))
  saveRDS(best_model, file.path(outputs_dir, "framework1_best_model.RDS"))
  
  cat("=== Framework 1 Analysis Completed Successfully ===\n")
  
  return(list(
    results_df = results_df,
    best_model = best_model,
    models_list = models_list
  ))
}

# --- Execute directly if called from terminal ---
run_framework1_analysis()
