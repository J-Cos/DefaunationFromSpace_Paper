# =============================================================================
# code/04_Framework2_Analysis.R
#
# Performs covariate model selection for Framework 2 (Spaceborne Biomass Prediction)
# using Leave-One-Basin-Out (LOBO) cross-validation across 37 candidate Tweedie GLMs
# (comparing strict vs possible elephant definitions and pure UOI), selects the best
# parsimonious model via the 1-SE rule.
# Figures are generated separately by 04b_Framework_Figures.R.
# =============================================================================

library(terra)
library(dplyr)
library(readr)
library(mgcv)
library(ggplot2)

source("code/functions/calibration_helpers.R")
source("code/functions/model_convergence.R")
source("code/functions/framework_helpers.R")

#' Run Framework 2 Analysis
#'
#' Fits a series of Tweedie GLMs predicting standing mammal biomass from GEDI understory
#' openness, conducts Leave-One-Basin-Out (LOBO) cross-validation, applies the 1-SE
#' parsimony strategy, and outputs diagnostic plots.
#'
#' @param scale_m Numeric. Spatial resolution in meters (default: 5000)
#' @param outputs_dir Character. Directory to save model outputs (default: "outputs")
#'
#' @return A list containing the model selection results data frame and the best fitted model object.
#' @export
run_framework2_analysis <- function(scale_m = 5000, outputs_dir = "outputs") {
  cat(sprintf("=== Starting Framework 2 Model Selection (%d m scale) ===\n\n", scale_m))
  
  # --- 1. Ingest and Calibrate Scale-Specific Cluster Data ---------------------
  joined_data <- extract_scale_data(scale_m)
  
  # Clean up any potential NA values in key environmental variables
  covs_to_check <- c("uoi", "elevation", "slope", "hnd", "precip", "clay", "forest_fraction", "B_H_index", "elephant_present_strict", "elephant_present_possible")
  joined_data <- joined_data %>% filter(complete.cases(joined_data[, covs_to_check]))
  
  cat("✓ Merged and calibrated data successfully. N =", nrow(joined_data), "clusters.\n")
  
  # --- 2. Define the 37 Candidate Covariate Models (Excluding Basin & MegaHx) ----
  formulas_list <- framework2_formulas
  
  expanded_formulas_list <- formulas_list

  num_models <- length(expanded_formulas_list)
  model_names <- names(expanded_formulas_list)
  
  cat(sprintf("Fitting all %d models on the full dataset...\n", num_models))
  full_models <- list()
  full_AICc <- numeric(num_models)
  full_edf <- numeric(num_models)
  full_dev_expl <- numeric(num_models)
  n_obs <- nrow(joined_data)
  
  for (m_idx in 1:num_models) {
    m_name <- model_names[m_idx]
    m_form <- expanded_formulas_list[[m_idx]]
    
    m_full <- gam(m_form, family = tw(), weights = w_combined_norm, data = joined_data, method = "ML")
    check_model_convergence(m_full, m_name)
    full_models[[m_name]] <- m_full
    
    aic_val <- AIC(m_full)
    k_m <- attr(logLik(m_full), "df")
    full_AICc[m_idx] <- aic_val + (2 * k_m * (k_m + 1)) / (n_obs - k_m - 1)
    full_edf[m_idx] <- sum(m_full$edf)
    full_dev_expl[m_idx] <- summary(m_full)$dev.expl
  }
  
  # --- 4. Perform Leave-One-Basin-Out (LOBO) Cross-Validation -----------------
  cat("Running Leave-One-Basin-Out (LOBO) Cross-Validation (Skipping Basin Dummies)...\n")
  basins <- unique(as.character(joined_data$basin))
  
  oos_predictions <- matrix(NA, nrow = nrow(joined_data), ncol = num_models)
  colnames(oos_predictions) <- model_names
  
  for (b in basins) {
    train_idx <- which(joined_data$basin != b)
    test_idx <- which(joined_data$basin == b)
    
    train_data <- joined_data[train_idx, ]
    test_data <- joined_data[test_idx, ]
    
    train_data$w_combined_norm <- train_data$w_combined / mean(train_data$w_combined)
    
    for (m_idx in 1:num_models) {
      m_form <- expanded_formulas_list[[m_idx]]
      m_name <- model_names[m_idx]
      
      # Skip LOBO cross-validation for any model containing "basin" to prevent unseen category factor errors
      if ("basin" %in% all.vars(m_form)) {
        next
      }
      
      fold_model <- tryCatch({
        m_fold <- gam(m_form, family = tw(), weights = w_combined_norm, data = train_data, method = "ML")
        check_model_convergence(m_fold, paste(m_name, "fold", b), raise_warning = FALSE)
        m_fold
      }, error = function(e) {
        NULL
      })
      
      if (!is.null(fold_model)) {
        pred <- tryCatch({
          predict(fold_model, newdata = test_data, type = "response")
        }, error = function(e) {
          rep(NA, nrow(test_data))
        })
        oos_predictions[test_idx, m_idx] <- pred
      }
    }
  }
  
  # --- 5. Calculate Out-of-Sample Metrics --------------------------------------
  y_obs <- joined_data$B_H_index
  log_y_obs <- log1p(y_obs)
  mean_log_y_obs <- mean(log_y_obs)
  
  oos_RMSE_log <- numeric(num_models)
  oos_MAE_log <- numeric(num_models)
  oos_RMSE_raw <- numeric(num_models)
  oos_MAE_raw <- numeric(num_models)
  oos_R2_log <- numeric(num_models)
  
  for (m_idx in 1:num_models) {
    y_pred <- oos_predictions[, m_idx]
    non_na_idx <- which(!is.na(y_pred))
    n_valid <- length(non_na_idx)
    
    if (n_valid > 0) {
      y_pred_valid <- y_pred[non_na_idx]
      y_obs_valid <- y_obs[non_na_idx]
      
      oos_RMSE_raw[m_idx] <- sqrt(mean((y_obs_valid - y_pred_valid)^2))
      oos_MAE_raw[m_idx] <- mean(abs(y_obs_valid - y_pred_valid))
      
      log_y_pred_valid <- log1p(y_pred_valid)
      log_y_obs_valid <- log_y_obs[non_na_idx]
      
      oos_RMSE_log[m_idx] <- sqrt(mean((log_y_obs_valid - log_y_pred_valid)^2))
      oos_MAE_log[m_idx] <- mean(abs(log_y_obs_valid - log_y_pred_valid))
      
      # NOTE: R² is computed relative to the GLOBAL mean (mean_log_y_obs from
      # the full dataset, L170), not the held-out fold's mean. This is a
      # deliberate design choice: since all models share the same denominator
      # (TSS), comparison across models within the same LOBO framework is valid.
      # Using fold-specific means would give a pure generalization metric but
      # would make R² values less comparable across different fold sizes.
      rss_log_model <- sum((log_y_obs_valid - log_y_pred_valid)^2)
      tss_log_valid <- sum((log_y_obs_valid - mean_log_y_obs)^2)
      oos_R2_log[m_idx] <- 1 - (rss_log_model / tss_log_valid)
    } else {
      oos_RMSE_raw[m_idx] <- NA; oos_MAE_raw[m_idx] <- NA
      oos_RMSE_log[m_idx] <- NA; oos_MAE_log[m_idx] <- NA
      oos_R2_log[m_idx] <- NA
    }
  }
  
  oos_MAE_log_SE <- sapply(1:num_models, function(m_idx) {
    abs_err <- abs(log_y_obs - log1p(oos_predictions[, m_idx]))
    if (all(is.na(abs_err))) return(NA)
    sd(abs_err, na.rm = TRUE) / sqrt(sum(!is.na(abs_err)))
  })
  
  # 1-SE Parsimony Selection Strategy (evaluating only valid OOS models)
  valid_idx <- which(!is.na(oos_MAE_log))
  if (length(valid_idx) > 0) {
    raw_best_idx <- valid_idx[which.min(oos_MAE_log[valid_idx])]
    best_mae <- oos_MAE_log[raw_best_idx]
    best_se <- oos_MAE_log_SE[raw_best_idx]
    threshold_1se <- best_mae + best_se
    
    compliant_indices <- valid_idx[which(oos_MAE_log[valid_idx] <= threshold_1se)]
    min_edf <- min(round(full_edf[compliant_indices], 4))
    best_parsimonious_indices <- compliant_indices[which(round(full_edf[compliant_indices], 4) == min_edf)]
    selected_idx <- best_parsimonious_indices[which.min(oos_MAE_log[best_parsimonious_indices])]
  } else {
    raw_best_idx <- NA
    selected_idx <- NA
    threshold_1se <- NA
  }
  
  parsimony_status <- sapply(1:num_models, function(m_idx) {
    if (is.na(oos_MAE_log[m_idx])) {
      return("Failed LOBOCV (Collinear/Unseen Levels)")
    }
    mae <- oos_MAE_log[m_idx]
    edf <- round(full_edf[m_idx], 4)
    raw_best_edf <- round(full_edf[raw_best_idx], 4)
    
    if (m_idx == selected_idx) {
      return("Parsimonious Selected Best")
    } else if (m_idx == raw_best_idx) {
      return("Raw Best (More Complex)")
    } else if (mae <= threshold_1se) {
      if (edf < raw_best_edf) {
        return("Parsimonious Candidate (Within 1-SE)")
      } else {
        return("Equivalent (Within 1-SE, More Complex)")
      }
    } else {
      return("Suboptimal (Outside 1-SE)")
    }
  })
  
  # Calculate AICc delta and status for all candidate models
  results_df_aic <- data.frame(
    Model = model_names,
    AICc = full_AICc,
    edf = full_edf,
    stringsAsFactors = FALSE
  )
  results_df_aic$delta_AICc <- results_df_aic$AICc - min(results_df_aic$AICc)
  best_aic_idx <- which.min(results_df_aic$AICc)
  
  AIC_status <- sapply(1:num_models, function(m_idx) {
    if (m_idx == best_aic_idx) {
      return("AICc Selected Best")
    } else if (results_df_aic$delta_AICc[m_idx] <= 2) {
      return("AICc Equivalent (delta <= 2)")
    } else if (results_df_aic$delta_AICc[m_idx] <= 7) {
      return("AICc Suboptimal (delta <= 7)")
    } else {
      return("AICc Poor (delta > 7)")
    }
  })
  
  # --- 6. Compile and Save Selection Table -------------------------------------
  results_df <- data.frame(
    Model = model_names,
    OOS_RMSE_log = oos_RMSE_log,
    OOS_MAE_log = oos_MAE_log,
    OOS_MAE_log_SE = oos_MAE_log_SE,
    OOS_R2_log = oos_R2_log,
    OOS_RMSE_raw = oos_RMSE_raw,
    OOS_MAE_raw = oos_MAE_raw,
    Full_AICc = full_AICc,
    delta_AICc = results_df_aic$delta_AICc,
    Full_edf = full_edf,
    Full_DevExpl = full_dev_expl,
    Parsimony_Status = parsimony_status,
    AICc_Status = AIC_status,
    stringsAsFactors = FALSE
  )
  
  # Sorting strategy centered on 1-SE rank
  results_df$sort_rank <- sapply(results_df$Parsimony_Status, function(status) {
    if (status == "Parsimonious Selected Best") return(1)
    if (status == "Raw Best (More Complex)") return(2)
    if (status == "Parsimonious Candidate (Within 1-SE)") return(3)
    if (status == "Equivalent (Within 1-SE, More Complex)") return(4)
    if (status == "Suboptimal (Outside 1-SE)") return(5)
    return(6)
  })
  
  results_df <- results_df %>%
    arrange(sort_rank, OOS_MAE_log) %>%
    select(-sort_rank)
  
  print(results_df)
  
  dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)

  
  write_csv(results_df, file.path(outputs_dir, "framework2_covariate_model_selection.csv"))
  
  # ───────────────────────────────────────────────────────────────────────────
  # RATIONALE FOR DUAL MODEL RETENTION (HI-7):
  # The pipeline evaluates models using two distinct criteria for different purposes:
  #
  # 1. Leave-One-Basin-Out Cross-Validation (LOBO-CV) / LOROCV:
  #    - Purpose: Identifies generalizability across unseen geographic basins.
  #    - Behavior: Selects the most parsimonious model (e.g. M2.1: UOI Only) 
  #      under a 1-SE rule. This represents the universal biophysical baseline.
  #
  # 2. Standard AICc Model Selection:
  #    - Purpose: Identifies the best fitting model locally within the sampled 
  #      basins, incorporating regional interactions (e.g. UOI * Basin + Elev).
  #    - Behavior: Captures biogeographical shifts and local ecological modifiers 
  #      (e.g., presence of elephants in Congo/SE_Asia vs absence in Amazon).
  #
  # Retaining both allows side-by-side comparison (e.g., Column 1 vs Column 2 
  # in predictive biomass mapping) of universal baseline vs local biogeographic 
  # calibration.
  # ───────────────────────────────────────────────────────────────────────────

  # Save RDS best model objects (LOBOCV parsimonious model) — refit with REML for inference
  best_model_name <- results_df$Model[which(results_df$Parsimony_Status == "Parsimonious Selected Best")]
  best_model <- gam(formula(full_models[[best_model_name]]), family = tw(),
                    weights = w_combined_norm, data = joined_data, method = "REML")
  check_model_convergence(best_model, paste(best_model_name, "(REML Refit)"))
  cat(sprintf("\n★ Parsimoniously Selected Best Model (1-SE Strategy): %s (OOS MAE_log: %.4f)\n", best_model_name, results_df$OOS_MAE_log[which(results_df$Model == best_model_name)]))
  cat("  (Selection via ML + LOBO-CV; coefficients from REML refit)\n\n")
  
  saveRDS(formula(best_model), file.path(outputs_dir, "framework2_best_formula.RDS"))
  saveRDS(best_model, file.path(outputs_dir, "framework2_best_model.RDS"))
  
  # Save RDS best model objects (Standard AIC model) — refit with REML for inference
  best_model_name_aic <- results_df_aic$Model[best_aic_idx]
  best_model_aic <- gam(formula(full_models[[best_model_name_aic]]), family = tw(),
                        weights = w_combined_norm, data = joined_data, method = "REML")
  check_model_convergence(best_model_aic, paste(best_model_name_aic, "(REML Refit)"))
  saveRDS(formula(best_model_aic), file.path(outputs_dir, "framework2_best_formula_aic.RDS"))
  saveRDS(best_model_aic, file.path(outputs_dir, "framework2_best_model_aic.RDS"))
  
  saveRDS(full_models, file.path(outputs_dir, "framework2_full_models.RDS"))
  
  cat("=== Framework 2 Analysis Completed Successfully ===\n")
  
  return(list(
    results_df = results_df,
    best_model = best_model,
    full_models = full_models
  ))
}

# --- Execute directly if called from terminal ---
run_framework2_analysis()
