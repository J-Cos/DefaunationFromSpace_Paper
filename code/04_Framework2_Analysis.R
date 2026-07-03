# =============================================================================
# code/04_Framework2_Analysis.R
#
# Performs covariate model selection for Framework 2 (Spaceborne Biomass Prediction)
# using Leave-One-Basin-Out (LOBO) cross-validation across 37 candidate Tweedie GLMs
# (comparing strict vs possible elephant definitions and pure UOI), selects the best
# parsimonious model via the 1-SE rule, and generates a publication-quality PNAS figure.
# =============================================================================

library(terra)
library(dplyr)
library(readr)
library(mgcv)
library(ggplot2)
library(scales)
library(cowplot)

#' Run Framework 2 Analysis
#'
#' Fits a series of Tweedie GLMs predicting standing mammal biomass from GEDI understory
#' openness, conducts Leave-One-Basin-Out (LOBO) cross-validation, applies the 1-SE
#' parsimony strategy, and outputs diagnostic plots.
#'
#' @param scale_m Numeric. Spatial resolution in meters (default: 5000)
#' @param outputs_dir Character. Directory to save model outputs (default: "outputs")
#' @param figures_dir Character. Directory to save figures (default: "figures")
#'
#' @return A list containing the model selection results data frame and the best fitted model object.
#' @export
run_framework2_analysis <- function(scale_m = 5000, outputs_dir = "outputs", figures_dir = "figures") {
  cat(sprintf("=== Starting Framework 2 Model Selection & Plotting (%d m scale) ===\n\n", scale_m))
  
  # --- 1. Ingest and Calibrate Scale-Specific Cluster Data ---------------------
  source("code/functions/calibration_helpers.R")
  source("code/functions/model_convergence.R")
  joined_data <- extract_scale_data(scale_m)
  
  # Clean up any potential NA values in key environmental variables
  covs_to_check <- c("uoi", "elevation", "slope", "hnd", "precip", "clay", "forest_fraction", "B_H_index", "elephant_present_strict", "elephant_present_possible")
  joined_data <- joined_data %>% filter(complete.cases(joined_data[, covs_to_check]))
  
  cat("✓ Merged and calibrated data successfully. N =", nrow(joined_data), "clusters.\n")
  
  # --- 2. Define the 37 Candidate Covariate Models (Excluding Basin & MegaHx) ----
  formulas_list <- list(
    # --- Base Models ---
    "M2.1: UOI Only"                           = B_H_index ~ uoi,
    "M2.2p: Elephant Possible Only"            = B_H_index ~ elephant_present_possible,
    "M2.2s: Elephant Strict Only"              = B_H_index ~ elephant_present_strict,
    
    # --- UOI + Environmental Covariate Set ---
    "M2.3: UOI + Elev"                         = B_H_index ~ uoi + elevation,
    "M2.4: UOI + Slope"                        = B_H_index ~ uoi + slope,
    "M2.5: UOI + HAND"                         = B_H_index ~ uoi + hnd,
    "M2.6: UOI + Precip"                       = B_H_index ~ uoi + precip,
    "M2.7: UOI + Clay"                         = B_H_index ~ uoi + clay,
    "M2.8: UOI + Forest"                       = B_H_index ~ uoi + forest_fraction,
 
    # --- Main Effect Backbone Set (Elephant Possible) ---
    "M2.9p: UOI + Elephant Possible"           = B_H_index ~ uoi + elephant_present_possible,
    "M2.10p: UOI + Elephant Possible + Elev"   = B_H_index ~ uoi + elephant_present_possible + elevation,
    "M2.11p: UOI + Elephant Possible + Slope"  = B_H_index ~ uoi + elephant_present_possible + slope,
    "M2.12p: UOI + Elephant Possible + HAND"   = B_H_index ~ uoi + elephant_present_possible + hnd,
    "M2.13p: UOI + Elephant Possible + Precip" = B_H_index ~ uoi + elephant_present_possible + precip,
    "M2.14p: UOI + Elephant Possible + Clay"   = B_H_index ~ uoi + elephant_present_possible + clay,
    "M2.15p: UOI + Elephant Possible + Forest" = B_H_index ~ uoi + elephant_present_possible + forest_fraction,
    
    # --- Main Effect Backbone Set (Elephant Strict) ---
    "M2.9s: UOI + Elephant Strict"             = B_H_index ~ uoi + elephant_present_strict,
    "M2.10s: UOI + Elephant Strict + Elev"     = B_H_index ~ uoi + elephant_present_strict + elevation,
    "M2.11s: UOI + Elephant Strict + Slope"     = B_H_index ~ uoi + elephant_present_strict + slope,
    "M2.12s: UOI + Elephant Strict + HAND"      = B_H_index ~ uoi + elephant_present_strict + hnd,
    "M2.13s: UOI + Elephant Strict + Precip"    = B_H_index ~ uoi + elephant_present_strict + precip,
    "M2.14s: UOI + Elephant Strict + Clay"      = B_H_index ~ uoi + elephant_present_strict + clay,
    "M2.15s: UOI + Elephant Strict + Forest"    = B_H_index ~ uoi + elephant_present_strict + forest_fraction,
    
    # --- Interaction Effect Backbone Set (Elephant Possible) ---
    "M2.16p: UOI * Elephant Possible"           = B_H_index ~ uoi * elephant_present_possible,
    "M2.17p: UOI * Elephant Possible + Elev"   = B_H_index ~ uoi * elephant_present_possible + elevation,
    "M2.18p: UOI * Elephant Possible + Slope"  = B_H_index ~ uoi * elephant_present_possible + slope,
    "M2.19p: UOI * Elephant Possible + HAND"   = B_H_index ~ uoi * elephant_present_possible + hnd,
    "M2.20p: UOI * Elephant Possible + Precip" = B_H_index ~ uoi * elephant_present_possible + precip,
    "M2.21p: UOI * Elephant Possible + Clay"   = B_H_index ~ uoi * elephant_present_possible + clay,
    "M2.22p: UOI * Elephant Possible + Forest" = B_H_index ~ uoi * elephant_present_possible + forest_fraction,
    
    # --- Interaction Effect Backbone Set (Elephant Strict) ---
    "M2.16s: UOI * Elephant Strict"             = B_H_index ~ uoi * elephant_present_strict,
    "M2.17s: UOI * Elephant Strict + Elev"     = B_H_index ~ uoi * elephant_present_strict + elevation,
    "M2.18s: UOI * Elephant Strict + Slope"     = B_H_index ~ uoi * elephant_present_strict + slope,
    "M2.19s: UOI * Elephant Strict + HAND"      = B_H_index ~ uoi * elephant_present_strict + hnd,
    "M2.20s: UOI * Elephant Strict + Precip"    = B_H_index ~ uoi * elephant_present_strict + precip,
    "M2.21s: UOI * Elephant Strict + Clay"      = B_H_index ~ uoi * elephant_present_strict + clay,
    "M2.22s: UOI * Elephant Strict + Forest"    = B_H_index ~ uoi * elephant_present_strict + forest_fraction
  )
  
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
    min_edf <- min(full_edf[compliant_indices])
    best_parsimonious_indices <- compliant_indices[which(full_edf[compliant_indices] == min_edf)]
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
    edf <- full_edf[m_idx]
    
    if (m_idx == selected_idx) {
      return("Parsimonious Selected Best")
    } else if (m_idx == raw_best_idx) {
      return("Raw Best (More Complex)")
    } else if (mae <= threshold_1se) {
      if (edf < full_edf[raw_best_idx]) {
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
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
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
  
  # --- 7. Generate Figure 4 Panels ---------------------------------------------
  cat("Generating PNAS-styled Figure 4 panels (LOBO and AIC variants)...\n")
  source("code/functions/theme_pnas.R")
  
  # A highly DRY and robust function to compile PNAS Figure 4 for any model
  generate_figure4_trio <- function(model_name, model_obj, results_df, joined_data, full_models, is_aic = FALSE) {
    # pal_region and pal_elephant_binary are defined in theme_pnas.R (sourced above)
    uoi_seq <- seq(from = 0.918, to = 0.970, length.out = 300)
    covs_to_fill <- c("elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
    
    best_formula_vars <- all.vars(formula(model_obj))
    uses_elephant_possible <- "elephant_present_possible" %in% best_formula_vars
    uses_elephant_strict   <- "elephant_present_strict" %in% best_formula_vars
    uses_elephant <- uses_elephant_possible || uses_elephant_strict || ("elephant_present" %in% best_formula_vars)
    uses_basin <- "basin" %in% best_formula_vars
    
    ele_col <- if (uses_elephant_strict) {
      "elephant_present_strict"
    } else if (uses_elephant_possible) {
      "elephant_present_possible"
    } else {
      "elephant_present"
    }
    
    # --- Panel A: Standing Mammal Biomass vs. GEDI Openness ---
    if (uses_basin) {
      pred_list <- lapply(levels(joined_data$basin), function(b) {
        nd <- data.frame(uoi = uoi_seq, basin = factor(b, levels = levels(joined_data$basin)))
        for (cv in covs_to_fill) nd[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
        nd$fit <- predict(model_obj, newdata = nd, type = "response")
        nd
      })
      pred_plot <- do.call(rbind, pred_list)
      
      p_a <- ggplot() +
        geom_point(data = joined_data, aes(x = uoi, y = B_H_index, fill = basin, size = trap_days, alpha = w_temp_cluster, shape = basin),
                   color = "black", stroke = 0.3) +
        geom_line(data = pred_plot, aes(x = uoi, y = fit, color = basin), linewidth = 0.75) +
        
        scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22), name = "Region") +
        scale_color_manual(values = pal_region, name = "Region") +
        scale_size_continuous(name = "Effort (Trap-days)", range = c(1.2, 4.0), breaks = c(100, 1000, 5000, 15000)) +
        scale_alpha_continuous(name = "Temporal Alignment Weight", range = c(0.25, 1.0), breaks = c(0.1, 0.5, 1.0), labels = c("Historical", "Intermediate", "Contemp.")) +
        scale_fill_manual(values = pal_region, name = "Region") +
        scale_x_continuous(breaks = seq(0.92, 0.97, by = 0.01), limits = c(0.918, 0.970)) +
        scale_y_continuous(trans = "log1p", labels = comma_format(), breaks = c(0, 10, 100, 1000, 3000), limits = c(0, 5000)) +
        labs(
          title = "A. Standing Mammal Biomass vs. GEDI Openness",
          x = "GEDI Understory Openness Index (UOI)",
          y = "Total Mammal Biomass Index (log1p scale)"
        )
    } else if (uses_elephant) {
      pred_df_absent <- data.frame(uoi = uoi_seq)
      pred_df_absent[[ele_col]] <- factor("Absent", levels = c("Absent", "Present"))
      for (cv in covs_to_fill) pred_df_absent[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
      pred_df_absent$fit <- predict(model_obj, newdata = pred_df_absent, type = "response")
      
      pred_df_present <- data.frame(uoi = uoi_seq)
      pred_df_present[[ele_col]] <- factor("Present", levels = c("Absent", "Present"))
      for (cv in covs_to_fill) pred_df_present[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
      pred_df_present$fit <- predict(model_obj, newdata = pred_df_present, type = "response")
      
      pred_total_plot <- rbind(pred_df_absent, pred_df_present)
      
      p_a <- ggplot() +
        geom_point(data = joined_data, aes(x = uoi, y = B_H_index, fill = .data[[ele_col]], size = trap_days, alpha = w_temp_cluster, shape = basin),
                   color = "black", stroke = 0.3) +
        geom_line(data = pred_total_plot, aes(x = uoi, y = fit, color = .data[[ele_col]]), linewidth = 0.75) +
        
        scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22), name = "Region") +
        scale_color_manual(values = pal_elephant_binary, name = "Elephant:") +
        scale_size_continuous(name = "Effort (Trap-days)", range = c(1.2, 4.0), breaks = c(100, 1000, 5000, 15000)) +
        scale_alpha_continuous(name = "Temporal Alignment Weight", range = c(0.25, 1.0), breaks = c(0.1, 0.5, 1.0), labels = c("Historical", "Intermediate", "Contemp.")) +
        scale_fill_manual(values = pal_elephant_binary, name = "Elephant:") +
        scale_x_continuous(breaks = seq(0.92, 0.97, by = 0.01), limits = c(0.918, 0.970)) +
        scale_y_continuous(trans = "log1p", labels = comma_format(), breaks = c(0, 10, 100, 1000, 3000), limits = c(0, 5000)) +
        labs(
          title = "A. Standing Mammal Biomass vs. GEDI Openness",
          x = "GEDI Understory Openness Index (UOI)",
          y = "Total Mammal Biomass Index (log1p scale)"
        )
    } else {
      pred_df <- data.frame(uoi = uoi_seq)
      for (cv in covs_to_fill) pred_df[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
      pred_df$fit <- predict(model_obj, newdata = pred_df, type = "response")
      
      p_a <- ggplot() +
        geom_point(data = joined_data, aes(x = uoi, y = B_H_index, fill = elephant_present_strict, size = trap_days, alpha = w_temp_cluster, shape = basin),
                   color = "black", stroke = 0.3) +
        geom_line(data = pred_df, aes(x = uoi, y = fit), color = pal_elephant_binary[["Present"]], linewidth = 0.75) +
        
        scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22), name = "Region") +
        scale_size_continuous(name = "Effort (Trap-days)", range = c(1.2, 4.0), breaks = c(100, 1000, 5000, 15000)) +
        scale_alpha_continuous(name = "Temporal Alignment Weight", range = c(0.25, 1.0), breaks = c(0.1, 0.5, 1.0), labels = c("Historical", "Intermediate", "Contemp.")) +
        scale_fill_manual(values = pal_elephant_binary, name = "Elephant:") +
        scale_x_continuous(breaks = seq(0.92, 0.97, by = 0.01), limits = c(0.918, 0.970)) +
        scale_y_continuous(trans = "log1p", labels = comma_format(), breaks = c(0, 10, 100, 1000, 3000), limits = c(0, 5000)) +
        labs(
          title = "A. Standing Mammal Biomass vs. GEDI Openness",
          x = "GEDI Understory Openness Index (UOI)",
          y = "Total Mammal Biomass Index (log1p scale)"
        )
    }
    
    p_a <- p_a +
      theme_pnas(base_size = 7.5) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 8.5, margin = margin(b = 6, t = 4)),
        axis.title.x = element_text(margin = margin(t = 4)),
        axis.title.y = element_text(margin = margin(r = 4)),
        plot.margin = margin(t = 6, r = 4, b = 6, l = 4, unit = "pt")
      )
    
    # --- Panel B: Residual stability and temporal independence ---
    joined_data$residuals <- log1p(joined_data$B_H_index) - log1p(fitted(model_obj))
    
    fill_col <- if (uses_basin) {
      "basin"
    } else if (uses_elephant) {
      ele_col
    } else {
      "elephant_present_strict"
    }
    
    p_b <- ggplot(joined_data, aes(x = w_temp_cluster, y = residuals)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "#555555", linewidth = 0.4) +
      geom_point(aes(fill = .data[[fill_col]], size = trap_days, alpha = w_temp_cluster, shape = basin), color = "black", stroke = 0.3) +
      geom_smooth(method = "lm", aes(weight = w_combined_norm), formula = y ~ x, color = pal_elephant_binary[["Present"]], linewidth = 0.6, se = TRUE, alpha = 0.1) +
      
      scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22), name = "Region") +
      scale_size_continuous(name = "Effort (Trap-days)", range = c(1.2, 4.0), breaks = c(100, 1000, 5000, 15000)) +
      scale_alpha_continuous(name = "Temporal Alignment Weight", range = c(0.25, 1.0), breaks = c(0.1, 0.5, 1.0), labels = c("Historical", "Intermediate", "Contemporaneous")) +
      scale_x_continuous(breaks = seq(0.1, 1.0, by = 0.2), limits = c(0.08, 1.02)) +
      scale_y_continuous(breaks = seq(-4, 4, by = 2), limits = c(-4.5, 4.5)) +
      labs(
        title = "B. Residual Independence & Temporal Stability",
        x = "Cluster Temporal Alignment Weight (W_temp)",
        y = "Best Model log1p Residuals"
      ) +
      theme_pnas(base_size = 7.5) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 8.5, margin = margin(b = 6, t = 4)),
        axis.title.x = element_text(margin = margin(t = 4)),
        axis.title.y = element_text(margin = margin(r = 4)),
        plot.margin = margin(t = 6, r = 4, b = 6, l = 4, unit = "pt")
      )
    
    if (fill_col == "basin") {
      p_b <- p_b + scale_fill_manual(values = pal_region, name = "Region:")
    } else {
      p_b <- p_b + scale_fill_manual(values = pal_elephant_binary, name = "Elephant:")
    }
    
    # --- Panel C: Model Selection Bar Plot ---
    format_model_label_lobo <- function(model_name, model_obj) {
      clean_name <- sub("^M2\\.[0-9\\.]+[a-z_]*: ", "", model_name)
      tokens <- strsplit(clean_name, "\\s+")[[1]]
      tokens <- tokens[tokens != ""]
      p_table <- summary(model_obj)$p.table
      
      var_map <- list(
        "ElephantPossible" = "elephant_present_possiblePresent",
        "ElephantStrict"   = "elephant_present_strictPresent",
        "Elephant"         = c("elephant_presentPresent", "elephant_present_possiblePresent", "elephant_present_strictPresent"),
        "Basin"            = c("basinCongo", "basinSE_Asia"),
        "UOI"              = "uoi",
        "Elevation"        = "elevation",
        "Elev"             = "elevation",
        "Slope"            = "slope",
        "HAND"             = "hnd",
        "Precipitation"    = "precip",
        "Precip"           = "precip",
        "Clay"             = "clay",
        "Forest"           = "forest_fraction",
        "UOI:Elephant"     = c("uoi:elephant_presentPresent", "uoi:elephant_present_possiblePresent", "uoi:elephant_present_strictPresent"),
        "UOI:Basin"        = c("uoi:basinCongo", "uoi:basinSE_Asia")
      )
      
      plotmath_tokens <- sapply(tokens, function(tok) {
        if (tok %in% c("+", "*", ":")) return(sprintf("plain(\" %s \")", tok))
        if (tok == "Only") return(sprintf("plain(\" %s\")", tok))
        
        matched_terms <- var_map[[tok]]
        if (!is.null(matched_terms)) {
          is_sig <- FALSE
          for (term in matched_terms) {
            if (term %in% rownames(p_table)) {
              p_val <- p_table[term, ncol(p_table)]
              if (!is.na(p_val) && p_val < 0.05) {
                is_sig <- TRUE
                break
              }
            }
          }
          return(ifelse(is_sig, sprintf("bold(\"%s\")", tok), sprintf("plain(\"%s\")", tok)))
        } else {
          return(sprintf("plain(\"%s\")", tok))
        }
      })
      
      paste(plotmath_tokens, collapse = " * ")
    }
    
    if (is_aic) {
      plot_df <- results_df %>%
        arrange(Full_AICc) %>%
        head(15)
    } else {
      plot_df <- results_df %>%
        head(15)
    }
    
    plot_df$plotmath_label <- sapply(1:nrow(plot_df), function(i) {
      format_model_label_lobo(plot_df$Model[i], full_models[[plot_df$Model[i]]])
    })
    
    plot_sel_df <- plot_df %>%
      mutate(
        CleanName = sub("^M2\\.[0-9\\.]+[a-z_]*: ", "", Model),
        CleanName = factor(CleanName, levels = rev(CleanName))
      )
    
    ordered_exprs <- plot_sel_df$plotmath_label[match(levels(plot_sel_df$CleanName), plot_sel_df$CleanName)]
    parsed_labels <- parse(text = ordered_exprs)
    
    plot_sel_df$Full_DevExpl_Pct <- plot_sel_df$Full_DevExpl * 100
    
    if (is_aic) {
      # AIC-based plot matches Figure 3's structure:
      # x-axis is Deviance Explained (goodness of fit, matching adjusted Pseudo-R2 in F1)
      # fill is delta_AIC (Delta AIC, matching Delta AIC fill in F1)
      p_c <- plot_model_selection_bars(
        plot_df = plot_sel_df,
        x_var = "Full_DevExpl_Pct",
        fill_var = "delta_AICc",
        fill_label = "Delta AICc",
        x_label = "Model Deviance Explained (%)",
        plot_title = sprintf("C. Standard AICc Model Selection (All %d Candidates)", nrow(plot_sel_df)),
        parsed_labels = parsed_labels
      )
    } else {
      # LOBO-based plot:
      # x-axis is Deviance Explained (goodness of fit)
      # fill is out-of-sample prediction error (OOS MAE)
      p_c <- plot_model_selection_bars(
        plot_df = plot_sel_df,
        x_var = "Full_DevExpl_Pct",
        fill_var = "OOS_MAE_log",
        fill_label = "OOS MAE (log1p)",
        x_label = "Model Deviance Explained (%)",
        plot_title = "C. Generalizability Model Selection (LOBO Out-of-Sample Validation)",
        parsed_labels = parsed_labels
      )
    }
    
    # --- 8. Construct Clean Shared Legend ---
    p_legend_obj <- ggplot(joined_data) +
      geom_point(aes(x = uoi, y = B_H_index, fill = elephant_present_strict, size = trap_days, alpha = w_temp_cluster, shape = basin), color = "black", stroke = 0.3) +
      scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22), name = "Region:") +
      scale_size_continuous(name = "Effort (Trap-days):", breaks = c(100, 1000, 5000, 15000), range = c(1.2, 4.0)) +
      scale_alpha_continuous(name = "Temporal Alignment:", range = c(0.25, 1.0), breaks = c(0.1, 0.5, 1.0), labels = c("Hist.", "Interm.", "Contemp.")) +
      scale_fill_manual(values = pal_elephant_binary, name = "Elephant:") +
      theme_pnas(base_size = 7.5) +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 7.0, face = "bold"),
        legend.text = element_text(size = 6.5)
      )
    
    shared_legend <- cowplot::get_legend(p_legend_obj)
    
    # --- 9. Assemble and Save Multipanel Figure ---
    row1 <- cowplot::plot_grid(
      p_a, p_b,
      ncol = 2,
      align = "h",
      axis = "tb",
      rel_widths = c(1.0, 1.0)
    )
    
    fig_final <- cowplot::plot_grid(
      row1,
      p_c,
      shared_legend,
      ncol = 1,
      rel_heights = c(1.0, 0.9, 0.12)
    )
    
    return(fig_final)
  }
  
  # 1. Generate and save the standard AICc-selected Figure 4 (alternate plot)
  cat("Generating Standard AICc-selected Figure 4...\n")
  fig_alt <- generate_figure4_trio(
    model_name = best_model_name_aic,
    model_obj = best_model_aic,
    results_df = results_df,
    joined_data = joined_data,
    full_models = full_models,
    is_aic = TRUE
  )
  
  fig4_png_path <- file.path(figures_dir, "figure4.png")
  save_pnas(plot = fig_alt, filename = fig4_png_path, type = "double", height_cm = 12.5)
  fig4_pdf_path <- sub("\\.png$", ".pdf", fig4_png_path)
  save_pnas(plot = fig_alt, filename = fig4_pdf_path, type = "double", height_cm = 12.5)
  
  # 2. Generate and save the LORO-CV/LOBO-selected Figure S2 (main plot)
  cat("Generating Generalizability LORO-CV-selected Figure S2...\n")
  fig_main <- generate_figure4_trio(
    model_name = best_model_name,
    model_obj = best_model,
    results_df = results_df,
    joined_data = joined_data,
    full_models = full_models,
    is_aic = FALSE
  )
  
  figS2_png_path <- file.path(figures_dir, "figureS2.png")
  save_pnas(plot = fig_main, filename = figS2_png_path, type = "double", height_cm = 12.5)
  figS2_pdf_path <- sub("\\.png$", ".pdf", figS2_png_path)
  save_pnas(plot = fig_main, filename = figS2_pdf_path, type = "double", height_cm = 12.5)
  

  
  cat(sprintf("✓ Saved Figure 4 to %s and Figure S2 to %s\n", fig4_png_path, figS2_png_path))
  cat("=== Framework 2 Integrated Analysis Completed Successfully ===\n")
  
  return(list(
    results_df = results_df,
    best_model = best_model
  ))
}

# --- Execute directly if called from terminal ---
run_framework2_analysis()
