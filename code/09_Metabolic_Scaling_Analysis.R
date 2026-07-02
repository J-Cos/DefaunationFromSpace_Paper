# =============================================================================
# code/09_Metabolic_Scaling_Analysis.R
#
# A rigorous, comparative, and reproducible scientific analysis evaluating the 
# impact of Metabolic Scaling Theory (MST, scaling exponent beta = 0.75) versus
# static standing biomass (beta = 1.0) on model performance, environmental 
# covariate selection, and generalizability across tropical forest basins.
#
# Inputs:
#   - outputs/camera_traps_joint_detections.csv
#   - outputs/camera_traps_robust_buffered_mcps.geojson
#   - outputs/EOdata/analysis_stack_5000_{Basin}.tif
#
# Outputs:
#   - outputs/metabolic_scaling_model_selection.csv
# =============================================================================

library(mgcv)
library(dplyr)
library(readr)

# --- 1. Load Custom Calibration Helpers & Core Data --------------------------
source("code/functions/theme_pnas.R")
source("code/functions/calibration_helpers.R")
source("code/functions/model_convergence.R")

cat("=====================================================================\n")
cat("=== code/09_Metabolic_Scaling_Analysis.R                          ===\n")
cat("=== Rigorous Comparative Evaluation of Metabolic Scaling Theory   ===\n")
cat("=====================================================================\n\n")

# Load unified cluster-level data at 5 km resolution
joined_data <- extract_scale_data(5000)
covs_to_check <- c("uoi", "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
joined_data <- joined_data %>% filter(complete.cases(joined_data[, covs_to_check]))

cat(sprintf("✓ Loaded %d cluster-level records at the 5 km spatial scale.\n\n", nrow(joined_data)))

# Initialize a results list to export as a CSV
comparison_records <- list()

# -----------------------------------------------------------------------------
# 2. Framework 1: Beta Regressions (Predicting UOI)
# -----------------------------------------------------------------------------
cat("--- [Framework 1] Beta Regressions Predicting Understory Openness (UOI) ---\n")

f1_templates <- list(
  "Index Only"                       = "uoi ~ {Index}",
  "Index + Elev"                     = "uoi ~ {Index} + elevation",
  "Index + Slope"                    = "uoi ~ {Index} + slope",
  "Index + HAND"                     = "uoi ~ {Index} + hnd",
  "Index + Precip"                   = "uoi ~ {Index} + precip",
  "Index + Clay"                     = "uoi ~ {Index} + clay",
  "Index + Forest"                   = "uoi ~ {Index} + forest_fraction",
  "Index + ElephantPossible"         = "uoi ~ {Index} + elephant_present_possible",
  "Index + ElephantPossible + Elev"  = "uoi ~ {Index} + elephant_present_possible + elevation",
  "Index + ElephantPossible + Slope" = "uoi ~ {Index} + elephant_present_possible + slope",
  "Index + ElephantPossible + HAND"  = "uoi ~ {Index} + elephant_present_possible + hnd",
  "Index + ElephantPossible + Precip"= "uoi ~ {Index} + elephant_present_possible + precip",
  "Index + ElephantPossible + Clay"  = "uoi ~ {Index} + elephant_present_possible + clay",
  "Index + ElephantPossible + Forest"= "uoi ~ {Index} + elephant_present_possible + forest_fraction",
  "Index * ElephantPossible"         = "uoi ~ {Index} * elephant_present_possible",
  "Index * ElephantPossible + Elev"  = "uoi ~ {Index} * elephant_present_possible + elevation",
  "Index * ElephantPossible + Slope" = "uoi ~ {Index} * elephant_present_possible + slope",
  "Index * ElephantPossible + HAND"  = "uoi ~ {Index} * elephant_present_possible + hnd",
  "Index * ElephantPossible + Precip"= "uoi ~ {Index} * elephant_present_possible + precip",
  "Index * ElephantPossible + Clay"  = "uoi ~ {Index} * elephant_present_possible + clay",
  "Index * ElephantPossible + Forest"= "uoi ~ {Index} * elephant_present_possible + forest_fraction",
  "Index + ElephantStrict"           = "uoi ~ {Index} + elephant_present_strict",
  "Index + ElephantStrict + Elev"     = "uoi ~ {Index} + elephant_present_strict + elevation",
  "Index + ElephantStrict + Slope"     = "uoi ~ {Index} + elephant_present_strict + slope",
  "Index + ElephantStrict + HAND"      = "uoi ~ {Index} + elephant_present_strict + hnd",
  "Index + ElephantStrict + Precip"    = "uoi ~ {Index} + elephant_present_strict + precip",
  "Index + ElephantStrict + Clay"      = "uoi ~ {Index} + elephant_present_strict + clay",
  "Index + ElephantStrict + Forest"    = "uoi ~ {Index} + elephant_present_strict + forest_fraction",
  "Index * ElephantStrict"             = "uoi ~ {Index} * elephant_present_strict",
  "Index * ElephantStrict + Elev"     = "uoi ~ {Index} * elephant_present_strict + elevation",
  "Index * ElephantStrict + Slope"     = "uoi ~ {Index} * elephant_present_strict + slope",
  "Index * ElephantStrict + HAND"      = "uoi ~ {Index} * elephant_present_strict + hnd",
  "Index * ElephantStrict + Precip"    = "uoi ~ {Index} * elephant_present_strict + precip",
  "Index * ElephantStrict + Clay"      = "uoi ~ {Index} * elephant_present_strict + clay",
  "Index * ElephantStrict + Forest"    = "uoi ~ {Index} * elephant_present_strict + forest_fraction"
)

run_f1_selection <- function(index_name) {
  best_aicc <- Inf
  best_name <- ""
  best_formula <- NULL
  n_obs <- nrow(joined_data)
  
  for (n in names(f1_templates)) {
    f_str <- gsub("\\{Index\\}", index_name, f1_templates[[n]])
    f <- as.formula(f_str)
    
    fit <- tryCatch({
      m <- gam(f, data = joined_data, family = betar(link = "logit"), weights = w_combined_norm, method = "ML")
      check_model_convergence(m, paste("F1", index_name, n))
      m
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      aic_val <- AIC(fit)
      k <- attr(logLik(fit), "df")
      aicc_val <- aic_val + (2 * k * (k + 1)) / (n_obs - k - 1)
      dev_expl <- summary(fit)$dev.expl
      
      # Record all model fits for complete tracking
      comparison_records[[length(comparison_records) + 1]] <<- list(
        Framework = "Framework 1",
        Metric = "AICc & DevExpl",
        ModelLabel = n,
        IndexUsed = index_name,
        Formula = f_str,
        FullAICc = aicc_val,
        DevianceExplained = dev_expl,
        OOS_MAE = NA
      )
      
      if (aicc_val < best_aicc) {
        best_aicc <- aicc_val
        best_name <- n
        best_formula <- f_str
      }
    }
  }
  return(list(name = best_name, formula = best_formula, aicc = best_aicc))
}

sel_f1_biomass <- run_f1_selection("B_H_index")
sel_f1_metabolism <- run_f1_selection("M_H_index")

cat(sprintf("  ★ Best Biomass Model:    \"%s\" (AICc = %.2f)\n  Formula: %s\n\n", 
            sel_f1_biomass$name, sel_f1_biomass$aicc, sel_f1_biomass$formula))
cat(sprintf("  ★ Best Metabolism Model: \"%s\" (AICc = %.2f)\n  Formula: %s\n\n", 
            sel_f1_metabolism$name, sel_f1_metabolism$aicc, sel_f1_metabolism$formula))

# -----------------------------------------------------------------------------
# 3. Framework 2: Tweedie GLMs
# -----------------------------------------------------------------------------
cat("--- [Framework 2] Tweedie GLMs (Biomass vs. Metabolism) ---\n")

f2_templates <- list(
  "UOI Only"                           = "{Index} ~ uoi",
  "UOI + Elev"                         = "{Index} ~ uoi + elevation",
  "UOI + Slope"                        = "{Index} ~ uoi + slope",
  "UOI + HAND"                         = "{Index} ~ uoi + hnd",
  "UOI + Precip"                       = "{Index} ~ uoi + precip",
  "UOI + Clay"                         = "{Index} ~ uoi + clay",
  "UOI + Forest"                       = "{Index} ~ uoi + forest_fraction",
  "UOI * Basin"                        = "{Index} ~ uoi * basin",
  "UOI * Basin + Elev"                 = "{Index} ~ uoi * basin + elevation",
  "UOI * Basin + Slope"                = "{Index} ~ uoi * basin + slope",
  "UOI * Basin + HAND"                 = "{Index} ~ uoi * basin + hnd",
  "UOI * Basin + Precip"               = "{Index} ~ uoi * basin + precip",
  "UOI * Basin + Clay"                 = "{Index} ~ uoi * basin + clay",
  "UOI * Basin + Forest"               = "{Index} ~ uoi * basin + forest_fraction",
  "UOI * ElephantPossible"             = "{Index} ~ uoi * elephant_present_possible",
  "UOI * ElephantPossible + Elev"      = "{Index} ~ uoi * elephant_present_possible + elevation",
  "UOI * ElephantPossible + Slope"      = "{Index} ~ uoi * elephant_present_possible + slope",
  "UOI * ElephantPossible + HAND"       = "{Index} ~ uoi * elephant_present_possible + hnd",
  "UOI * ElephantPossible + Precip"     = "{Index} ~ uoi * elephant_present_possible + precip",
  "UOI * ElephantPossible + Clay"       = "{Index} ~ uoi * elephant_present_possible + clay",
  "UOI * ElephantPossible + Forest"     = "{Index} ~ uoi * elephant_present_possible + forest_fraction",
  "UOI * ElephantStrict"               = "{Index} ~ uoi * elephant_present_strict",
  "UOI * ElephantStrict + Elev"         = "{Index} ~ uoi * elephant_present_strict + elevation",
  "UOI * ElephantStrict + Slope"         = "{Index} ~ uoi * elephant_present_strict + slope",
  "UOI * ElephantStrict + HAND"          = "{Index} ~ uoi * elephant_present_strict + hnd",
  "UOI * ElephantStrict + Precip"        = "{Index} ~ uoi * elephant_present_strict + precip",
  "UOI * ElephantStrict + Clay"          = "{Index} ~ uoi * elephant_present_strict + clay",
  "UOI * ElephantStrict + Forest"        = "{Index} ~ uoi * elephant_present_strict + forest_fraction"
)

# A. Standard AICc Selection Pathway (Full Sample Fit)
run_f2_aic_selection <- function(index_name) {
  best_aicc <- Inf
  best_name <- ""
  best_formula <- NULL
  n_obs <- nrow(joined_data)
  
  for (n in names(f2_templates)) {
    f_str <- gsub("\\{Index\\}", index_name, f2_templates[[n]])
    f <- as.formula(f_str)
    
    fit <- tryCatch({
      m <- gam(f, data = joined_data, family = tw(), weights = w_combined_norm, method = "ML")
      check_model_convergence(m, paste("F2 AICc", index_name, n))
      m
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      aic_val <- AIC(fit)
      k <- attr(logLik(fit), "df")
      aicc_val <- aic_val + (2 * k * (k + 1)) / (n_obs - k - 1)
      dev_expl <- summary(fit)$dev.expl
      
      comparison_records[[length(comparison_records) + 1]] <<- list(
        Framework = "Framework 2 AICc",
        Metric = "AICc & DevExpl",
        ModelLabel = n,
        IndexUsed = index_name,
        Formula = f_str,
        FullAICc = aicc_val,
        DevianceExplained = dev_expl,
        OOS_MAE = NA
      )
      
      if (aicc_val < best_aicc) {
        best_aicc <- aicc_val
        best_name <- n
        best_formula <- f_str
      }
    }
  }
  return(list(name = best_name, formula = best_formula, aicc = best_aicc))
}

sel_f2_aic_biomass <- run_f2_aic_selection("B_H_index")
sel_f2_aic_metabolism <- run_f2_aic_selection("M_H_index")

cat("A. Standard AICc Selection Pathway (Full Sample Fit):\n")
cat(sprintf("  ★ Best Biomass Model:    \"%s\" (AICc = %.2f)\n  Formula: %s\n\n", 
            sel_f2_aic_biomass$name, sel_f2_aic_biomass$aicc, sel_f2_aic_biomass$formula))
cat(sprintf("  ★ Best Metabolism Model: \"%s\" (AICc = %.2f)\n  Formula: %s\n\n", 
            sel_f2_aic_metabolism$name, sel_f2_aic_metabolism$aicc, sel_f2_aic_metabolism$formula))

# B. LORO-CV Generalizability Selection Pathway (Out-of-Sample CV MAE)
run_f2_lobo_selection <- function(index_name) {
  basins <- unique(joined_data$basin)
  
  valid_templates <- list()
  for (n in names(f2_templates)) {
    f_str <- gsub("\\{Index\\}", index_name, f2_templates[[n]])
    if (grepl("basin", f_str)) next # Skip basin terms due to LOBO folding
    valid_templates[[n]] <- f_str
  }
  
  num_models <- length(valid_templates)
  model_names <- names(valid_templates)
  
  oos_predictions <- matrix(NA, nrow = nrow(joined_data), ncol = num_models)
  colnames(oos_predictions) <- model_names
  
  full_edf <- numeric(num_models)
  full_aicc <- numeric(num_models)
  full_devexpl <- numeric(num_models)
  
  for (m_idx in 1:num_models) {
    n <- model_names[m_idx]
    f_str <- valid_templates[[n]]
    f <- as.formula(f_str)
    
    # 1. LOBO cross-validation
    for (b in basins) {
      train_idx <- which(joined_data$basin != b)
      test_idx  <- which(joined_data$basin == b)
      
      if (length(train_idx) == 0 || length(test_idx) == 0) next
      
      fit_fold <- tryCatch({
        m_fold <- gam(f, data = joined_data[train_idx, ], family = tw(), weights = w_combined_norm, method = "ML")
        check_model_convergence(m_fold, paste("F2 LORO-CV", index_name, n, "fold", b), raise_warning = FALSE)
        m_fold
      }, error = function(e) NULL)
      
      if (!is.null(fit_fold)) {
        pred <- predict(fit_fold, newdata = joined_data[test_idx, ], type = "response")
        oos_predictions[test_idx, m_idx] <- pred
      }
    }
    
    # 2. Fit full model to get EDF, AICc and Deviance Explained
    fit_full <- tryCatch({
      m_full <- gam(f, data = joined_data, family = tw(), weights = w_combined_norm, method = "ML")
      check_model_convergence(m_full, paste("F2 LORO-CV Full", index_name, n))
      m_full
    }, error = function(e) NULL)
    
    if (!is.null(fit_full)) {
      full_edf[m_idx] <- sum(fit_full$edf)
      full_devexpl[m_idx] <- summary(fit_full)$dev.expl
      
      aic_val <- AIC(fit_full)
      k <- attr(logLik(fit_full), "df")
      full_aicc[m_idx] <- aic_val + (2 * k * (k + 1)) / (nrow(joined_data) - k - 1)
    } else {
      full_edf[m_idx] <- NA
      full_devexpl[m_idx] <- NA
      full_aicc[m_idx] <- NA
    }
  }
  
  # 3. Calculate micro-averaged OOS MAE and SE
  log_y_obs <- log1p(joined_data[[index_name]])
  oos_MAE_log <- sapply(1:num_models, function(m_idx) {
    abs_err <- abs(log_y_obs - log1p(oos_predictions[, m_idx]))
    if (all(is.na(abs_err))) return(NA)
    mean(abs_err, na.rm = TRUE)
  })
  
  oos_MAE_log_SE <- sapply(1:num_models, function(m_idx) {
    abs_err <- abs(log_y_obs - log1p(oos_predictions[, m_idx]))
    if (all(is.na(abs_err))) return(NA)
    sd(abs_err, na.rm = TRUE) / sqrt(sum(!is.na(abs_err)))
  })
  
  # 4. Apply 1-SE parsimony rule
  valid_idx <- which(!is.na(oos_MAE_log))
  if (length(valid_idx) > 0) {
    raw_best_idx <- valid_idx[which.min(oos_MAE_log[valid_idx])]
    best_mae <- oos_MAE_log[raw_best_idx]
    best_se <- oos_MAE_log_SE[raw_best_idx]
    threshold_1se <- best_mae + best_se
    
    compliant_indices <- valid_idx[which(oos_MAE_log[valid_idx] <= threshold_1se)]
    min_edf <- min(full_edf[compliant_indices], na.rm = TRUE)
    best_parsimonious_indices <- compliant_indices[which(full_edf[compliant_indices] == min_edf)]
    selected_idx <- best_parsimonious_indices[which.min(oos_MAE_log[best_parsimonious_indices])]
    
    best_name <- model_names[selected_idx]
    best_mae_val <- oos_MAE_log[selected_idx]
  } else {
    best_name <- model_names[1]
    best_mae_val <- NA
  }
  
  # Record the comparison logs
  for (m_idx in 1:num_models) {
    n <- model_names[m_idx]
    comparison_records[[length(comparison_records) + 1]] <<- list(
      Framework = "Framework 2 LORO-CV",
      Metric = "LORO-CV MAE",
      ModelLabel = n,
      IndexUsed = index_name,
      Formula = valid_templates[[n]],
      FullAICc = full_aicc[m_idx],
      DevianceExplained = full_devexpl[m_idx],
      OOS_MAE = oos_MAE_log[m_idx]
    )
  }
  
  best_formula <- gsub("\\{Index\\}", index_name, f2_templates[[best_name]])
  return(list(name = best_name, formula = best_formula, mae = best_mae_val))
}

sel_f2_lobo_biomass <- run_f2_lobo_selection("B_H_index")
sel_f2_lobo_metabolism <- run_f2_lobo_selection("M_H_index")

cat("B. LORO-CV Generalizability Selection Pathway:\n")
cat(sprintf("  ★ Best Biomass Model:    \"%s\" (LORO-CV MAE = %.4f)\n  Formula: %s\n\n", 
            sel_f2_lobo_biomass$name, sel_f2_lobo_biomass$mae, sel_f2_lobo_biomass$formula))
cat(sprintf("  ★ Best Metabolism Model: \"%s\" (LORO-CV MAE = %.4f)\n  Formula: %s\n\n", 
            sel_f2_lobo_metabolism$name, sel_f2_lobo_metabolism$mae, sel_f2_lobo_metabolism$formula))

# -----------------------------------------------------------------------------
# 4. Save Quantitative Results to Permanent CSV Deliverable
# -----------------------------------------------------------------------------
results_df <- do.call(rbind, lapply(comparison_records, as.data.frame))
write_csv(results_df, "outputs/metabolic_scaling_model_selection.csv")
cat("✓ Successfully exported all comparative model selection data to:\n")
cat("  - outputs/metabolic_scaling_model_selection.csv\n\n")

# --- 5. Robust Integrity Verification Check ----------------------------------
if (sel_f1_biomass$name == sel_f1_metabolism$name &&
    sel_f2_aic_biomass$name == sel_f2_aic_metabolism$name &&
    sel_f2_lobo_biomass$name == sel_f2_lobo_metabolism$name) {
  cat("=====================================================================\n")
  cat("✓ INTEGRITY CHECK PASSED: Formulations are 100% stable and identical.\n")
  cat("=====================================================================\n")
} else {
  cat("=====================================================================\n")
  cat("⚠ NOTE: There is a divergence in selected covariate formulations.\n")
  cat("=====================================================================\n")
}
