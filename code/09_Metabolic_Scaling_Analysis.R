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
  best_aic <- Inf
  best_name <- ""
  best_formula <- NULL
  
  for (n in names(f1_templates)) {
    f_str <- gsub("\\{Index\\}", index_name, f1_templates[[n]])
    f <- as.formula(f_str)
    
    fit <- tryCatch({
      gam(f, data = joined_data, family = betar(link = "logit"), weights = w_combined_norm)
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      a <- AIC(fit)
      dev_expl <- summary(fit)$dev.expl
      
      # Record all model fits for complete tracking
      comparison_records[[length(comparison_records) + 1]] <<- list(
        Framework = "Framework 1",
        Metric = "AIC & DevExpl",
        ModelLabel = n,
        IndexUsed = index_name,
        Formula = f_str,
        FullAIC = a,
        DevianceExplained = dev_expl,
        OOS_MAE = NA
      )
      
      if (a < best_aic) {
        best_aic <- a
        best_name <- n
        best_formula <- f_str
      }
    }
  }
  return(list(name = best_name, formula = best_formula, aic = best_aic))
}

sel_f1_biomass <- run_f1_selection("B_H_index")
sel_f1_metabolism <- run_f1_selection("M_H_index")

cat(sprintf("  ★ Best Biomass Model:    \"%s\" (AIC = %.2f)\n  Formula: %s\n\n", 
            sel_f1_biomass$name, sel_f1_biomass$aic, sel_f1_biomass$formula))
cat(sprintf("  ★ Best Metabolism Model: \"%s\" (AIC = %.2f)\n  Formula: %s\n\n", 
            sel_f1_metabolism$name, sel_f1_metabolism$aic, sel_f1_metabolism$formula))

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

# A. Standard AIC Selection Pathway (Full Sample Fit)
run_f2_aic_selection <- function(index_name) {
  best_aic <- Inf
  best_name <- ""
  best_formula <- NULL
  
  for (n in names(f2_templates)) {
    f_str <- gsub("\\{Index\\}", index_name, f2_templates[[n]])
    f <- as.formula(f_str)
    
    fit <- tryCatch({
      gam(f, data = joined_data, family = tw(), weights = w_combined_norm)
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      a <- AIC(fit)
      dev_expl <- summary(fit)$dev.expl
      
      comparison_records[[length(comparison_records) + 1]] <<- list(
        Framework = "Framework 2 AIC",
        Metric = "AIC & DevExpl",
        ModelLabel = n,
        IndexUsed = index_name,
        Formula = f_str,
        FullAIC = a,
        DevianceExplained = dev_expl,
        OOS_MAE = NA
      )
      
      if (a < best_aic) {
        best_aic <- a
        best_name <- n
        best_formula <- f_str
      }
    }
  }
  return(list(name = best_name, formula = best_formula, aic = best_aic))
}

sel_f2_aic_biomass <- run_f2_aic_selection("B_H_index")
sel_f2_aic_metabolism <- run_f2_aic_selection("M_H_index")

cat("A. Standard AIC Selection Pathway (Full Sample Fit):\n")
cat(sprintf("  ★ Best Biomass Model:    \"%s\" (AIC = %.2f)\n  Formula: %s\n\n", 
            sel_f2_aic_biomass$name, sel_f2_aic_biomass$aic, sel_f2_aic_biomass$formula))
cat(sprintf("  ★ Best Metabolism Model: \"%s\" (AIC = %.2f)\n  Formula: %s\n\n", 
            sel_f2_aic_metabolism$name, sel_f2_aic_metabolism$aic, sel_f2_aic_metabolism$formula))

# B. LORO-CV Generalizability Selection Pathway (Out-of-Sample CV MAE)
run_f2_lobo_selection <- function(index_name) {
  basins <- unique(joined_data$basin)
  
  results <- list()
  
  for (n in names(f2_templates)) {
    f_str <- gsub("\\{Index\\}", index_name, f2_templates[[n]])
    if (grepl("basin", f_str)) next # Skip basin terms due to LOBO folding
    
    f <- as.formula(f_str)
    errors <- c()
    
    for (b in basins) {
      train <- joined_data %>% filter(basin != b)
      test <- joined_data %>% filter(basin == b)
      
      if (nrow(train) == 0 || nrow(test) == 0) next
      
      fit <- tryCatch({
        gam(f, data = train, family = tw(), weights = w_combined_norm)
      }, error = function(e) NULL)
      
      if (!is.null(fit)) {
        pred <- predict(fit, newdata = test, type = "response")
        errors <- c(errors, mean(abs(log1p(test[[index_name]]) - log1p(pred))))
      }
    }
    
    if (length(errors) > 0) {
      avg_mae <- mean(errors)
      results[[n]] <- avg_mae
      
      # Also fit full model to get DevExpl
      fit_full <- tryCatch({
        gam(f, data = joined_data, family = tw(), weights = w_combined_norm)
      }, error = function(e) NULL)
      dev_expl <- if (!is.null(fit_full)) summary(fit_full)$dev.expl else NA
      full_aic <- if (!is.null(fit_full)) AIC(fit_full) else NA
      
      comparison_records[[length(comparison_records) + 1]] <<- list(
        Framework = "Framework 2 LORO-CV",
        Metric = "LORO-CV MAE",
        ModelLabel = n,
        IndexUsed = index_name,
        Formula = f_str,
        FullAIC = full_aic,
        DevianceExplained = dev_expl,
        OOS_MAE = avg_mae
      )
    }
  }
  
  best_name <- names(results)[which.min(unlist(results))]
  best_mae <- results[[best_name]]
  best_formula <- gsub("\\{Index\\}", index_name, f2_templates[[best_name]])
  
  return(list(name = best_name, formula = best_formula, mae = best_mae))
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
