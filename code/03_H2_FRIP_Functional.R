# =============================================================================
# 03_H2_FRIP_Functional.R
#
# H2 Hypothesis Testing: Functional defaunation signal via FRIP.
#
# Prediction: Amazon > Congo for FRIP (nutrient pump broken where megafauna
# are depleted → flooding more strongly predicts productivity).
#
# Analyses:
#   1. FRIP ANOVA: basin + country + protection
#   2. FRIP ~ Defaunation Index (OLS)
#   3. DI index comparison (AIC model selection)
#   4. Spatial cross-validation (3×3 tile leave-one-out)
#   5. Basin-split denoising validation
#   6. Multi-scale full loop
#
# Input:
#   - outputs/rds/loaded_data.rds
#
# Output:
#   - outputs/rds/h2_results.rds
#
# Dependencies:
#   Loaded via source("code/functions/frip_analysis.R")
# =============================================================================

# --- Setup -------------------------------------------------------------------

library(terra)
library(dplyr)
library(tibble)

source("code/functions/frip_analysis.R")

cat("=== 03: H2 — FRIP Functional Analysis ===\n\n")


# --- Load data ---------------------------------------------------------------

cat("Loading data...\n")
loaded <- readRDS("outputs/rds/loaded_data.rds")
native_stacks     <- lapply(loaded$native_stacks, unwrap)
multiscale_stacks <- lapply(loaded$multiscale_stacks, function(scale_list) {
  lapply(scale_list, unwrap)
})
pa_rast           <- unwrap(loaded$pa_rast)
basins_r          <- unwrap(loaded$basins_r)
countries_r       <- unwrap(loaded$countries_r)
basins_v          <- unwrap(loaded$basins_v)

RDS_DIR <- file.path("outputs", "rds")
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

# Defaunation index rasters (loaded using the helper function from load_data.R)
source("code/functions/load_data.R")

# Pick a representative multi-scale stack (e.g. 25 km)
stack_25k <- merge(multiscale_stacks[["25000"]][["Congo"]],
                   multiscale_stacks[["25000"]][["Amazon"]])

di_indices <- load_defaunation_indices(stack_25k)
di_bl <- di_indices[["DI"]]
di_bogoni <- di_indices[["DI_Bogoni"]]
di_rast <- di_bl  # primary index


# --- Analysis 1: FRIP ANOVA -------------------------------------------------

cat("\n--- 1. FRIP ANOVA: basin + country + protection ---\n")

frip_anova <- test_frip_by_basin_country_pa(stack_25k, basins_r,
                                             countries_r, pa_rast)
print(summary(frip_anova$aov))
cat("  Compact letters:\n")
print(frip_anova$letters)


# --- Analysis 2: FRIP ~ DI (OLS) --------------------------------------------

cat("\n--- 2. FRIP ~ Defaunation Index (OLS) ---\n")

frip_di_ols <- test_frip_vs_di(stack_25k, di_rast)
cat(sprintf("  R² = %.4f, p = %.2e\n",
            frip_di_ols$glance$r.squared,
            frip_di_ols$tidy$p.value[2]))


# --- Analysis 3: DI index comparison ----------------------------------------

cat("\n--- 3. AIC model selection: BL vs Bogoni DI ---\n")

di_comparison <- compare_di_indices(stack_25k, di_bl, di_bogoni)
print(di_comparison)


# --- Analysis 4: Spatial CV (3×3 tile LOO) -----------------------------------

cat("\n--- 4. Spatial cross-validation (3×3 tile LOO) ---\n")

cv_results <- run_denoising_cv(stack_25k, di_rast, basins_v)
cat(sprintf("  Mean R² ratio: %.3f\n", mean(cv_results$ratio, na.rm = TRUE)))
print(cv_results)


# --- Analysis 5: Basin-split denoising --------------------------------------

cat("\n--- 5. Basin-split denoising validation ---\n")

split_results <- run_basin_split_denoising(stack_25k, di_rast, basins_v)
print(split_results)


# --- Analysis 6: Multi-scale loop -------------------------------------------

cat("\n--- 6. Full multi-scale H2 loop ---\n")

h2_multiscale <- run_h2_multiscale(
  multiscale_stacks,
  basins_r   = basins_r,
  countries_r = countries_r,
  pa_r       = pa_rast,
  di_rast    = di_rast,
  basins_v   = basins_v
)
print(h2_multiscale)


# --- Save results ------------------------------------------------------------

cat("\nSaving H2 results...\n")

h2_results <- list(
  frip_anova    = frip_anova,
  frip_di_ols   = frip_di_ols,
  di_comparison = di_comparison,
  cv_results    = cv_results,
  split_results = split_results,
  h2_multiscale = h2_multiscale
)
saveRDS(h2_results, file.path(RDS_DIR, "h2_results.rds"))

cat("=== 03: H2 Done ===\n")
