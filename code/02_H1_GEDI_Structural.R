# =============================================================================
# 02_H1_GEDI_Structural.R
#
# H1 Hypothesis Testing: Structural defaunation signal via GEDI UOI.
#
# Prediction: Congo > Amazon for Understory Openness Index (UOI).
# Forest understories are more open where megafauna are intact.
#
# Analyses:
#   1. Regional t-test (Congo vs Amazon UOI) at native scale
#   2. Two-way ANOVA: UOI ~ Region × Protection
#   3. Per-PA ANOVA + TukeyHSD + compact letters
#   4. Multi-scale regional t-test (20 scales)
#
# Input:
#   - outputs/rds/loaded_data.rds
#
# Output:
#   - outputs/rds/h1_results.rds
#
# Dependencies:
#   Loaded via source("code/functions/gedi_analysis.R")
# =============================================================================

# --- Setup -------------------------------------------------------------------

library(terra)
library(dplyr)
library(tibble)

source("code/functions/theme_pnas.R")
source("code/functions/gedi_analysis.R")

cat("=== 02: H1 — GEDI Structural Analysis ===\n\n")


# --- Load data ---------------------------------------------------------------

cat("Loading data...\n")
loaded <- readRDS("outputs/rds/loaded_data.rds")
native_stacks     <- lapply(loaded$native_stacks, unwrap)
multiscale_stacks <- lapply(loaded$multiscale_stacks, function(scale_list) {
  lapply(scale_list, unwrap)
})
pa_rast           <- unwrap(loaded$pa_rast)

RDS_DIR <- file.path("outputs", "rds")
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)


# --- Analysis 1: Regional t-test (native scale) -----------------------------

cat("\n--- 1. Regional t-test: Congo UOI > Amazon UOI ---\n")

regional_ttest <- test_regional_uoi(native_stacks)
cat(sprintf("  t = %.3f, p = %.2e, d = %.3f\n",
            regional_ttest$t_stat,
            regional_ttest$p_value,
            regional_ttest$cohens_d))
cat(sprintf("  Congo mean = %.4f, Amazon mean = %.4f\n",
            regional_ttest$congo_mean,
            regional_ttest$amazon_mean))


# --- Analysis 2: Two-way ANOVA (Region × Protection) ------------------------

cat("\n--- 2. Two-way ANOVA: UOI ~ Region × Protection ---\n")

region_prot_aov <- test_uoi_region_protection(native_stacks, pa_rast)
print(summary(region_prot_aov))


# --- Analysis 3: Per-PA ANOVA + TukeyHSD ------------------------------------

cat("\n--- 3. Per-PA ANOVA + TukeyHSD ---\n")

pa_results <- test_uoi_by_pa(native_stacks, pa_rast)
cat("  Compact letter display:\n")
print(pa_results$letters)


# --- Analysis 4: Multi-scale t-test loop -------------------------------------

cat("\n--- 4. Multi-scale regional t-test (20 scales) ---\n")

h1_multiscale <- run_h1_multiscale(multiscale_stacks, pa_rast)
cat(sprintf("  Significant at %d / %d scales\n",
            sum(h1_multiscale$signif),
            nrow(h1_multiscale)))
print(h1_multiscale)


# --- Save results ------------------------------------------------------------

cat("\nSaving H1 results...\n")

h1_results <- list(
  regional_ttest  = regional_ttest,
  region_prot_aov = region_prot_aov,
  pa_results      = pa_results,
  h1_multiscale   = h1_multiscale
)
saveRDS(h1_results, file.path(RDS_DIR, "h1_results.rds"))

cat("=== 02: H1 Done ===\n")
