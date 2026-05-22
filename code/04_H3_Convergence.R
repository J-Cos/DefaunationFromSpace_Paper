# =============================================================================
# 04_H3_Convergence.R
#
# H3 Hypothesis Testing: Structural–functional convergence.
#
# Prediction: UOI and FRIP spatially converge. Protected areas and pixels
# with high UOI (intact structure) should have low FRIP (intact function).
#
# Analyses:
#   1. PA-scale convergence (Spearman r: mean UOI vs mean FRIP per PA)
#   2. Pixel-scale convergence (per-basin Spearman UOI × FRIP)
#   3. Bivariate 2×2 classification map
#   4. Multi-scale convergence (UOI–FRIP r across 20 scales)
#   5. Dual-signal protection test (UOI and FRIP: protected vs unprotected)
#
# Input:
#   - outputs/rds/loaded_data.rds
#
# Output:
#   - outputs/rds/h3_results.rds
#
# Dependencies:
#   Loaded via source("code/functions/convergence_analysis.R")
# =============================================================================

# --- Setup -------------------------------------------------------------------

library(terra)
library(dplyr)
library(tibble)

source("code/functions/convergence_analysis.R")

cat("=== 04: H3 — Convergence Analysis ===\n\n")


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


# --- Analysis 1: PA-scale convergence ---------------------------------------

cat("\n--- 1. PA-scale convergence: Spearman r (mean UOI vs mean FRIP) ---\n")

pa_convergence <- test_pa_convergence(native_stacks, pa_rast)
cat(sprintf("  Spearman rho = %.3f, p = %.2e\n",
            pa_convergence$rho,
            pa_convergence$p_value))
cat(sprintf("  N protected areas = %d\n", nrow(pa_convergence$pa_df)))


# --- Analysis 2: Pixel-scale convergence ------------------------------------

cat("\n--- 2. Pixel-scale convergence: per-basin Spearman UOI × FRIP ---\n")

pixel_convergence <- test_pixel_convergence(native_stacks)
for (basin in names(pixel_convergence)) {
  res <- pixel_convergence[[basin]]
  cat(sprintf("  %s: rho = %.3f, p = %.2e, n = %d\n",
              basin, res$rho, res$p_value, res$n_pixels))
}


# --- Analysis 3: Bivariate 2×2 classification --------------------------------

cat("\n--- 3. Bivariate classification (2×2: UOI × FRIP) ---\n")

bivariate_rast <- classify_bivariate(native_stacks)
cat("  Class frequencies:\n")
print(freq(bivariate_rast))

# Save bivariate raster for figure generation
writeRaster(bivariate_rast,
            file.path("outputs", "bivariate_classification.tif"),
            overwrite = TRUE)


# --- Analysis 4: Multi-scale convergence ------------------------------------

cat("\n--- 4. Multi-scale convergence (20 scales) ---\n")

ms_convergence <- test_convergence_multiscale(multiscale_stacks, pa_rast)
cat(sprintf("  Significant pixel-level convergence at %d / %d scales\n",
            sum(ms_convergence$p_pixel < 0.05),
            nrow(ms_convergence)))
print(ms_convergence)


# --- Analysis 5: Dual-signal protection test --------------------------------

cat("\n--- 5. Protection dual-signal test ---\n")

dual_signal <- test_protection_dual_signal(native_stacks, pa_rast)
print(dual_signal)


# --- Save results ------------------------------------------------------------

cat("\nSaving H3 results...\n")

h3_results <- list(
  pa_convergence    = pa_convergence,
  pixel_convergence = pixel_convergence,
  bivariate_rast    = "outputs/bivariate_classification.tif",
  ms_convergence    = ms_convergence,
  dual_signal       = dual_signal
)
saveRDS(h3_results, file.path(RDS_DIR, "h3_results.rds"))

cat("=== 04: H3 Done ===\n")
