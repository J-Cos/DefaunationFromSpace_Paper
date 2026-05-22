# =============================================================================
# 05_H4_Temporal_Trends.R
#
# H4 Hypothesis Testing: Temporal trends in the functional signal.
#
# Prediction: FRIP has strengthened over time where defaunation is increasing.
# Uses the GEE pre-computed `frip_mk_tau` band (Mann-Kendall trend τ of
# annual FRIP across 2001–2023).
#
# Analyses:
#   1. Summarise MK-tau by basin
#   2. MK-tau intercept-free model by PA
#   3. Protection effect on MK-tau (protected vs unprotected)
#   4. Multi-scale MK-tau loop (20 scales)
#
# Input:
#   - outputs/rds/loaded_data.rds
#
# Output:
#   - outputs/rds/h4_results.rds
#
# Dependencies:
#   Loaded via source("code/functions/temporal_analysis.R")
# =============================================================================

# --- Setup -------------------------------------------------------------------

library(terra)
library(dplyr)
library(tibble)

source("code/functions/temporal_analysis.R")

cat("=== 05: H4 — Temporal Trends Analysis ===\n\n")


# --- Load data ---------------------------------------------------------------

cat("Loading data...\n")
loaded <- readRDS("outputs/rds/loaded_data.rds")
multiscale_stacks <- lapply(loaded$multiscale_stacks, function(scale_list) {
  lapply(scale_list, unwrap)
})
pa_rast <- unwrap(loaded$pa_rast)

RDS_DIR <- file.path("outputs", "rds")
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

# Use a representative scale for native-scale analyses (25 km)
stacks_25k <- list(
  Congo  = multiscale_stacks[["25000"]][["Congo"]],
  Amazon = multiscale_stacks[["25000"]][["Amazon"]]
)


# --- Analysis 1: MK-tau summary by basin ------------------------------------

cat("\n--- 1. MK-tau summary by basin ---\n")

tau_summary <- summarise_mk_tau_by_basin(stacks_25k)
print(tau_summary)


# --- Analysis 2: MK-tau by PA (intercept-free model) ------------------------

cat("\n--- 2. MK-tau intercept-free model by PA ---\n")

tau_by_pa <- test_mk_tau_by_pa(stacks_25k, pa_rast)
cat("  Trend classification:\n")
tau_by_pa$coefficients %>%
  count(trend_class) %>%
  print()


# --- Analysis 3: Protection effect on MK-tau --------------------------------

cat("\n--- 3. Protection effect: unprotected tau > protected tau? ---\n")

tau_protection <- test_mk_tau_protection(stacks_25k, pa_rast)
cat(sprintf("  t = %.3f, p = %.2e\n",
            tau_protection$t_stat,
            tau_protection$p_value))
cat(sprintf("  Protected mean tau  = %.4f\n", tau_protection$mean_protected))
cat(sprintf("  Unprotected mean tau = %.4f\n", tau_protection$mean_unprotected))


# --- Analysis 4: Multi-scale MK-tau loop ------------------------------------

cat("\n--- 4. Multi-scale MK-tau analysis (20 scales) ---\n")

h4_multiscale <- run_h4_multiscale(multiscale_stacks, pa_rast)
cat(sprintf("  Protection significant at %d / %d scales\n",
            sum(h4_multiscale$protection_p < 0.05, na.rm = TRUE),
            nrow(h4_multiscale)))
print(h4_multiscale)


# --- Save results ------------------------------------------------------------

cat("\nSaving H4 results...\n")

h4_results <- list(
  tau_summary    = tau_summary,
  tau_by_pa      = tau_by_pa,
  tau_protection = tau_protection,
  h4_multiscale  = h4_multiscale
)
saveRDS(h4_results, file.path(RDS_DIR, "h4_results.rds"))

cat("=== 05: H4 Done ===\n")
