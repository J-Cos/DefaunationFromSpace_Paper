# =============================================================================
# code/10_Collect_Results_Stats.R
#
# Post-hoc harvester: reads all pipeline outputs and produces a long-form CSV
# (outputs/results_statistics.csv) containing every numeric result used in the
# manuscript.  No existing scripts are modified.
#
# Excluded: GEE-derived stats (GEDI coverage %, median shots per cell, basin-
# wide UOI raster ranges).
# Excluded: parameters we set (clustering threshold, trap-day filter, body-mass
# thresholds, temporal weight brackets).
#
# Usage:
#   Rscript code/10_Collect_Results_Stats.R
# =============================================================================

library(readr)
library(dplyr)
library(mgcv)
library(terra)
library(jsonlite)
library(statmod)

config_path <- if (file.exists("code/config.json")) "code/config.json" else "config.json"
config <- jsonlite::read_json(config_path)
day_range_coeff <- config$allometric_scaling$day_range_coeff
day_range_exp   <- config$allometric_scaling$day_range_exp

cat("=== 10: Collecting Results Statistics ===\n\n")

# ---------------------------------------------------------------------------
# Helper: accumulator
# ---------------------------------------------------------------------------
stats <- tibble(
  section    = character(),
  stat_id    = character(),
  description = character(),
  value      = numeric(),
  ci_lower   = numeric(),
  ci_upper   = numeric(),
  p_value    = numeric(),
  unit       = character(),
  source_file = character()
)

add <- function(sec, id, desc, val,
                ci_lo = NA_real_, ci_hi = NA_real_,
                pval  = NA_real_, unit = "", src = "") {
  stats <<- bind_rows(stats, tibble(
    section     = sec,
    stat_id     = id,
    description = desc,
    value       = as.numeric(val),
    ci_lower    = as.numeric(ci_lo),
    ci_upper    = as.numeric(ci_hi),
    p_value     = as.numeric(pval),
    unit        = unit,
    source_file = src
  ))
}

# Helper to compute profile likelihood confidence intervals (CIs) for Tweedie GLMs.
# Profile-likelihood CIs are statistically superior to standard Wald-based CIs 
# (estimate +/- 1.96 * SE) under log links and moderate sample sizes (N ~ 100), 
# preventing coverage distortion.
#
# Arguments:
#   model_obj: A fitted gam object containing a Tweedie family.
#
# Returns:
#   A matrix containing the 2.5% and 97.5% confidence intervals. If profiling 
#   fails to converge, falls back gracefully to Wald intervals.
get_profile_cis <- function(model_obj) {
  tryCatch({
    p_val <- model_obj$family$getTheta(TRUE)
    df <- model_obj$model
    df$weights_var <- model_obj$prior.weights
    glm_fit <- glm(
      formula(model_obj),
      data = df,
      family = tweedie(var.power = p_val, link.power = 0),
      weights = weights_var,
      start = coef(model_obj),
      control = glm.control(maxit = 500)
    )
    suppressMessages(confint(glm_fit))
  }, error = function(e) {
    ptab <- summary(model_obj)$p.table
    ci <- cbind(
      ptab[, "Estimate"] - 1.96 * ptab[, "Std. Error"],
      ptab[, "Estimate"] + 1.96 * ptab[, "Std. Error"]
    )
    colnames(ci) <- c("2.5 %", "97.5 %")
    ci
  })
}

# Computes the Pearson r and Spearman rho correlations between the 
# LOBO-selected (parsimonious) and AICc-selected model predictions.
# Replicates the 20 km prediction grids comparison (symmetrical layout)
# across all tropical forest basin regions.
#
# Arguments:
#   ptab_lobo: Coefficient table (summary p.table) of the LOBO model.
#   ptab_aic: Coefficient table (summary p.table) of the AIC model.
#
# Returns:
#   A list containing the overall Pearson/Spearman coefficients and 
#   basin-specific lists, or NULL if files are missing.
compute_prediction_correlations <- function(ptab_lobo, ptab_aic) {
  basins <- c("Congo", "Amazon", "SE_Asia")
  all_lobo <- c()
  all_aic  <- c()
  per_basin <- list()

  for (b in basins) {
    r_path <- file.path("outputs", "EOdata",
                        sprintf("analysis_stack_5000_%s.tif", b))
    if (!file.exists(r_path)) {
      r_path <- file.path("outputs", "synthetic_EOdata",
                          sprintf("analysis_stack_5000_%s.tif", b))
    }
    if (!file.exists(r_path)) next

    r <- rast(r_path)
    # Aggregate 5 km → 20 km (factor 4)
    r_20 <- aggregate(r, fact = 4, fun = "mean", na.rm = TRUE)

    # Assign band names for reliable access (band 3 = uoi in 12-band stacks)
    band_names <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                    "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
    if (nlyr(r_20) == 12) names(r_20) <- band_names
    if (nlyr(r_20) == 11) names(r_20) <- band_names[-4]  # no uoi_sd
    uoi_vals <- values(r_20[["uoi"]])
    valid    <- !is.na(uoi_vals)
    uoi_v    <- uoi_vals[valid]

    if (length(uoi_v) == 0) next

    # LOBO: log(BH) = intercept + slope * UOI
    pred_lobo <- ptab_lobo["(Intercept)", "Estimate"] +
                 ptab_lobo["uoi", "Estimate"] * uoi_v

    # AIC model: need elephant presence.  Use continent-level assignment.
    ele <- ifelse(b %in% c("Congo", "SE_Asia"), 1, 0)
    b_int_aic  <- ptab_aic["(Intercept)", "Estimate"]
    b_uoi_aic  <- ptab_aic["uoi", "Estimate"]
    b_ele_aic  <- ptab_aic["elephant_present_possiblePresent", "Estimate"]
    b_uxe_aic  <- ptab_aic["uoi:elephant_present_possiblePresent", "Estimate"]

    pred_aic <- b_int_aic + b_uoi_aic * uoi_v +
                b_ele_aic * ele + b_uxe_aic * uoi_v * ele

    all_lobo <- c(all_lobo, pred_lobo)
    all_aic  <- c(all_aic,  pred_aic)

    per_basin[[b]] <- list(
      r_pearson  = cor(pred_lobo, pred_aic, method = "pearson"),
      rho_spearman = cor(pred_lobo, pred_aic, method = "spearman")
    )
  }

  if (length(all_lobo) == 0) return(NULL)

  list(
    overall_r   = cor(all_lobo, all_aic, method = "pearson"),
    overall_rho = cor(all_lobo, all_aic, method = "spearman"),
    per_basin   = per_basin
  )
}


# ───────────────────────────────────────────────────────────────────────────
# SECTION 2 — Camera trap datasets
# ───────────────────────────────────────────────────────────────────────────
cat("Section 2: Camera trap sampling descriptors...\n")

det  <- read_csv("outputs/camera_traps_joint_detections.csv",  show_col_types = FALSE)
dep  <- read_csv("outputs/camera_traps_joint_metrics.csv",     show_col_types = FALSE)
clus <- read_csv("outputs/camera_traps_cluster_level_metrics.csv", show_col_types = FALSE)

SRC_CT  <- "process_camera_traps.py"
SRC_VIS <- "visualise_camera_traps.py"
SRC_CLU <- "camera_traps_cluster_level_metrics.csv"

# --- Projects & deployments ---
add("2", "n_wi_projects",      "Number of Wildlife Insights projects",
    n_distinct(dep$project_name), src = SRC_CT, unit = "count")

add("2", "n_deployments_total", "Total camera deployments",
    nrow(dep), src = SRC_CT, unit = "count")

for (r in c("Amazon", "Congo", "SE_Asia")) {
  add("2", paste0("n_deployments_", tolower(r)),
      paste0("Camera deployments – ", r),
      sum(dep$region == r), src = SRC_CT, unit = "count")
}

# --- Trap-days ---
add("2", "total_trap_days_all_deployments",
    "Total survey effort across all deployments (trap-days)",
    sum(dep$trap_days, na.rm = TRUE), src = SRC_CT, unit = "trap-days")

add("2", "total_trap_days_retained_clusters",
    "Total survey effort in retained clusters (trap-days)",
    sum(clus$trap_days, na.rm = TRUE), src = SRC_CLU, unit = "trap-days")

# --- Clusters (pre-filter) ---
n_clusters_prefilter <- n_distinct(det$cluster_id)

add("2", "n_clusters_prefilter",
    "Spatially independent camera-trap clusters (before GEDI overlap filter)",
    n_clusters_prefilter, src = SRC_VIS, unit = "count")

# --- Clusters (post-filter) ---
add("2", "n_clusters_retained",
    "Clusters retained after GEDI overlap and trap-day filters",
    nrow(clus), src = SRC_CLU, unit = "count")

for (r in c("Amazon", "Congo", "SE_Asia")) {
  add("2", paste0("n_clusters_", tolower(r)),
      paste0("Retained clusters – ", r),
      sum(clus$region == r), src = SRC_CLU, unit = "count")
}

# Clusters with mammal detections (n_detections_total > 0)
add("2", "n_clusters_with_detections",
    "Clusters yielding mammal detections",
    sum(clus$n_detections_total > 0), src = SRC_CLU, unit = "count")

# --- Detection events ---
add("2", "n_detection_events",
    "Total independent detection events across retained clusters",
    sum(clus$n_detections_total), src = SRC_CLU, unit = "count")

# --- Body mass range ---
bm <- det$body_mass_kg[!is.na(det$body_mass_kg) & det$body_mass_kg > 0]
add("2", "body_mass_min_kg",
    "Minimum detected body mass",
    min(bm), src = SRC_CT, unit = "kg")
add("2", "body_mass_max_kg",
    "Maximum detected body mass",
    max(bm), src = SRC_CT, unit = "kg")
add("2", "body_mass_orders_of_magnitude",
    "Orders of magnitude spanned by detected body masses",
    log10(max(bm) / min(bm)), src = SRC_CT, unit = "log10 ratio")

# --- Biomass index ranges per basin ---
for (r in c("Amazon", "Congo", "SE_Asia")) {
  sub <- clus[clus$region == r, ]
  add("2", paste0("bh_min_", tolower(r)),
      paste0("Minimum B_H – ", r),
      min(sub$B_H_index), src = SRC_CLU, unit = "index")
  add("2", paste0("bh_max_", tolower(r)),
      paste0("Maximum B_H – ", r),
      max(sub$B_H_index), src = SRC_CLU, unit = "index")
}

# Congo clusters with megaherbivore biomass > 1000
add("2", "congo_clusters_bh_gt1000",
    "Congo clusters with B_H >1000 (megaherbivore)",
    sum(clus$region == "Congo" & clus$B_H_gt1000 > 0),
    src = SRC_CLU, unit = "count")

# Megaherbivore fraction (>50 kg body mass, as stored in megafauna_fraction_gt50) by basin
for (r in c("Amazon", "Congo", "SE_Asia")) {
  sub <- clus[clus$region == r & clus$B_H_index > 0, ]
  if (nrow(sub) > 0) {
    mean_frac_50 <- mean(sub$megafauna_fraction_gt50, na.rm = TRUE)
  } else {
    mean_frac_50 <- 0
  }
  add("2", paste0("megaherbivore_gt50_pct_", tolower(r)),
      paste0("Mean megaherbivore biomass fraction (>50 kg) – ", r),
      mean_frac_50, src = SRC_CLU, unit = "percent")
}

# Megaherbivore fraction (>100 kg body mass, as stored in megafauna_fraction_gt100) by basin
# This is what the manuscript reports as "megaherbivore biomass fraction"
for (r in c("Amazon", "Congo", "SE_Asia")) {
  sub <- clus[clus$region == r & clus$B_H_index > 0, ]
  if (nrow(sub) > 0) {
    mean_frac_100 <- mean(sub$megafauna_fraction_gt100, na.rm = TRUE)
  } else {
    mean_frac_100 <- 0
  }
  add("2", paste0("megaherbivore_gt100_pct_", tolower(r)),
      paste0("Mean megaherbivore biomass fraction (>100 kg) – ", r),
      mean_frac_100, src = SRC_CLU, unit = "percent")
}

# Proboscidean biomass fraction by basin (cluster-level aggregation)
# Uses robust detections with cluster_id to match visualise_camera_traps.py logic
robust_det <- read_csv("outputs/camera_traps_robust_detections.csv", show_col_types = FALSE)

for (r in c("Amazon", "Congo", "SE_Asia")) {
  sub_det <- robust_det[robust_det$region == r &
                          !is.na(robust_det$body_mass_kg) &
                          robust_det$body_mass_kg > 0 &
                          robust_det$taxon_quality %in% c("species", "genus", "family"), ]

  if (nrow(sub_det) == 0) {
    add("2", paste0("proboscidean_biomass_pct_", tolower(r)),
        paste0("Proboscidean fraction of total standing biomass – ", r),
        0, src = SRC_VIS, unit = "percent")
    next
  }

  sub_det$day_range_km <- day_range_coeff * (sub_det$body_mass_kg ^ day_range_exp)
  # cluster-level trap_days from the cluster file
  sub_det <- sub_det %>%
    left_join(
      clus %>% select(cluster_id, cluster_trap_days = trap_days),
      by = "cluster_id"
    )
  sub_det$corrected_rai <- (sub_det$n_detections / sub_det$cluster_trap_days * 100) / sub_det$day_range_km
  sub_det$biomass_contrib <- sub_det$corrected_rai * sub_det$body_mass_kg

  # Aggregate to cluster level, then compute mean Proboscidean fraction
  cluster_totals <- sub_det %>%
    group_by(cluster_id) %>%
    summarise(
      total_bm = sum(biomass_contrib, na.rm = TRUE),
      probo_bm = sum(biomass_contrib[order == "Proboscidea"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(probo_pct = ifelse(total_bm > 0, probo_bm / total_bm * 100, 0))

  add("2", paste0("proboscidean_biomass_pct_", tolower(r)),
      paste0("Mean Proboscidean biomass fraction across clusters – ", r),
      mean(cluster_totals$probo_pct, na.rm = TRUE), src = SRC_VIS, unit = "percent")
}

cat("  ✓ Section 2 done.\n")

# ───────────────────────────────────────────────────────────────────────────
# SECTION 3 — Framework 1: mammal biomass predicts UOI (Beta Regression)
# ───────────────────────────────────────────────────────────────────────────
cat("Section 3: Framework 1 model results...\n")

fw1_sel <- read_csv("outputs/framework1_covariate_model_selection.csv", show_col_types = FALSE)
fw1_mod <- readRDS("outputs/framework1_best_model.RDS")
SRC_F1  <- "03_Framework1_Analysis.R"

add("3", "fw1_n_candidate_models",
    "Number of candidate models evaluated",
    nrow(fw1_sel), src = SRC_F1, unit = "count")

# Best model
add("3", "fw1_best_aicc",
    "Best model AICc",
    fw1_sel$AICc[1], src = SRC_F1, unit = "AICc")
add("3", "fw1_best_dev_expl",
    "Best model deviance explained",
    fw1_sel$dev_expl[1], src = SRC_F1, unit = "proportion")

# Coefficients with CIs and p-values from the saved GAM object
ptab1  <- summary(fw1_mod)$p.table
vcov1  <- vcov(fw1_mod)
coef_names1 <- rownames(ptab1)

for (cn in coef_names1) {
  est  <- ptab1[cn, "Estimate"]
  se   <- ptab1[cn, "Std. Error"]
  pv   <- ptab1[cn, "Pr(>|z|)"]
  ci_l <- est - 1.96 * se
  ci_u <- est + 1.96 * se

  clean <- gsub("elephant_present_possiblePresent", "ElephantPossible",
           gsub("\\(Intercept\\)", "intercept", cn))

  add("3", paste0("fw1_coef_", tolower(gsub("[^a-zA-Z0-9]", "_", clean))),
      paste0("FW1 best model coefficient – ", clean, " (logit scale)"),
      est, ci_lo = ci_l, ci_hi = ci_u, pval = pv,
      src = SRC_F1, unit = "logit")
}

# Biomass-only model (simplest formulation)
bm_only_row <- fw1_sel[grepl("Biomass Only", fw1_sel$Model), ]
if (nrow(bm_only_row) > 0) {
  add("3", "fw1_biomass_only_dev_expl",
      "Biomass-only model deviance explained",
      bm_only_row$dev_expl[1], src = SRC_F1, unit = "proportion")
}

# --- Weighting sensitivity (Framework 1) ---
wsens1 <- read_csv("outputs/weighting_sensitivity_framework1.csv", show_col_types = FALSE)
SRC_WS <- "weighting_sensitivity_framework1.csv"

for (i in seq_len(nrow(wsens1))) {
  row <- wsens1[i, ]
  tag <- tolower(gsub("[^a-zA-Z0-9]", "_", row$Weighting))
  add("3", paste0("fw1_wsens_slope_", tag),
      paste0("FW1 sensitivity biomass slope – ", row$Weighting),
      row$Slope, pval = row$P_Value,
      src = SRC_WS, unit = "logit")
  add("3", paste0("fw1_wsens_devexpl_", tag),
      paste0("FW1 sensitivity deviance explained – ", row$Weighting),
      row$DevExpl, src = SRC_WS, unit = "percent")
}

cat("  ✓ Section 3 done.\n")

# ───────────────────────────────────────────────────────────────────────────
# SECTION 4 — Framework 2: UOI predicts mammal biomass (Tweedie GLM)
# ───────────────────────────────────────────────────────────────────────────
cat("Section 4: Framework 2 model results...\n")

fw2_sel  <- read_csv("outputs/framework2_covariate_model_selection.csv",  show_col_types = FALSE)
fw2_mod  <- readRDS("outputs/framework2_best_model.RDS")       # LOBO-selected
fw2_aic  <- readRDS("outputs/framework2_best_model_aic.RDS")   # AIC-selected
SRC_F2   <- "04_Framework2_Analysis.R"

add("4", "fw2_n_candidate_models",
    "Number of candidate models evaluated",
    nrow(fw2_sel), src = SRC_F2, unit = "count")

# --- AICc-selected model ---
# Find the AICc-best row (lowest Full_AICc in the full selection table)
aicc_best_row <- fw2_sel[which.min(fw2_sel$Full_AICc), ]

add("4", "fw2_aic_best_dev_expl",
    "AICc-selected model deviance explained",
    aicc_best_row$Full_DevExpl, src = SRC_F2, unit = "proportion")
add("4", "fw2_aic_best_aicc",
    "AICc-selected model AICc value",
    aicc_best_row$Full_AICc, src = SRC_F2, unit = "AICc")



# AIC model coefficients
ptab_aic  <- summary(fw2_aic)$p.table
vcov_aic  <- vcov(fw2_aic)
prof_cis_aic <- get_profile_cis(fw2_aic)

for (cn in rownames(ptab_aic)) {
  est  <- ptab_aic[cn, "Estimate"]
  pv   <- ptab_aic[cn, "Pr(>|t|)"]
  ci_l <- prof_cis_aic[cn, 1]
  ci_u <- prof_cis_aic[cn, 2]

  clean <- gsub("elephant_present_possiblePresent", "ElephantPossible",
           gsub("uoi:elephant_present_possiblePresent", "UOI_x_ElephantPossible",
           gsub("\\(Intercept\\)", "intercept", cn)))

  add("4", paste0("fw2_aic_coef_", tolower(gsub("[^a-zA-Z0-9]", "_", clean))),
      paste0("FW2 AIC model coefficient – ", clean, " (log-link scale)"),
      est, ci_lo = ci_l, ci_hi = ci_u, pval = pv,
      src = SRC_F2, unit = "log-link")
}

# Slope within elephant range = beta_uoi + beta_uoi:elephant
# (requires variance-covariance for combined CI)
if ("uoi" %in% rownames(ptab_aic) &&
    "uoi:elephant_present_possiblePresent" %in% rownames(ptab_aic)) {

  b_uoi  <- ptab_aic["uoi", "Estimate"]
  b_int  <- ptab_aic["uoi:elephant_present_possiblePresent", "Estimate"]
  slope_within <- b_uoi + b_int

  var_sum <- vcov_aic["uoi", "uoi"] +
             vcov_aic["uoi:elephant_present_possiblePresent",
                      "uoi:elephant_present_possiblePresent"] +
             2 * vcov_aic["uoi", "uoi:elephant_present_possiblePresent"]
  se_sum <- sqrt(var_sum)
  z_sum  <- slope_within / se_sum
  p_sum  <- 2 * pnorm(-abs(z_sum))

  add("4", "fw2_aic_slope_within_elephant",
      "UOI–biomass slope within possible elephant range (log-link scale)",
      slope_within,
      ci_lo = slope_within - 1.96 * se_sum,
      ci_hi = slope_within + 1.96 * se_sum,
      pval  = p_sum, src = SRC_F2, unit = "log-link")
}

# --- LOBO-selected (parsimonious) model ---
lobo_best_row <- fw2_sel[grepl("Parsimonious Selected", fw2_sel$Parsimony_Status), ]
if (nrow(lobo_best_row) == 0) lobo_best_row <- fw2_sel[1, ]

add("4", "fw2_lobo_best_dev_expl",
    "LORO-CV parsimonious model deviance explained",
    lobo_best_row$Full_DevExpl, src = SRC_F2, unit = "proportion")
add("4", "fw2_lobo_oos_rmse_log",
    "LORO-CV out-of-sample RMSE (log scale)",
    lobo_best_row$OOS_RMSE_log, src = SRC_F2, unit = "log-scale RMSE")
add("4", "fw2_lobo_oos_mae_log",
    "LORO-CV out-of-sample MAE (log scale)",
    lobo_best_row$OOS_MAE_log, src = SRC_F2, unit = "log-scale MAE")
add("4", "fw2_lobo_oos_r2_log",
    "LORO-CV out-of-sample R2 (log scale)",
    lobo_best_row$OOS_R2_log, src = SRC_F2, unit = "R2")

# LOBO model coefficients
ptab_lobo <- summary(fw2_mod)$p.table
prof_cis_lobo <- get_profile_cis(fw2_mod)

for (cn in rownames(ptab_lobo)) {
  est  <- ptab_lobo[cn, "Estimate"]
  pv   <- ptab_lobo[cn, "Pr(>|t|)"]
  ci_l <- prof_cis_lobo[cn, 1]
  ci_u <- prof_cis_lobo[cn, 2]

  clean <- gsub("\\(Intercept\\)", "intercept", cn)

  add("4", paste0("fw2_lobo_coef_", tolower(gsub("[^a-zA-Z0-9]", "_", clean))),
      paste0("FW2 LOBO model coefficient – ", clean, " (log-link scale)"),
      est, ci_lo = ci_l, ci_hi = ci_u, pval = pv,
      src = SRC_F2, unit = "log-link")
}

# Fold-change per 0.01 UOI increase (Tweedie log-link: exp(beta * delta))
if ("uoi" %in% rownames(ptab_lobo)) {
  b_uoi_lobo <- ptab_lobo["uoi", "Estimate"]
  delta      <- 0.01

  fold_change    <- exp(b_uoi_lobo * delta)
  fold_change_lo <- exp(prof_cis_lobo["uoi", 1] * delta)
  fold_change_hi <- exp(prof_cis_lobo["uoi", 2] * delta)

  add("4", "fw2_fold_change_per_001_uoi",
      "Fold-change in expected biomass per 0.01 UOI increase",
      fold_change,
      ci_lo = fold_change_lo, ci_hi = fold_change_hi,
      src = SRC_F2, unit = "fold")
}

# Deviance explained difference (AIC vs LOBO)
add("4", "fw2_dev_expl_difference",
    "Additional deviance explained by AIC model vs LOBO model",
    as.numeric(aicc_best_row$Full_DevExpl) - as.numeric(lobo_best_row$Full_DevExpl),
    src = SRC_F2, unit = "proportion")

# --- Weighting sensitivity (Framework 2) ---
if (file.exists("outputs/weighting_sensitivity_framework2.csv")) {
  wsens2 <- read_csv("outputs/weighting_sensitivity_framework2.csv", show_col_types = FALSE)
  SRC_WS2 <- "weighting_sensitivity_framework2.csv"

  for (i in seq_len(nrow(wsens2))) {
    row <- wsens2[i, ]
    tag <- tolower(gsub("[^a-zA-Z0-9]", "_", row$Weighting))
    add("4", paste0("fw2_wsens_slope_amazon_", tag),
        paste0("FW2 sensitivity UOI slope (Amazon) – ", row$Weighting),
        row$Slope_Amazon, pval = row$P_Amazon,
        src = SRC_WS2, unit = "log-link")
    if ("Slope_Congo_Diff" %in% names(row)) {
      add("4", paste0("fw2_wsens_slope_congo_diff_", tag),
          paste0("FW2 sensitivity Congo slope difference – ", row$Weighting),
          row$Slope_Congo_Diff, pval = row$P_Congo_Diff,
          src = SRC_WS2, unit = "log-link")
    }
  }
}

cat("  ✓ Section 4 done.\n")

# ───────────────────────────────────────────────────────────────────────────
# SECTION 5 — Predictive biomass maps
# ───────────────────────────────────────────────────────────────────────────
cat("Section 5: Predictive map statistics...\n")

SRC_MAP <- "05_Predictive_Biomass_Maps.R"
SRC_COR <- "05b_Predictive_Columns_Correlation.R"

# Predictive equation from LOBO model: log(BH) = slope * UOI + intercept
if ("uoi" %in% rownames(ptab_lobo)) {
  add("5", "fw2_equation_slope",
      "Predictive equation UOI slope (log-link)",
      ptab_lobo["uoi", "Estimate"],
      src = SRC_MAP, unit = "log-link")
}
add("5", "fw2_equation_intercept",
    "Predictive equation intercept (log-link)",
    ptab_lobo["(Intercept)", "Estimate"],
    src = SRC_MAP, unit = "log-link")

# --- Pearson r and Spearman rho between LOBO and AIC model predictions ---
# Replicate the 05b logic: load 5 km rasters, aggregate to 20 km,
# predict with both models, correlate.
cors <- tryCatch(compute_prediction_correlations(ptab_lobo, ptab_aic), error = function(e) {
  cat("  ⚠ Could not compute prediction correlations:", e$message, "\n")
  NULL
})

if (!is.null(cors)) {
  add("5", "pred_pearson_r_overall",
      "Pearson r between LOBO and AIC predictions (all basins, 20 km)",
      cors$overall_r, src = SRC_COR, unit = "correlation")
  add("5", "pred_spearman_rho_overall",
      "Spearman rho between LOBO and AIC predictions (all basins, 20 km)",
      cors$overall_rho, src = SRC_COR, unit = "correlation")

  for (b in names(cors$per_basin)) {
    add("5", paste0("pred_pearson_r_", tolower(b)),
        paste0("Pearson r between predictions – ", b),
        cors$per_basin[[b]]$r_pearson, src = SRC_COR, unit = "correlation")
    add("5", paste0("pred_spearman_rho_", tolower(b)),
        paste0("Spearman rho between predictions – ", b),
        cors$per_basin[[b]]$rho_spearman, src = SRC_COR, unit = "correlation")
  }

  # Within-basin rank-order preservation (Spearman rho of predicted BH ranks)
  # This should be 1.0 since both models are monotonic in UOI within a basin
  # (the interaction only shifts slope/intercept between elephant present/absent)
  for (b in names(cors$per_basin)) {
    add("5", paste0("pred_rank_spearman_", tolower(b)),
        paste0("Within-basin rank-order Spearman rho – ", b),
        cors$per_basin[[b]]$rho_spearman, src = SRC_COR, unit = "correlation")
  }
}

cat("  ✓ Section 5 done.\n")

# ───────────────────────────────────────────────────────────────────────────
# Write output
# ───────────────────────────────────────────────────────────────────────────
out_path <- "outputs/results_statistics.csv"
write_csv(stats, out_path)

cat(sprintf("\n=== Done. Wrote %d statistics to %s ===\n", nrow(stats), out_path))
