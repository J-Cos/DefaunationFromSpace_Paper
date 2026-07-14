# =============================================================================
# code/04b_Framework_Figures.R
#
# Generates all framework analysis figures by loading saved model outputs from
# 03_Framework1_Analysis.R and 04_Framework2_Analysis.R.
#
# Outputs:
#   - figures/figure3.png (and .pdf) — Framework 1: scatter + model selection + residuals
#   - figures/figure4.png (and .pdf) — Framework 2: synthesis scatter + AICc bars + LOBO bars
#
# This script does NOT re-run any model fitting. It loads pre-fitted models
# and selection tables from the outputs/ directory.
# =============================================================================

library(ggplot2)
library(dplyr)
library(mgcv)
library(scales)
library(patchwork)

source("code/functions/calibration_helpers.R")
source("code/functions/theme_pnas.R")
source("code/functions/framework_helpers.R")

# =============================================================================
# SHARED DATA LOADING
# =============================================================================

cat("=== Loading Data & Models for Framework Figures ===\n\n")

# Load calibrated cluster data (shared by both frameworks)
joined_data <- extract_scale_data(5000)

# --- Framework 1 models ---
f1_best_model  <- readRDS("outputs/framework1_best_model.RDS")
f1_full_models <- readRDS("outputs/framework1_full_models.RDS")
f1_selection   <- read.csv("outputs/framework1_covariate_model_selection.csv",
                           stringsAsFactors = FALSE)

# --- Framework 2 models ---
f2_model_lobo  <- readRDS("outputs/framework2_best_model.RDS")
f2_model_aic   <- readRDS("outputs/framework2_best_model_aic.RDS")
f2_full_models <- readRDS("outputs/framework2_full_models.RDS")
f2_selection   <- read.csv("outputs/framework2_covariate_model_selection.csv",
                           stringsAsFactors = FALSE)

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

cat("✓ All data and models loaded successfully.\n\n")


# =============================================================================
# FIGURE 3 — Framework 1: GEDI Openness vs. Standing Biomass
# =============================================================================

generate_figure3 <- function() {
  cat("--- Generating Figure 3 (Framework 1) ---\n")

  best_model <- f1_best_model
  models_list <- f1_full_models
  results_df <- f1_selection

  # --- Panel A: Scatter Plot with Best Model Fit ---
  biomass_seq <- seq(0, 5000, length.out = 300)
  
  # Predict for 10th, 50th, and 90th percentiles of slope
  slope_vals <- quantile(joined_data$slope, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
  slope_range <- range(joined_data$slope, na.rm = TRUE)
  
  nd <- do.call(rbind, lapply(names(slope_vals), function(pct_name) {
    slp <- slope_vals[[pct_name]]
    df <- data.frame(B_H_index = biomass_seq, slope = slp, percentile = pct_name)
    best_formula_vars <- all.vars(formula(best_model))
    # Fill every predictor except B_H_index and slope with its median
    for (v in setdiff(best_formula_vars[-1], c("B_H_index", "slope"))) {
      df[[v]] <- median(joined_data[[v]], na.rm = TRUE)
    }
    preds <- predict(best_model, newdata = df, type = "link", se.fit = TRUE)
    df$fit   <- plogis(preds$fit)
    df$lower <- plogis(preds$fit - 1.96 * preds$se.fit)
    df$upper <- plogis(preds$fit + 1.96 * preds$se.fit)
    df
  }))

  p_a <- ggplot() +
    geom_ribbon(data = nd, aes(x = B_H_index, ymin = lower, ymax = upper, group = percentile, fill = slope),
                alpha = 0.08) +
    geom_point(data = joined_data,
               aes(x = B_H_index, y = uoi, fill = slope, shape = basin, size = trap_days, alpha = w_temp_cluster),
               color = "black", stroke = 0.3) +
    geom_line(data = nd, aes(x = B_H_index, y = fit, group = percentile, color = slope),
              linetype = "dashed", linewidth = 0.75) +
    scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22), labels = c("Amazon" = "Amazon", "Congo" = "Congo", "SE_Asia" = "SE Asia"), name = "Region:") +
    scale_fill_viridis_c(option = "viridis", name = "Slope:", limits = slope_range) +
    scale_color_viridis_c(option = "viridis", name = "Slope:", limits = slope_range) +
    scale_size_continuous(name = "Effort (Trap-days):", range = c(1.2, 4.0), breaks = c(100, 1000, 5000, 15000)) +
    scale_alpha_continuous(name = "Temporal Alignment:", range = c(0.25, 1.0), breaks = c(0.1, 0.5, 1.0), labels = c("Hist.", "Interm.", "Contemp.")) +
    scale_x_continuous(trans = "log1p", labels = comma_format(),
                       breaks = c(0, 10, 100, 1000, 3000), limits = c(0, 5000)) +
    scale_y_continuous(breaks = seq(0.92, 0.97, by = 0.01), limits = c(0.921, 0.969), expand = c(0, 0)) +
    labs(
      title = "A",
      x = "Mammal Standing Biomass Index (log1p scale)",
      y = "GEDI Understory Openness Index (UOI)"
    ) +
    guides(
      size = "none",
      alpha = "none",
      shape = guide_legend(override.aes = list(size = 2.5, fill = "gray50")),
      fill  = guide_colorbar(barheight = unit(1.8, "cm"), barwidth = unit(0.3, "cm")),
      color = guide_colorbar(barheight = unit(1.8, "cm"), barwidth = unit(0.3, "cm"))
    ) +
    theme_pnas(base_size = 7.5) +
    theme(
      legend.position = "right",
      aspect.ratio = 1.0,
      plot.title = element_text(face = "bold", size = 8.5, margin = margin(b = 6, t = 4)),
      axis.title.x = element_text(margin = margin(t = 4)),
      axis.title.y = element_text(margin = margin(r = 4)),
      plot.margin = margin(t = 1.5, r = 4, b = 1.5, l = 4, unit = "pt"),
      legend.title = element_text(size = 7.0, face = "bold"),
      legend.text = element_text(size = 6.5)
    )

  # --- Panel B: Model Selection Bar Plot ---
  results_df$plotmath_label <- sapply(1:nrow(results_df), function(i) {
    format_model_label(results_df$Model[i], models_list[[results_df$Model[i]]], VAR_MAP_FRAMEWORK1)
  })

  plot_sel_df <- results_df %>%
    head(15) %>%
    mutate(
      CleanName = gsub("^M[0-9\\.]+[a-z_]*: ", "", Model),
      CleanName = factor(CleanName, levels = rev(CleanName)),
      dev_expl_pct = dev_expl * 100
    )

  ordered_exprs <- plot_sel_df$plotmath_label[match(levels(plot_sel_df$CleanName), plot_sel_df$CleanName)]
  parsed_labels <- parse(text = ordered_exprs)

  p_b <- plot_model_selection_bars(
    plot_df = plot_sel_df,
    x_var = "dev_expl_pct",
    fill_var = "delta_AICc",
    fill_label = "Delta AICc",
    x_label = "Model Deviance Explained (%)",
    plot_title = "B",
    parsed_labels = parsed_labels
  ) + labs(y = NULL) + theme(aspect.ratio = 0.25, plot.margin = margin(t = 1.5, r = 4, b = 1.5, l = 4, unit = "pt"))

  # --- Assemble Figure 3 via Patchwork ---
  fig3_final <- p_a + p_b + plot_layout(ncol = 1)

  save_pnas(plot = fig3_final, filename = "figures/figure3.png", type = "double", height_cm = 16.0)
  save_pnas(plot = fig3_final, filename = "figures/figure3.pdf", type = "double", height_cm = 16.0)

  cat("✓ Saved Figure 3 to figures/figure3.png and .pdf\n\n")
}


# =============================================================================
# FIGURE 4 — Framework 2: Synthesis (Scatter + AICc Bars + LOBO Bars)
# =============================================================================

generate_figure4 <- function() {
  cat("--- Generating Figure 4 (Framework 2 Synthesis) ---\n")

  # --- Posterior simulation helper ---
  get_sim_predictions <- function(model_obj, newdata, n_sim = 1000) {
    beta <- coef(model_obj)
    Vp <- vcov(model_obj)
    X <- predict(model_obj, newdata = newdata, type = "lpmatrix")
    set.seed(42)
    beta_sim <- mgcv::rmvn(n_sim, beta, Vp)
    link_preds <- X %*% t(beta_sim)
    resp_preds <- exp(link_preds)
    newdata$fit <- apply(resp_preds, 1, median)
    newdata$lower <- apply(resp_preds, 1, quantile, probs = 0.025)
    newdata$upper <- apply(resp_preds, 1, quantile, probs = 0.975)
    return(newdata)
  }

  # --- Panel A: Scatter with Dual-Model Predictions & CIs ---
  uoi_seq <- seq(from = 0.918, to = 0.970, length.out = 300)
  covs_to_fill <- c("elevation", "slope", "hnd", "precip", "clay", "forest_fraction")

  # Predict LOBO (UOI-only) with posterior simulation CIs
  pred_df_lobo <- data.frame(uoi = uoi_seq)
  for (cv in covs_to_fill) pred_df_lobo[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
  pred_df_lobo <- get_sim_predictions(f2_model_lobo, pred_df_lobo)

  # Predict AIC (UOI * Elephant Strict) - Absent
  pred_df_aic_absent <- data.frame(uoi = uoi_seq)
  pred_df_aic_absent$elephant_present_strict <- factor("Absent", levels = c("Absent", "Present"))
  for (cv in covs_to_fill) pred_df_aic_absent[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
  pred_df_aic_absent <- get_sim_predictions(f2_model_aic, pred_df_aic_absent)

  # Predict AIC (UOI * Elephant Strict) - Present
  pred_df_aic_present <- data.frame(uoi = uoi_seq)
  pred_df_aic_present$elephant_present_strict <- factor("Present", levels = c("Absent", "Present"))
  for (cv in covs_to_fill) pred_df_aic_present[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
  pred_df_aic_present <- get_sim_predictions(f2_model_aic, pred_df_aic_present)

  custom_colors <- c(
    "UOI Only (LOBO)" = "black",
    "Elephants Absent (AICc)" = "#D55E00",
    "Elephants Present (AICc)" = "#0072B2"
  )

  custom_linetypes <- c(
    "UOI Only (LOBO)" = "dashed",
    "Elephants Absent (AICc)" = "solid",
    "Elephants Present (AICc)" = "solid"
  )

  p_a <- ggplot() +
    geom_point(data = joined_data,
               aes(x = uoi, y = B_H_index, fill = elephant_present_strict, size = trap_days, alpha = w_temp_cluster, shape = basin),
               color = "black", stroke = 0.3) +
    geom_ribbon(data = pred_df_lobo, aes(x = uoi, ymin = lower, ymax = upper), fill = "black", alpha = 0.08) +
    geom_ribbon(data = pred_df_aic_absent, aes(x = uoi, ymin = lower, ymax = upper), fill = "#D55E00", alpha = 0.08) +
    geom_ribbon(data = pred_df_aic_present, aes(x = uoi, ymin = lower, ymax = upper), fill = "#0072B2", alpha = 0.08) +

    geom_line(data = pred_df_lobo, aes(x = uoi, y = fit, color = "UOI Only (LOBO)", linetype = "UOI Only (LOBO)"), linewidth = 0.75) +
    geom_line(data = pred_df_aic_absent, aes(x = uoi, y = fit, color = "Elephants Absent (AICc)", linetype = "Elephants Absent (AICc)"), linewidth = 0.75) +
    geom_line(data = pred_df_aic_present, aes(x = uoi, y = fit, color = "Elephants Present (AICc)", linetype = "Elephants Present (AICc)"), linewidth = 0.75) +

    scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22), labels = c("Amazon" = "Amazon", "Congo" = "Congo", "SE_Asia" = "SE Asia"), name = "Region:") +
    scale_size_continuous(name = "Effort (Trap-days):", range = c(1.2, 4.0), breaks = c(100, 1000, 5000, 15000)) +
    scale_alpha_continuous(name = "Temporal Alignment:", range = c(0.25, 1.0), breaks = c(0.1, 0.5, 1.0), labels = c("Hist.", "Interm.", "Contemp.")) +
    scale_fill_manual(values = pal_elephant_binary, name = "Elephant Status:") +
    scale_color_manual(values = custom_colors, name = "Model Predictions:") +
    scale_linetype_manual(values = custom_linetypes, name = "Model Predictions:") +
    scale_x_continuous(breaks = seq(0.92, 0.97, by = 0.01)) +
    scale_y_continuous(trans = "log1p", labels = comma_format(), breaks = c(0, 10, 100, 1000, 3000)) +
    coord_cartesian(xlim = c(0.918, 0.970), ylim = c(0, 5000)) +
    labs(
      title = "A",
      x = "GEDI Understory Openness Index (UOI)",
      y = "Total Mammal Biomass Index (log1p scale)"
    ) +
    guides(
      size = "none",
      alpha = "none",
      fill = "none",
      shape = guide_legend(override.aes = list(size = 2, fill = "gray50")),
      color = guide_legend(override.aes = list(linewidth = 0.8)),
      linetype = guide_legend()
    ) +
    theme_pnas(base_size = 7.5) +
    theme(
      legend.position = "right",
      aspect.ratio = 1.0,
      plot.title = element_text(face = "bold", size = 8.5, margin = margin(b = 6, t = 4)),
      axis.title.x = element_text(margin = margin(t = 4)),
      axis.title.y = element_text(margin = margin(r = 4)),
      plot.margin = margin(t = 1.5, r = 4, b = 1.5, l = 4, unit = "pt"),
      legend.title = element_text(size = 7.0, face = "bold"),
      legend.text = element_text(size = 6.5),
      legend.key.width = unit(0.4, "cm")
    )

  # --- Panel B: Standard AICc Model Selection (Top 15) ---
  plot_df_aic <- f2_selection %>%
    arrange(Full_AICc) %>%
    head(15)

  plot_df_aic$plotmath_label <- sapply(1:nrow(plot_df_aic), function(i) {
    format_model_label(plot_df_aic$Model[i], f2_full_models[[plot_df_aic$Model[i]]], VAR_MAP_FRAMEWORK2)
  })

  plot_sel_df_aic <- plot_df_aic %>%
    mutate(
      CleanName = sub("^M2\\.[0-9\\.]+[a-z_]*: ", "", Model),
      CleanName = factor(CleanName, levels = rev(CleanName))
    )

  ordered_exprs_aic <- plot_sel_df_aic$plotmath_label[match(levels(plot_sel_df_aic$CleanName), plot_sel_df_aic$CleanName)]
  parsed_labels_aic <- parse(text = ordered_exprs_aic)

  plot_sel_df_aic$Full_DevExpl_Pct <- plot_sel_df_aic$Full_DevExpl * 100

  p_b <- plot_model_selection_bars(
    plot_df = plot_sel_df_aic,
    x_var = "Full_DevExpl_Pct",
    fill_var = "delta_AICc",
    fill_label = "Delta AICc",
    x_label = "Model Deviance Explained (%)",
    plot_title = "B",
    parsed_labels = parsed_labels_aic
  ) + labs(y = NULL) + theme(aspect.ratio = 0.25, plot.margin = margin(t = 1.5, r = 4, b = 1.5, l = 4, unit = "pt"))

  # --- Panel C: Generalizability Model Selection (Top 15) ---
  plot_df_lobo <- f2_selection %>%
    head(15)

  plot_df_lobo$plotmath_label <- sapply(1:nrow(plot_df_lobo), function(i) {
    format_model_label(plot_df_lobo$Model[i], f2_full_models[[plot_df_lobo$Model[i]]], VAR_MAP_FRAMEWORK2)
  })

  plot_sel_df_lobo <- plot_df_lobo %>%
    mutate(
      CleanName = sub("^M2\\.[0-9\\.]+[a-z_]*: ", "", Model),
      CleanName = factor(CleanName, levels = rev(CleanName))
    )

  ordered_exprs_lobo <- plot_sel_df_lobo$plotmath_label[match(levels(plot_sel_df_lobo$CleanName), plot_sel_df_lobo$CleanName)]
  parsed_labels_lobo <- parse(text = ordered_exprs_lobo)

  plot_sel_df_lobo$Full_DevExpl_Pct <- plot_sel_df_lobo$Full_DevExpl * 100

  p_c <- plot_model_selection_bars(
    plot_df = plot_sel_df_lobo,
    x_var = "Full_DevExpl_Pct",
    fill_var = "OOS_MAE_log",
    fill_label = "OOS MAE\n(log1p)",
    x_label = "Model Deviance Explained (%)",
    plot_title = "C",
    parsed_labels = parsed_labels_lobo
  ) + labs(y = NULL) + theme(aspect.ratio = 0.25, plot.margin = margin(t = 1.5, r = 4, b = 1.5, l = 4, unit = "pt"))

  # --- Assemble Figure 4 via Patchwork ---
  fig4_final <- p_a + p_b + p_c + plot_layout(ncol = 1)

  save_pnas(plot = fig4_final, filename = "figures/figure4.png", type = "double", height_cm = 20.7)
  save_pnas(plot = fig4_final, filename = "figures/figure4.pdf", type = "double", height_cm = 20.7)

  cat("✓ Saved Figure 4 to figures/figure4.png and .pdf\n\n")
}


# =============================================================================
# EXECUTE
# =============================================================================

generate_figure3()
generate_figure4()

cat("=== All Framework Figures Generated Successfully ===\n")
