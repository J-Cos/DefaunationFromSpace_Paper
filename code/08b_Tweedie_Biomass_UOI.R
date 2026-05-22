# =============================================================================
# 08b_Tweedie_Biomass_UOI.R
#
# Alternative version of the biomass ~ UOI relationship figure using Tweedie
# GLMs instead of weighted linear regression.
#
# Rationale:
#   - Biomass is non-negative and right-skewed with exact zeros
#   - Tweedie(p ~ 1.5, log link) naturally handles zero-inflated positive
#     continuous data without arbitrary transformations
#   - Log link ensures predictions are always >= 0
#   - GAM smooth s(UOI) allows nonlinearity detection with small N
#
# Outputs:
#   outputs/congo_ct_gee_biomass_uoi_tweedie.png  (3-panel: total, >50kg, >100kg)
# =============================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(cowplot)
library(dplyr)
library(readr)
library(mgcv)      # GAM with Tweedie family

cat("=== Tweedie GLM: Biomass ~ UOI Alternative Modelling ===\n\n")

# --- 1. Load Theme & Data ---------------------------------------------------
source("code/functions/theme_pnas.R")

# Load the pre-joined dataset from the main pipeline
joined_data <- read_csv("outputs/congo_camera_trap_gee_joined.csv", show_col_types = FALSE)
cat(sprintf("Loaded %d clusters with valid UOI + biomass data.\n", nrow(joined_data)))

# Normalize weights so they sum to N (mean = 1) to ensure proper statistical scale.
# Otherwise, raw precision weights are treated as frequency weights in GAM/GLM,
# which inflates sample size, underestimating standard errors and p-values.
joined_data$w_uoi <- joined_data$w_uoi / mean(joined_data$w_uoi)

# Recap zero structure
cat(sprintf("  B_H_index  zeros: %d / %d\n", sum(joined_data$B_H_index == 0), nrow(joined_data)))
cat(sprintf("  B_H_gt50   zeros: %d / %d\n", sum(joined_data$B_H_gt50 == 0), nrow(joined_data)))
cat(sprintf("  B_H_gt100  zeros: %d / %d\n", sum(joined_data$B_H_gt100 == 0), nrow(joined_data)))

# Retrieve base font family from theme_pnas
base_family <- if (requireNamespace("showtext", quietly = TRUE)) "Helvetica Neue" else "sans"


# --- 2. Fit Tweedie GLMs & Generate Prediction Curves -----------------------

# We fit two models per response variable:
#   1. Tweedie GLM (linear on log scale): gam(y ~ UOI, family = tw())
#   2. Tweedie GAM (smooth on log scale): gam(y ~ s(UOI, k=4), family = tw())
# Compare via AIC. k=4 limits wiggliness appropriate for N~22.

fit_tweedie_models <- function(df, y_col, weight_col) {
  # Tweedie needs response > 0 OR exactly 0. Our data satisfies this.
  y <- df[[y_col]]
  w <- df[[weight_col]]
  
  # Model 1: Linear on log scale (Tweedie GLM)
  m_linear <- gam(
    as.formula(paste(y_col, "~ uoi")),
    family = tw(),
    weights = w,
    data = df,
    method = "REML"
  )
  
  # Model 2: Smooth nonlinear (Tweedie GAM, k=4 for small N)
  m_smooth <- tryCatch({
    gam(
      as.formula(paste(y_col, "~ s(uoi, k = 4)")),
      family = tw(),
      weights = w,
      data = df,
      method = "REML"
    )
  }, error = function(e) NULL)
  
  list(linear = m_linear, smooth = m_smooth)
}

# Generate prediction data frame with confidence intervals
predict_tweedie <- function(model, df, n_pred = 200) {
  uoi_seq <- seq(min(df$uoi) - 0.002, max(df$uoi) + 0.002, length.out = n_pred)
  newdata <- data.frame(uoi = uoi_seq)
  
  # Predict on link (log) scale with SE, then transform to response scale
  pred <- predict(model, newdata = newdata, type = "link", se.fit = TRUE)
  
  newdata$fit <- exp(pred$fit)  # Back-transform from log link
  newdata$lower <- exp(pred$fit - 1.96 * pred$se.fit)
  newdata$upper <- exp(pred$fit + 1.96 * pred$se.fit)
  
  newdata
}

# Extract model summary statistics
summarise_tweedie <- function(model, y_col) {
  s <- summary(model)
  
  # For parametric (linear) term, extract p-value for UOI slope
  if ("uoi" %in% rownames(s$p.table)) {
    p_val <- s$p.table["uoi", "Pr(>|t|)"]
    coef_uoi <- s$p.table["uoi", "Estimate"]
    # Percentage change in response per 0.01 UOI increase
    pct_change <- (exp(coef_uoi * 0.01) - 1) * 100
    type <- "GLM"
  } else {
    # Smooth term
    p_val <- s$s.table[1, "p-value"]
    edf <- s$s.table[1, "edf"]
    pct_change <- NA
    type <- sprintf("GAM (edf=%.1f)", edf)
  }
  
  dev_expl <- s$dev.expl * 100
  p_est <- s$family$getTheta(TRUE)  # Tweedie power parameter
  
  list(
    type = type,
    dev_explained = dev_expl,
    p_value = p_val,
    tweedie_p = p_est,
    pct_change = pct_change,
    aic = AIC(model)
  )
}


# --- 3. Fit all models -------------------------------------------------------

responses <- list(
  list(y = "B_H_index",  label = "Total Biomass Index",                    title_letter = "A"),
  list(y = "B_H_gt50",   label = "Large Fauna Biomass (>50 kg)",           title_letter = "B"),
  list(y = "B_H_gt100",  label = "Megafauna Biomass (>100 kg)",            title_letter = "C")
)

cat("\n--- Tweedie GLM Results ---\n\n")

all_models <- list()
all_preds <- list()
all_stats <- list()

for (resp in responses) {
  cat(sprintf("  %s. %s (%s):\n", resp$title_letter, resp$label, resp$y))
  
  models <- fit_tweedie_models(joined_data, resp$y, "w_uoi")
  
  # Summarise linear model
  stats_lin <- summarise_tweedie(models$linear, resp$y)
  cat(sprintf("    GLM: deviance explained = %.1f%%, p = %.4f, Tweedie p = %.2f\n",
              stats_lin$dev_explained, stats_lin$p_value, stats_lin$tweedie_p))
  if (!is.na(stats_lin$pct_change)) {
    cat(sprintf("    → %.1f%% change in biomass per 0.01 UOI increase\n", stats_lin$pct_change))
  }
  
  # Summarise smooth model if it converged
  if (!is.null(models$smooth)) {
    stats_smo <- summarise_tweedie(models$smooth, resp$y)
    cat(sprintf("    GAM: deviance explained = %.1f%%, p = %.4f, edf in label\n",
                stats_smo$dev_explained, stats_smo$p_value))
    cat(sprintf("    AIC comparison: GLM = %.1f, GAM = %.1f → %s preferred\n",
                stats_lin$aic, stats_smo$aic,
                if (stats_lin$aic <= stats_smo$aic) "GLM (simpler)" else "GAM"))
    
    # Use whichever has lower AIC
    if (stats_smo$aic < stats_lin$aic - 2) {
      best_model <- models$smooth
      best_stats <- stats_smo
      cat("    ★ Using GAM (>2 AIC improvement)\n")
    } else {
      best_model <- models$linear
      best_stats <- stats_lin
      cat("    ★ Using GLM (parsimony / <2 AIC difference)\n")
    }
  } else {
    best_model <- models$linear
    best_stats <- stats_lin
    cat("    GAM did not converge, using GLM.\n")
  }
  
  # Generate predictions
  pred_df <- predict_tweedie(best_model, joined_data)
  
  all_models[[resp$y]] <- best_model
  all_preds[[resp$y]] <- pred_df
  all_stats[[resp$y]] <- best_stats
  
  cat("\n")
}


# --- 4. Build Publication Figure (3-panel, 1 row) ----------------------------

cat("Generating Tweedie GLM figure (3-panel)...\n")

make_tweedie_panel <- function(df, pred_df, stats, y_col, y_label, title, line_color) {
  
  # Format statistics annotation
  p_text <- if (stats$p_value < 0.001) "p < 0.001" else sprintf("p = %.3f", stats$p_value)
  stat_label <- sprintf(
    "Tweedie %s\nDev. expl. = %.0f%%\n%s (N = %d)",
    stats$type, stats$dev_explained, p_text, nrow(df)
  )
  
  ggplot() +
    # Confidence ribbon from Tweedie predictions (always >= 0)
    geom_ribbon(
      data = pred_df,
      aes(x = uoi, ymin = lower, ymax = upper),
      fill = line_color,
      alpha = 0.12
    ) +
    # Fitted curve
    geom_line(
      data = pred_df,
      aes(x = uoi, y = fit),
      color = line_color,
      linewidth = 0.7
    ) +
    # Data points
    geom_point(
      data = df,
      aes(x = uoi, y = !!sym(y_col), size = trap_days,
          fill = uoi_se, shape = basin),
      color = "black",
      stroke = 0.3,
      alpha = 0.85
    ) +
    scale_shape_manual(
      values = c("Congo" = 21, "Amazon" = 24),
      name = "Basin"
    ) +
    scale_size_continuous(
      range = c(1.5, 6.0),
      guide = "none"
    ) +
    scale_fill_viridis_c(option = "magma", name = "UOI\nSpatial SE") +
    # Annotation
    annotate(
      "text",
      x = Inf, y = Inf,
      label = stat_label,
      hjust = 1.1, vjust = 1.2,
      family = base_family,
      size = 2.2,
      fontface = "bold",
      color = "grey20"
    ) +
    # y-axis starts at zero
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    theme_pnas(base_size = 8) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 5.5, face = "bold"),
      legend.text = element_text(size = 5.0),
      legend.key.height = unit(0.25, "cm"),
      legend.key.width = unit(0.12, "cm"),
      legend.margin = margin(l = -2, r = 0, t = 0, b = 0, unit = "pt")
    ) +
    labs(
      title = title,
      x = "GEDI Understory Openness Index (UOI)",
      y = y_label
    )
}

# Build 3 panels
panels <- list()

panel_configs <- list(
  list(y = "B_H_index", label = "Total Biomass Index",          title = "A. Total Biomass vs. Understory Openness",     color = "#311B92"),
  list(y = "B_H_gt50",  label = "Large Fauna Biomass (>50 kg)", title = "B. Large Fauna Biomass vs. Understory Openness", color = "#4527A0"),
  list(y = "B_H_gt100", label = "Megafauna Biomass (>100 kg)",  title = "C. Megafauna Biomass vs. Understory Openness",  color = "#1565C0")
)

for (cfg in panel_configs) {
  panels[[cfg$y]] <- make_tweedie_panel(
    df = joined_data,
    pred_df = all_preds[[cfg$y]],
    stats = all_stats[[cfg$y]],
    y_col = cfg$y,
    y_label = cfg$label,
    title = cfg$title,
    line_color = cfg$color
  )
}

# Assemble 1x3 layout
fig_tweedie <- cowplot::plot_grid(
  panels[["B_H_index"]],
  panels[["B_H_gt50"]],
  panels[["B_H_gt100"]],
  ncol = 3,
  align = "vh"
)

save_pnas(
  plot = fig_tweedie,
  filename = "outputs/congo_ct_gee_biomass_uoi_tweedie.png",
  type = "double",
  height_cm = 8.0
)
cat("✓ Saved: outputs/congo_ct_gee_biomass_uoi_tweedie.png\n")


# --- 5. Also make a comparison figure: Linear vs Tweedie side by side --------
cat("\nGenerating comparison figure (Linear vs Tweedie)...\n")

# Refit linear models for overlay
make_comparison_panel <- function(df, pred_df_tw, y_col, y_label, title, line_color) {
  
  # Fit weighted linear model for comparison
  w_fit <- lm(as.formula(paste(y_col, "~ uoi")), weights = df$w_uoi, data = df)
  
  # Generate linear predictions
  uoi_seq <- seq(min(df$uoi) - 0.002, max(df$uoi) + 0.002, length.out = 200)
  lm_pred <- predict(w_fit, newdata = data.frame(uoi = uoi_seq), 
                      interval = "confidence", level = 0.95)
  lm_df <- data.frame(uoi = uoi_seq, fit = lm_pred[,1], 
                       lower = lm_pred[,2], upper = lm_pred[,3])
  
  ggplot() +
    # Linear model (red, dashed)
    geom_ribbon(
      data = lm_df,
      aes(x = uoi, ymin = lower, ymax = upper),
      fill = "#C62828",
      alpha = 0.08
    ) +
    geom_line(
      data = lm_df,
      aes(x = uoi, y = fit),
      color = "#C62828",
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    # Tweedie model (solid)
    geom_ribbon(
      data = pred_df_tw,
      aes(x = uoi, ymin = lower, ymax = upper),
      fill = line_color,
      alpha = 0.12
    ) +
    geom_line(
      data = pred_df_tw,
      aes(x = uoi, y = fit),
      color = line_color,
      linewidth = 0.7
    ) +
    # Zero reference line
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60", linetype = "dotted") +
    # Data points
    geom_point(
      data = df,
      aes(x = uoi, y = !!sym(y_col), size = trap_days,
          fill = uoi_se, shape = basin),
      color = "black",
      stroke = 0.3,
      alpha = 0.85
    ) +
    scale_shape_manual(
      values = c("Congo" = 21, "Amazon" = 24),
      name = "Basin"
    ) +
    scale_size_continuous(range = c(1.5, 6.0), guide = "none") +
    scale_fill_viridis_c(option = "magma", name = "UOI\nSpatial SE") +
    annotate(
      "text", x = -Inf, y = Inf,
      label = "—— Tweedie GLM\n- - - Linear",
      hjust = -0.05, vjust = 1.2,
      family = base_family,
      size = 2.0,
      color = "grey40"
    ) +
    theme_pnas(base_size = 8) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 5.5, face = "bold"),
      legend.text = element_text(size = 5.0),
      legend.key.height = unit(0.25, "cm"),
      legend.key.width = unit(0.12, "cm"),
      legend.margin = margin(l = -2, r = 0, t = 0, b = 0, unit = "pt")
    ) +
    labs(
      title = title,
      x = "GEDI Understory Openness Index (UOI)",
      y = y_label
    )
}

comp_panels <- list()
for (cfg in panel_configs) {
  comp_panels[[cfg$y]] <- make_comparison_panel(
    df = joined_data,
    pred_df_tw = all_preds[[cfg$y]],
    y_col = cfg$y,
    y_label = cfg$label,
    title = cfg$title,
    line_color = cfg$color
  )
}

fig_comparison <- cowplot::plot_grid(
  comp_panels[["B_H_index"]],
  comp_panels[["B_H_gt50"]],
  comp_panels[["B_H_gt100"]],
  ncol = 3,
  align = "vh"
)

save_pnas(
  plot = fig_comparison,
  filename = "outputs/congo_ct_gee_biomass_uoi_tweedie_vs_linear.png",
  type = "double",
  height_cm = 8.0
)
cat("✓ Saved: outputs/congo_ct_gee_biomass_uoi_tweedie_vs_linear.png\n")

cat("\n=== Tweedie GLM analysis complete ===\n")
