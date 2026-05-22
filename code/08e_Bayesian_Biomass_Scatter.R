# =============================================================================
# 08e_Bayesian_Biomass_Scatter.R
#
# Generates a publication-quality scatter plot showing the Bayesian Hurdle Gamma
# model fit (posterior expected value and 95% Credible Interval) against the
# empirical camera trap data points.
#
# Saves the figure as:
#   outputs/congo_ct_gee_biomass_uoi_bayesian_fit.png
# =============================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(brms)
library(cowplot)
library(scales)

cat("=== Generating Bayesian Model vs. Empirical Data Scatter Plot ===\n\n")

# --- 1. Load Theme, Data & Model ---------------------------------------------
source("code/functions/theme_pnas.R")

# Load empirical joined data
joined_data <- read_csv("outputs/congo_camera_trap_gee_joined.csv", show_col_types = FALSE)

# Load pre-compiled Bayesian Hurdle Gamma model fit
if (!file.exists("outputs/bayesian_biomass_model.rds")) {
  stop("Bayesian model RDS file not found. Please run code/08d_Bayesian_Biomass_Model.R first.")
}
fit_bayesian <- readRDS("outputs/bayesian_biomass_model.rds")

# --- 2. Generate Model Fit Predictions ----------------------------------------
# Spans the empirical range of UOI in the dataset
uoi_seq <- seq(from = 0.915, to = 0.975, length.out = 200)
pred_data <- data.frame(
  uoi = uoi_seq,
  uoi_c = uoi_seq - 0.95
)

cat("Extracting posterior expected values and credible intervals...\n")
fitted_draws <- fitted(fit_bayesian, newdata = pred_data)

pred_data$Estimate <- fitted_draws[, "Estimate"]
pred_data$Lower <- fitted_draws[, "Q2.5"]
pred_data$Upper <- fitted_draws[, "Q97.5"]

# --- 3. Build Plots -----------------------------------------------------------

uoi_label <- "GEDI Understory Openness Index (UOI)"
biomass_label <- "Mammal Standing Biomass Index"

# Set up the base plot structure
base_plot <- function(log_scale = FALSE) {
  if (log_scale) {
    # For log scale, we add +1 to biomass to display zero-biomass points (avoiding log(0))
    p_log <- ggplot() +
      geom_ribbon(
        data = pred_data,
        aes(x = uoi, ymin = Lower + 1, ymax = Upper + 1),
        fill = "#2E7D32",
        alpha = 0.15
      ) +
      geom_line(
        data = pred_data,
        aes(x = uoi, y = Estimate + 1),
        color = "#1B5E20",
        linewidth = 0.8
      ) +
      geom_point(
        data = joined_data,
        aes(
          x = uoi,
          y = B_H_index + 1,
          shape = region,
          color = region,
          size = trap_days
        ),
        alpha = 0.85,
        stroke = 0.4
      ) +
      scale_shape_manual(
        values = c(Congo = 21, Amazon = 24),
        name = "Basin"
      ) +
      scale_color_manual(
        values = pal_basin,
        name = "Basin"
      ) +
      scale_size_continuous(
        name = "Effort (Trap-days)",
        breaks = c(100, 1000, 5000, 10000, 20000),
        range = c(1.5, 5.0)
      ) +
      scale_y_log10(
        labels = trans_format("log10", math_format(10^.x)),
        breaks = c(1, 10, 100, 1000, 5000)
      ) +
      labs(
        x = uoi_label,
        y = paste(biomass_label, "(Log Scale, +1)"),
        title = "B. Bayesian Model Fit (Log Scale)",
        subtitle = "Expected value & 95% Credible Interval"
      ) +
      theme_pnas(base_size = 8) +
      theme(
        legend.position = "none"
      )
    return(p_log)
  } else {
    p_lin <- ggplot() +
      geom_ribbon(
        data = pred_data,
        aes(x = uoi, ymin = Lower, ymax = Upper),
        fill = "#2E7D32",
        alpha = 0.15
      ) +
      geom_line(
        data = pred_data,
        aes(x = uoi, y = Estimate),
        color = "#1B5E20",
        linewidth = 0.8
      ) +
      geom_point(
        data = joined_data,
        aes(
          x = uoi,
          y = B_H_index,
          shape = region,
          color = region,
          size = trap_days
        ),
        alpha = 0.85,
        stroke = 0.4
      ) +
      scale_shape_manual(
        values = c(Congo = 21, Amazon = 24),
        name = "Basin"
      ) +
      scale_color_manual(
        values = pal_basin,
        name = "Basin"
      ) +
      scale_size_continuous(
        name = "Effort (Trap-days)",
        breaks = c(100, 1000, 5000, 10000, 20000),
        range = c(1.5, 5.0)
      ) +
      scale_y_continuous(
        labels = scales::comma_format(),
        limits = c(0, 5000),
        oob = scales::squish
      ) +
      labs(
        x = uoi_label,
        y = biomass_label,
        title = "A. Bayesian Model Fit (Linear Scale)",
        subtitle = "Expected value & 95% Credible Interval"
      ) +
      theme_pnas(base_size = 8) +
      theme(
        legend.position = "none"
      )
    return(p_lin)
  }
}

p_lin <- base_plot(log_scale = FALSE)
p_log <- base_plot(log_scale = TRUE)

# --- 4. Extract Legend --------------------------------------------------------
p_legend_obj <- ggplot(joined_data) +
  geom_point(
    aes(
      x = uoi,
      y = B_H_index,
      shape = region,
      color = region,
      size = trap_days
    ),
    alpha = 0.85,
    stroke = 0.4
  ) +
  scale_shape_manual(
    values = c(Congo = 21, Amazon = 24),
    name = "Basin"
  ) +
  scale_color_manual(
    values = pal_basin,
    name = "Basin"
  ) +
  scale_size_continuous(
    name = "Effort (Trap-days)",
    breaks = c(100, 1000, 5000, 10000, 20000),
    range = c(1.5, 5.0)
  ) +
  theme_pnas(base_size = 8) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6.5)
  )

shared_legend <- cowplot::get_legend(p_legend_obj)

# --- 5. Assemble Grid & Save --------------------------------------------------
fig_panels <- cowplot::plot_grid(
  p_lin, p_log,
  ncol = 2,
  align = "h"
)

fig_final <- cowplot::plot_grid(
  fig_panels,
  shared_legend,
  ncol = 1,
  rel_heights = c(1.0, 0.1)
)

save_pnas(
  plot = fig_final,
  filename = "outputs/congo_ct_gee_biomass_uoi_bayesian_fit.png",
  type = "double",
  height_cm = 8.5
)

cat("✓ Saved: outputs/congo_ct_gee_biomass_uoi_bayesian_fit.png\n")
cat("=== Scatter Plot Generation Complete ===\n")
