# =============================================================================
# code/05b_Predictive_Columns_Correlation.R
#
# Generates a regional faceted scatter plot correlating Column 1 (Biophysical
# Template, LOBO M2.1: UOI Only) and Column 2 (Basin-Calibrated Model, AIC
# M2.17b: UOI * Basin + Elev) standing mammal biomass predictions at 20 km.
# Computes and displays both Pearson r and Spearman rank rho coefficients.
#
# Inputs:
#   - outputs/EOdata/analysis_stack_5000_{Basin}.tif
#   - outputs/framework2_best_model.RDS
#   - outputs/framework2_best_model_aic.RDS
#
# Outputs:
#   - figures/figureS6.png (and .pdf)
# =============================================================================

library(terra)
library(ggplot2)
library(dplyr)
library(readr)

# 1. Load models
outputs_dir <- "outputs"
m_best <- readRDS(file.path(outputs_dir, "framework2_best_model.RDS"))
m_best_aic <- readRDS(file.path(outputs_dir, "framework2_best_model_aic.RDS"))

# Load global calibration parameters
sel_path <- file.path(outputs_dir, "framework2_covariate_model_selection.csv")
lobo_sel <- read_csv(sel_path, show_col_types = FALSE)
best_row <- lobo_sel %>% filter(Model == "M2.1: UOI Only")
if (nrow(best_row) == 0) best_row <- lobo_sel[1, ]
OOS_MAE_log <- best_row$OOS_MAE_log[1]

# Baseline mean log1p
source("code/functions/calibration_helpers.R")
joined_data <- extract_scale_data(5000)
mean_log_y_obs <- mean(log1p(joined_data$B_H_index))

# Define scale
scale_m <- 20000

rename_stack <- function(r) {
  if (is.null(r)) return(NULL)
  n_bands <- nlyr(r)
  if (n_bands == 12) {
    names(r) <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                  "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  } else if (n_bands == 11) {
    names(r) <- c("frip", "frip_mk_tau", "uoi", "rh98", "gedi_n",
                  "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  }
  return(r)
}

load_and_aggregate_real <- function(basin_name) {
  r_5000_path <- file.path(outputs_dir, "EOdata", sprintf("analysis_stack_5000_%s.tif", basin_name))
  if (!file.exists(r_5000_path)) {
    warning(sprintf("Real 5,000m stack for %s missing, trying synthetic fallback...", basin_name))
    r_synth_path <- file.path(outputs_dir, "synthetic_EOdata", sprintf("analysis_stack_%d_%s.tif", scale_m, basin_name))
    if (!file.exists(r_synth_path)) return(NULL)
    r <- rast(r_synth_path)
    return(rename_stack(r))
  }
  
  cat(sprintf("Loading real 5,000m stack for %s...\n", basin_name))
  r_5000 <- rast(r_5000_path)
  r_5000 <- rename_stack(r_5000)
  
  cat(sprintf("Aggregating real 5,000m stack for %s by factor of 4 to %d m...\n", basin_name, scale_m))
  r_agg <- terra::aggregate(r_5000, fact = 4, fun = "mean", na.rm = TRUE)
  return(r_agg)
}

r_congo <- load_and_aggregate_real("Congo")
r_amazon <- load_and_aggregate_real("Amazon")
r_seasia <- load_and_aggregate_real("SE_Asia")

# Helper to predict any model on a SpatRaster stack and dynamically inject required category covariates
predict_model_on_raster <- function(model, r_stack, basin_name) {
  model_vars <- all.vars(formula(model))
  raster_vars <- intersect(model_vars, names(r_stack))
  df <- as.data.frame(r_stack[[raster_vars]], cells = TRUE, xy = TRUE, na.rm = TRUE)
  if (nrow(df) == 0) return(NULL)
  
  # Inject categories if expected by the model formulas
  if ("basin" %in% model_vars) {
    df$basin <- factor(basin_name, levels = c("Amazon", "Congo", "SE_Asia"))
  }
  if ("elephant_present_possible" %in% model_vars) {
    df$elephant_present_possible <- factor(ifelse(basin_name == "Amazon", "Absent", "Present"), levels = c("Absent", "Present"))
  }
  if ("elephant_present_strict" %in% model_vars) {
    df$elephant_present_strict <- factor(ifelse(basin_name == "Amazon", "Absent", "Present"), levels = c("Absent", "Present"))
  }
  
  df$pred <- predict(model, newdata = df, type = "response")
  return(df)
}

# Collect data frames
predict_basin_df <- function(r_stack, basin_name) {
  if (is.null(r_stack)) return(NULL)
  
  # Predict Column 1 (LOBO: B_H_index ~ uoi)
  df_lobo <- predict_model_on_raster(m_best, r_stack, basin_name)
  if (is.null(df_lobo) || nrow(df_lobo) == 0) return(NULL)
  df_lobo$z_lobo <- (log1p(df_lobo$pred) - mean_log_y_obs) / OOS_MAE_log
  
  # Predict Column 2 (AIC: B_H_index ~ uoi * basin + elevation)
  df_aic <- predict_model_on_raster(m_best_aic, r_stack, basin_name)
  if (is.null(df_aic) || nrow(df_aic) == 0) return(NULL)
  df_aic$z_aic <- (log1p(df_aic$pred) - mean_log_y_obs) / OOS_MAE_log
  
  # Join them by cells/coordinates to get side-by-side predictions
  res_df <- df_lobo %>%
    select(x, y, z_lobo) %>%
    inner_join(df_aic %>% select(x, y, z_aic), by = c("x", "y"))
  
  res_df$basin <- factor(basin_name, levels = c("Congo", "Amazon", "SE_Asia"))
  return(res_df[, c("basin", "z_lobo", "z_aic")])
}

congo_df <- predict_basin_df(r_congo, "Congo")
amazon_df <- predict_basin_df(r_amazon, "Amazon")
seasia_df <- if (!is.null(r_seasia)) predict_basin_df(r_seasia, "SE_Asia") else NULL

all_data <- bind_rows(congo_df, amazon_df, seasia_df)

# Calculate correlations (Pearson and Spearman)
r_overall <- cor(all_data$z_lobo, all_data$z_aic, method = "pearson", use = "complete.obs")
rho_overall <- cor(all_data$z_lobo, all_data$z_aic, method = "spearman", use = "complete.obs")

r_amazon <- cor(all_data$z_lobo[all_data$basin == "Amazon"], all_data$z_aic[all_data$basin == "Amazon"], method = "pearson", use = "complete.obs")
rho_amazon <- cor(all_data$z_lobo[all_data$basin == "Amazon"], all_data$z_aic[all_data$basin == "Amazon"], method = "spearman", use = "complete.obs")

r_congo <- cor(all_data$z_lobo[all_data$basin == "Congo"], all_data$z_aic[all_data$basin == "Congo"], method = "pearson", use = "complete.obs")
rho_congo <- cor(all_data$z_lobo[all_data$basin == "Congo"], all_data$z_aic[all_data$basin == "Congo"], method = "spearman", use = "complete.obs")

has_seasia <- !is.null(seasia_df) && nrow(seasia_df) > 0
r_seasia <- if (has_seasia) cor(all_data$z_lobo[all_data$basin == "SE_Asia"], all_data$z_aic[all_data$basin == "SE_Asia"], method = "pearson", use = "complete.obs") else NA
rho_seasia <- if (has_seasia) cor(all_data$z_lobo[all_data$basin == "SE_Asia"], all_data$z_aic[all_data$basin == "SE_Asia"], method = "spearman", use = "complete.obs") else NA

# Create annotation data frame for regional correlations (placed in top-left corner)
anno_df <- data.frame(
  basin = factor(c("Congo", "Amazon", "SE_Asia"), levels = c("Congo", "Amazon", "SE_Asia")),
  z_lobo = c(-2.3, -2.3, -2.3),
  z_aic = c(2.3, 2.3, 2.3),
  label = c(
    sprintf("Pearson r = %.2f\nSpearman rho = %.2f", r_congo, rho_congo),
    sprintf("Pearson r = %.2f\nSpearman rho = %.2f", r_amazon, rho_amazon),
    sprintf("Pearson r = %.2f\nSpearman rho = %.2f", if (is.na(r_seasia)) 0.0 else r_seasia, if (is.na(rho_seasia)) 0.0 else rho_seasia)
  )
)
if (!has_seasia) {
  anno_df <- anno_df %>% filter(basin != "SE_Asia")
}

# Plot a simple, unified correlation scatter plot faceted by region
p <- ggplot(all_data, aes(x = z_lobo, y = z_aic, color = basin)) +
  geom_point(alpha = 0.25, size = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#e74c3c", se = FALSE, linewidth = 0.8) +
  geom_text(
    data = anno_df,
    aes(x = z_lobo, y = z_aic, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    color = "black",
    fontface = "bold",
    size = 2.8
  ) +
  scale_color_manual(values = c("Amazon" = "#E65100", "Congo" = "#1B5E20", "SE_Asia" = "#0D47A1")) +
  facet_wrap(~ basin, ncol = 3) +
  labs(
    title = "Correlation Between Biophysical Template & Basin-Calibrated Models",
    subtitle = sprintf("Z-score Standing Mammal Biomass predictions (at %d m predictive scale) | Overall: Pearson r = %.2f, Spearman rho = %.2f", scale_m, r_overall, rho_overall),
    x = "Column 1: Biophysical Template (Z-score)",
    y = "Column 2: Basin-Calibrated Model (Z-score)"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, color = "grey40"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "none"
  )

# Save supplementary figures
dir.create("figures", showWarnings = FALSE)
ggsave("figures/figureS4.png", plot = p, width = 18, height = 8, units = "cm", dpi = 600, bg = "white")
ggsave("figures/figureS4.pdf", plot = p, width = 18, height = 8, units = "cm", dpi = 600, bg = "white")
cat("✓ Successfully saved supplementary figure to figures/figureS4.png and .pdf\n")


