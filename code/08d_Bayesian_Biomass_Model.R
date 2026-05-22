# =============================================================================
# 08d_Bayesian_Biomass_Model.R
#
# Fits a mathematically rigorous Bayesian Hurdle Gamma model with ecologically
# justified regularizing priors to model standing mammal biomass index (n < 100).
# Extrapolates expected biomass and spatial uncertainty (95% Credible Interval)
# across the Congo and Amazon study regions.
#
# Generates a premium 4-panel (2x2) PNAS-style publication figure:
#   outputs/congo_ct_gee_biomass_predictions_and_uncertainty_bayesian.png
# =============================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(cowplot)
library(dplyr)
library(readr)
library(brms)
library(scales)

cat("=== Formally Fitting Bayesian Hurdle Gamma Model (N = 22) ===\n\n")

# --- 1. Load Theme, Data & Model ---------------------------------------------
source("code/functions/theme_pnas.R")

# Ensure outputs directory exists
dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

# Load the pre-joined dataset
joined_data <- read_csv("outputs/congo_camera_trap_gee_joined.csv", show_col_types = FALSE)

# Normalize weights so they sum to N (mean = 1) for proper Bayesian statistical scale
joined_data$w_uoi <- joined_data$w_uoi / mean(joined_data$w_uoi)

# Center the GEDI UOI predictor at its mean (0.95) to stabilize intercept priors
# and improve Hamiltonian Monte Carlo (HMC) chain mixing.
joined_data$uoi_c <- joined_data$uoi - 0.95

# Formulate the Hurdle Gamma model:
#   1. Logistic regression for probability of hurdle zero (hu ~ 1)
#   2. Gamma regression with a log link for positive continuous biomass (B_H_index ~ uoi_c)
# Both weighted by normalized effort-uncertainty weights.
b_formula <- bf(
  B_H_index | weights(w_uoi) ~ uoi_c,
  hu ~ 1
)

# Ecologically and statistically justified regularizing priors:
#   - Intercept: normal(5, 2). At uoi_c = 0 (mean UOI of 0.95), expected biomass index
#     is exp(5) ~ 150. Allows broad biological variation from exp(1)~3 to exp(9)~8100.
#   - Slope: normal(0, 30). Strongly regularizing. Center of 0 represents no relationship.
#     A standard deviation of 30 allows steep slopes up to 60-80 if strongly supported,
#     but regularizes extreme slopes (>100) to structurally prevent explosive predictions.
#   - Hurdle Intercept: logistic(0, 1). Weakly informative prior on zero probability.
#   - Shape: gamma(0.01, 0.01). Standard conjugate prior for Gamma shape.
b_priors <- c(
  prior(normal(5, 2), class = "Intercept"),
  prior(normal(0, 30), class = "b", coef = "uoi_c"),
  prior(logistic(0, 1), class = "Intercept", dpar = "hu"),
  prior(gamma(0.01, 0.01), class = "shape")
)

cat("Compiling and sampling from Bayesian Hurdle Gamma model...\n")
fit_bayesian <- brm(
  formula = b_formula,
  data = joined_data,
  family = hurdle_gamma(link = "log"),
  prior = b_priors,
  chains = 4,
  cores = 4,
  iter = 5000,
  warmup = 2500,
  control = list(adapt_delta = 0.98),
  backend = "rstan",
  seed = 42
)

cat("\nBayesian Model Estimation Summary:\n")
print(summary(fit_bayesian))
cat("\n")

# Save model object
saveRDS(fit_bayesian, "outputs/bayesian_biomass_model.rds")

# --- 2. Load Spatial GIS Layers -----------------------------------------------
geojson_path <- "outputs/camera_traps_robust_buffered_mcps.geojson"
if (!file.exists(geojson_path)) {
  stop("Camera trap buffered MCPs missing: ", geojson_path)
}
mcps <- terra::vect(geojson_path)
mcps_congo <- mcps[mcps$region == "Congo", ]
mcps_amazon <- mcps[mcps$region == "Amazon", ]

r_congo_path <- "outputs/EOdata/analysis_stack_5000_Congo.tif"
r_amazon_path <- "outputs/EOdata/analysis_stack_5000_Amazon.tif"

r_congo <- rast(r_congo_path)
r_amazon <- rast(r_amazon_path)

names(r_congo) <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                    "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
names(r_amazon) <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                     "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")

# Crop rasters to study extents with 1.0 degree padding
study_extent_congo <- ext(
  xmin(ext(mcps_congo)) - 1.0, xmax(ext(mcps_congo)) + 1.0,
  ymin(ext(mcps_congo)) - 1.0, ymax(ext(mcps_congo)) + 1.0
)
study_extent_amazon <- ext(
  xmin(ext(mcps_amazon)) - 1.0, xmax(ext(mcps_amazon)) + 1.0,
  ymin(ext(mcps_amazon)) - 1.0, ymax(ext(mcps_amazon)) + 1.0
)

r_congo_cropped <- crop(r_congo, study_extent_congo)
r_amazon_cropped <- crop(r_amazon, study_extent_amazon)

# --- 3. Run Spatial Predictions -----------------------------------------------
cat("Predicting expected biomass and uncertainty on Congo raster...\n")
congo_cells <- as.data.frame(r_congo_cropped[["uoi"]], cells = TRUE, xy = TRUE, na.rm = TRUE)
names(congo_cells)[names(congo_cells) == "uoi"] <- "uoi"
# Create centered predictor variable matching model training
congo_cells$uoi_c <- congo_cells$uoi - 0.95

# Predict expected value of the posterior distribution (mu) and credible bounds
fitted_congo <- fitted(fit_bayesian, newdata = congo_cells)

congo_cells$pred <- fitted_congo[, "Estimate"]
congo_cells$lower <- fitted_congo[, "Q2.5"]
congo_cells$upper <- fitted_congo[, "Q97.5"]
congo_cells$ci_range <- congo_cells$upper - congo_cells$lower

# Assign predicted values to SpatRasters
r_pred_congo <- rast(r_congo_cropped[["uoi"]])
names(r_pred_congo) <- "pred"
values(r_pred_congo) <- NA
r_pred_congo[congo_cells$cell] <- as.vector(congo_cells$pred)

r_unc_congo <- rast(r_congo_cropped[["uoi"]])
names(r_unc_congo) <- "ci_range"
values(r_unc_congo) <- NA
r_unc_congo[congo_cells$cell] <- as.vector(congo_cells$ci_range)

cat("Predicting expected biomass and uncertainty on Amazon raster...\n")
amazon_cells <- as.data.frame(r_amazon_cropped[["uoi"]], cells = TRUE, xy = TRUE, na.rm = TRUE)
names(amazon_cells)[names(amazon_cells) == "uoi"] <- "uoi"
amazon_cells$uoi_c <- amazon_cells$uoi - 0.95

fitted_amazon <- fitted(fit_bayesian, newdata = amazon_cells)

amazon_cells$pred <- fitted_amazon[, "Estimate"]
amazon_cells$lower <- fitted_amazon[, "Q2.5"]
amazon_cells$upper <- fitted_amazon[, "Q97.5"]
amazon_cells$ci_range <- amazon_cells$upper - amazon_cells$lower

# Assign predicted values to SpatRasters
r_pred_amazon <- rast(r_amazon_cropped[["uoi"]])
names(r_pred_amazon) <- "pred"
values(r_pred_amazon) <- NA
r_pred_amazon[amazon_cells$cell] <- as.vector(amazon_cells$pred)

r_unc_amazon <- rast(r_amazon_cropped[["uoi"]])
names(r_unc_amazon) <- "ci_range"
values(r_unc_amazon) <- NA
r_unc_amazon[amazon_cells$cell] <- as.vector(amazon_cells$ci_range)

# Cat distribution metrics to verify regularization effect
cat("\nCongo Bayesian Prediction Summary (Response Scale):\n")
print(summary(congo_cells$pred))
cat("Amazon Bayesian Prediction Summary (Response Scale):\n")
print(summary(amazon_cells$pred))
cat("\n")

# --- 4. Build ggplot Panels ---------------------------------------------------
cat("Building 4-panel Bayesian predictive maps...\n")

# Retrieve font
base_family <- if (requireNamespace("showtext", quietly = TRUE)) "Helvetica Neue" else "sans"

# Define visual ceiling based on our empirical training maximum plus margin
# (Using 5000 captures all detailed gradients cleanly)
visual_cap <- 5000

# Panel A: Congo Prediction
p_pred_congo <- ggplot() +
  geom_spatraster(data = r_pred_congo) +
  scale_fill_viridis_c(
    option = "inferno",
    name = "Expected\nBiomass\nIndex",
    limits = c(0, visual_cap),
    oob = scales::squish,
    na.value = "transparent"
  ) +
  geom_spatvector(data = mcps_congo, fill = NA, color = "black", linewidth = 0.4, alpha = 0.8) +
  theme_pnas(base_size = 8) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 5.5, face = "bold"),
    legend.text = element_text(size = 5.0),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.12, "cm"),
    legend.margin = margin(l = -2, r = 0, t = 0, b = 0, unit = "pt"),
    axis.title = element_blank()
  ) +
  labs(
    title = "A. Congo Basin Expected Mammal Biomass Index",
    subtitle = "Bayesian Hurdle Gamma posterior expected value (mean)"
  )

# Panel B: Congo Uncertainty (95% Credible Interval Range)
p_unc_congo <- ggplot() +
  geom_spatraster(data = r_unc_congo) +
  scale_fill_viridis_c(
    option = "mako",
    name = "95% Cred.\nInterval\nRange",
    limits = c(0, visual_cap),
    oob = scales::squish,
    na.value = "transparent"
  ) +
  geom_spatvector(data = mcps_congo, fill = NA, color = "black", linewidth = 0.4, alpha = 0.8) +
  theme_pnas(base_size = 8) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 5.5, face = "bold"),
    legend.text = element_text(size = 5.0),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.12, "cm"),
    legend.margin = margin(l = -2, r = 0, t = 0, b = 0, unit = "pt"),
    axis.title = element_blank()
  ) +
  labs(
    title = "B. Congo Basin Prediction Uncertainty",
    subtitle = "Posterior expected width (Upper 95% CI - Lower 95% CI)"
  )

# Panel C: Amazon Prediction
p_pred_amazon <- ggplot() +
  geom_spatraster(data = r_pred_amazon) +
  scale_fill_viridis_c(
    option = "inferno",
    name = "Expected\nBiomass\nIndex",
    limits = c(0, visual_cap),
    oob = scales::squish,
    na.value = "transparent"
  ) +
  geom_spatvector(data = mcps_amazon, fill = NA, color = "black", linewidth = 0.4, alpha = 0.8) +
  theme_pnas(base_size = 8) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 5.5, face = "bold"),
    legend.text = element_text(size = 5.0),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.12, "cm"),
    legend.margin = margin(l = -2, r = 0, t = 0, b = 0, unit = "pt"),
    axis.title = element_blank()
  ) +
  labs(
    title = "C. Amazon Basin Expected Mammal Biomass Index",
    subtitle = "Bayesian Hurdle Gamma posterior expected value (mean)"
  )

# Panel D: Amazon Uncertainty (95% Credible Interval Range)
p_unc_amazon <- ggplot() +
  geom_spatraster(data = r_unc_amazon) +
  scale_fill_viridis_c(
    option = "mako",
    name = "95% Cred.\nInterval\nRange",
    limits = c(0, visual_cap),
    oob = scales::squish,
    na.value = "transparent"
  ) +
  geom_spatvector(data = mcps_amazon, fill = NA, color = "black", linewidth = 0.4, alpha = 0.8) +
  theme_pnas(base_size = 8) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 5.5, face = "bold"),
    legend.text = element_text(size = 5.0),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.12, "cm"),
    legend.margin = margin(l = -2, r = 0, t = 0, b = 0, unit = "pt"),
    axis.title = element_blank()
  ) +
  labs(
    title = "D. Amazon Basin Prediction Uncertainty",
    subtitle = "Posterior expected width (Upper 95% CI - Lower 95% CI)"
  )

# --- 5. Assemble Grid & Save --------------------------------------------------
fig_spatial_pred <- cowplot::plot_grid(
  p_pred_congo, p_unc_congo,
  p_pred_amazon, p_unc_amazon,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

save_pnas(
  plot = fig_spatial_pred,
  filename = "outputs/congo_ct_gee_biomass_predictions_and_uncertainty_bayesian.png",
  type = "double",
  height_cm = 16.0
)

cat("✓ Saved: outputs/congo_ct_gee_biomass_predictions_and_uncertainty_bayesian.png\n")
cat("=== Bayesian Spatial Prediction Mapping Complete ===\n")
