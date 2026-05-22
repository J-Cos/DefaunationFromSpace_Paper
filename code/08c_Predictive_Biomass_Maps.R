# =============================================================================
# 08c_Predictive_Biomass_Maps.R
#
# Generates landscape-scale predictive maps of standing mammal biomass index
# across the Congo and Amazon study areas using the best Tweedie GLM.
# Fits a weight-normalized Tweedie model and estimates:
#   1. Expected Biomass Index (fit)
#   2. 95% Confidence Interval Range (width) as our spatial uncertainty metric.
#
# Saves a 4-panel premium PNAS-style publication figure:
#   outputs/congo_ct_gee_biomass_predictions_and_uncertainty.png
# =============================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(cowplot)
library(dplyr)
library(readr)
library(mgcv)

cat("=== Spatial Biomass Prediction & Uncertainty Mapping ===\n\n")

# --- 1. Load Theme, Data & Model ---------------------------------------------
source("code/functions/theme_pnas.R")

# Ensure outputs directory exists
dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

# Load the pre-joined dataset from the main pipeline
joined_data <- read_csv("outputs/congo_camera_trap_gee_joined.csv", show_col_types = FALSE)

# Normalize weights so they sum to N (mean = 1) to ensure proper statistical scale
joined_data$w_uoi <- joined_data$w_uoi / mean(joined_data$w_uoi)

# Fit the most predictive Tweedie GLM (Total Biomass Index vs UOI, deviance explained = 33.1%)
m_best <- gam(
  B_H_index ~ uoi,
  family = tw(),
  weights = w_uoi,
  data = joined_data,
  method = "REML"
)

cat("Fitted Best Tweedie GLM (Total Biomass Index ~ UOI):\n")
print(summary(m_best))
cat("\n")

# --- 2. Load MCPs and GEE Raster Stacks ---------------------------------------
geojson_path <- "outputs/camera_traps_robust_buffered_mcps.geojson"
if (!file.exists(geojson_path)) {
  stop("Python robust cluster buffered MCPs GeoJSON missing: ", geojson_path, 
       "\nPlease run python3 code/visualise_camera_traps.py first to generate it.")
}
mcps <- terra::vect(geojson_path)
mcps_congo <- mcps[mcps$region == "Congo", ]
mcps_amazon <- mcps[mcps$region == "Amazon", ]

r_congo_path <- "outputs/EOdata/analysis_stack_5000_Congo.tif"
r_amazon_path <- "outputs/EOdata/analysis_stack_5000_Amazon.tif"

if (!file.exists(r_congo_path)) stop("GEE 5000m Congo GeoTIFF missing: ", r_congo_path)
if (!file.exists(r_amazon_path)) stop("GEE 5000m Amazon GeoTIFF missing: ", r_amazon_path)

r_congo <- rast(r_congo_path)
r_amazon <- rast(r_amazon_path)

names(r_congo) <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                    "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
names(r_amazon) <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                     "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")

# Crop rasters using study extents with 1.0 degree padding
study_extent_congo <- ext(mcps_congo)
study_extent_congo <- ext(
  xmin(study_extent_congo) - 1.0,
  xmax(study_extent_congo) + 1.0,
  ymin(study_extent_congo) - 1.0,
  ymax(study_extent_congo) + 1.0
)

study_extent_amazon <- ext(mcps_amazon)
study_extent_amazon <- ext(
  xmin(study_extent_amazon) - 1.0,
  xmax(study_extent_amazon) + 1.0,
  ymin(study_extent_amazon) - 1.0,
  ymax(study_extent_amazon) + 1.0
)

r_congo_cropped <- crop(r_congo, study_extent_congo)
r_amazon_cropped <- crop(r_amazon, study_extent_amazon)

# --- 3. Run Predictions on Rasters --------------------------------------------
cat("Predicting biomass and uncertainty on Congo raster...\n")
congo_cells <- as.data.frame(r_congo_cropped[["uoi"]], cells = TRUE, xy = TRUE, na.rm = TRUE)
names(congo_cells)[names(congo_cells) == "uoi"] <- "uoi"

# Predict link scale and standard error
pred_congo <- predict(m_best, newdata = congo_cells, type = "link", se.fit = TRUE)

congo_cells$pred <- exp(pred_congo$fit)
congo_cells$lower <- exp(pred_congo$fit - 1.96 * pred_congo$se.fit)
congo_cells$upper <- exp(pred_congo$fit + 1.96 * pred_congo$se.fit)
congo_cells$ci_range <- congo_cells$upper - congo_cells$lower

# Create blank SpatRasters and assign values
r_pred_congo <- rast(r_congo_cropped[["uoi"]])
names(r_pred_congo) <- "pred"
values(r_pred_congo) <- NA
r_pred_congo[congo_cells$cell] <- as.vector(congo_cells$pred)

r_unc_congo <- rast(r_congo_cropped[["uoi"]])
names(r_unc_congo) <- "ci_range"
values(r_unc_congo) <- NA
r_unc_congo[congo_cells$cell] <- as.vector(congo_cells$ci_range)

cat("Predicting biomass and uncertainty on Amazon raster...\n")
amazon_cells <- as.data.frame(r_amazon_cropped[["uoi"]], cells = TRUE, xy = TRUE, na.rm = TRUE)
names(amazon_cells)[names(amazon_cells) == "uoi"] <- "uoi"

# Predict link scale and standard error
pred_amazon <- predict(m_best, newdata = amazon_cells, type = "link", se.fit = TRUE)

amazon_cells$pred <- exp(pred_amazon$fit)
amazon_cells$lower <- exp(pred_amazon$fit - 1.96 * pred_amazon$se.fit)
amazon_cells$upper <- exp(pred_amazon$fit + 1.96 * pred_amazon$se.fit)
amazon_cells$ci_range <- amazon_cells$upper - amazon_cells$lower

# Create blank SpatRasters and assign values
r_pred_amazon <- rast(r_amazon_cropped[["uoi"]])
names(r_pred_amazon) <- "pred"
values(r_pred_amazon) <- NA
r_pred_amazon[amazon_cells$cell] <- as.vector(amazon_cells$pred)

r_unc_amazon <- rast(r_amazon_cropped[["uoi"]])
names(r_unc_amazon) <- "ci_range"
values(r_unc_amazon) <- NA
r_unc_amazon[amazon_cells$cell] <- as.vector(amazon_cells$ci_range)

# --- 4. Build ggplot Panels ---------------------------------------------------
cat("Building 4-panel predictive maps...\n")

# Retrieve font
base_family <- if (requireNamespace("showtext", quietly = TRUE)) "Helvetica Neue" else "sans"

# Panel A: Congo Prediction
p_pred_congo <- ggplot() +
  geom_spatraster(data = r_pred_congo) +
  scale_fill_viridis_c(
    option = "inferno",
    name = "Predicted\nBiomass\nIndex",
    limits = c(0, 5000),
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
    title = "A. Congo Basin Predicted Mammal Biomass Index",
    subtitle = "Landscape-scale prediction from weight-normalized Tweedie GLM"
  )

# Panel B: Congo Uncertainty (95% CI Range)
p_unc_congo <- ggplot() +
  geom_spatraster(data = r_unc_congo) +
  scale_fill_viridis_c(
    option = "mako",
    name = "95% CI\nRange\n(Uncertainty)",
    limits = c(0, 5000),
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
    subtitle = "Spatial standard error width (Upper 95% CI - Lower 95% CI)"
  )

# Panel C: Amazon Prediction
p_pred_amazon <- ggplot() +
  geom_spatraster(data = r_pred_amazon) +
  scale_fill_viridis_c(
    option = "inferno",
    name = "Predicted\nBiomass\nIndex",
    limits = c(0, 5000),
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
    title = "C. Amazon Basin Predicted Mammal Biomass Index",
    subtitle = "Landscape-scale prediction from weight-normalized Tweedie GLM"
  )

# Panel D: Amazon Uncertainty (95% CI Range)
p_unc_amazon <- ggplot() +
  geom_spatraster(data = r_unc_amazon) +
  scale_fill_viridis_c(
    option = "mako",
    name = "95% CI\nRange\n(Uncertainty)",
    limits = c(0, 5000),
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
    subtitle = "Spatial standard error width (Upper 95% CI - Lower 95% CI)"
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
  filename = "outputs/congo_ct_gee_biomass_predictions_and_uncertainty.png",
  type = "double",
  height_cm = 16.0
)

cat("✓ Saved: outputs/congo_ct_gee_biomass_predictions_and_uncertainty.png\n")
cat("=== Spatial Prediction Mapping Complete ===\n")
