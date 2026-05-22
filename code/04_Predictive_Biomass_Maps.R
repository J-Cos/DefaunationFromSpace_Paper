# =============================================================================
# code/04_Predictive_Biomass_Maps.R
#
# Consolidates biomass predictive mapping for both the 5 km core scale and
# the 20 km peak predictive scale. Implements functional programming and the
# DRY (Don't Repeat Yourself) principle.
#
# Inputs:
#   - outputs/EOdata/analysis_stack_5000_{Basin}.tif
#   - outputs/EOdata/analysis_stack_20000_{Basin}.tif
#   - outputs/framework2_best_model.RDS
#
# Outputs:
#   - outputs/framework2_biomass_predictions_5km.png (and .pdf)
#   - outputs/best_model_biomass_predictions_20km.png (and .pdf)
# =============================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(cowplot)
library(dplyr)
library(readr)
library(mgcv)

cat("=== Generating Standing Mammal Biomass Predictive Maps (5 km & 20 km) ===\n\n")

# --- 1. Load Theme, Data & Vectors --------------------------------------------
source("code/functions/theme_pnas.R")

dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

geojson_path <- "outputs/camera_traps_robust_buffered_mcps.geojson"
if (!file.exists(geojson_path)) {
  stop("MCP GeoJSON missing: ", geojson_path)
}
mcps <- terra::vect(geojson_path)
mcps_congo <- mcps[mcps$region == "Congo", ]
mcps_amazon <- mcps[mcps$region == "Amazon", ]

# Load country outlines
countries_v <- terra::vect("data/world-administrative-boundaries")

# Load Saved Best Framework 2 Model
model_path <- "outputs/framework2_best_model.RDS"
if (!file.exists(model_path)) {
  stop("Framework 2 best model RDS file missing: ", model_path)
}
m_best <- readRDS(model_path)
s_best <- summary(m_best)

cat("★ Loaded Best Framework 2 Model from RDS successfully.\n")
cat(sprintf("Deviance Explained = %.2f%%\n\n", s_best$dev.expl * 100))


# --- 2. Functional Mapping Routine (DRY Implementation) ---------------------
generate_predictive_map <- function(scale_m, output_base_name, map_title_suffix) {
  cat(sprintf("--- Generating map at scale: %d m ---\n", scale_m))
  
  # Load specific scale TIFF stacks
  r_congo_path <- sprintf("outputs/EOdata/analysis_stack_%d_Congo.tif", scale_m)
  r_amazon_path <- sprintf("outputs/EOdata/analysis_stack_%d_Amazon.tif", scale_m)
  
  # Fallback to synthetic if needed
  if (!file.exists(r_congo_path)) r_congo_path <- sprintf("outputs/synthetic_EOdata/analysis_stack_%d_Congo.tif", scale_m)
  if (!file.exists(r_amazon_path)) r_amazon_path <- sprintf("outputs/synthetic_EOdata/analysis_stack_%d_Amazon.tif", scale_m)
  
  if (!file.exists(r_congo_path) || !file.exists(r_amazon_path)) {
    stop(sprintf("GeoTIFF stacks for scale %d m are missing.", scale_m))
  }
  
  r_congo <- rast(r_congo_path)
  r_amazon <- rast(r_amazon_path)
  
  aggregate_names <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                       "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  names(r_congo) <- aggregate_names
  names(r_amazon) <- aggregate_names
  
  # Crop using camera trap buffered MCP extents with 1.0 degree padding
  study_extent_congo <- ext(mcps_congo)
  study_extent_congo <- ext(
    xmin(study_extent_congo) - 1.0, xmax(study_extent_congo) + 1.0,
    ymin(study_extent_congo) - 1.0, ymax(study_extent_congo) + 1.0
  )
  
  study_extent_amazon <- ext(mcps_amazon)
  study_extent_amazon <- ext(
    xmin(study_extent_amazon) - 1.0, xmax(study_extent_amazon) + 1.0,
    ymin(study_extent_amazon) - 1.0, ymax(study_extent_amazon) + 1.0
  )
  
  r_congo_cropped <- crop(r_congo, study_extent_congo)
  r_amazon_cropped <- crop(r_amazon, study_extent_amazon)
  
  # Project and crop country borders for neatness
  countries_congo <- crop(project(countries_v, crs(r_congo_cropped)), study_extent_congo)
  countries_amazon <- crop(project(countries_v, crs(r_amazon_cropped)), study_extent_amazon)
  
  # Run Predictions on Rasters
  covariate_bands <- c("uoi", "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  
  # Congo predictions
  congo_cells <- as.data.frame(r_congo_cropped[[covariate_bands]], cells = TRUE, xy = TRUE, na.rm = TRUE)
  congo_cells$basin <- factor("Congo", levels = c("Amazon", "Congo"))
  congo_cells$pred <- exp(predict(m_best, newdata = congo_cells, type = "link"))
  
  r_pred_congo <- rast(r_congo_cropped[["uoi"]])
  names(r_pred_congo) <- "pred"
  values(r_pred_congo) <- NA
  r_pred_congo[congo_cells$cell] <- as.vector(congo_cells$pred)
  
  # Amazon predictions
  amazon_cells <- as.data.frame(r_amazon_cropped[[covariate_bands]], cells = TRUE, xy = TRUE, na.rm = TRUE)
  amazon_cells$basin <- factor("Amazon", levels = c("Amazon", "Congo"))
  amazon_cells$pred <- exp(predict(m_best, newdata = amazon_cells, type = "link"))
  
  r_pred_amazon <- rast(r_amazon_cropped[["uoi"]])
  names(r_pred_amazon) <- "pred"
  values(r_pred_amazon) <- NA
  r_pred_amazon[amazon_cells$cell] <- as.vector(amazon_cells$pred)
  
  # Build ggplot Panels
  t_theme <- theme_pnas(base_size = 8)
  
  # Panel A: Congo Basin Map (Stretched 0-3000)
  p_congo <- ggplot() +
    geom_spatraster(data = r_pred_congo, aes(fill = pred)) +
    scale_fill_viridis_c(
      option = "plasma",
      name = "Standing Mammal Biomass Index (Congo)",
      limits = c(0, 3000),
      oob = scales::squish,
      na.value = "transparent"
    ) +
    geom_spatvector(data = countries_congo, fill = NA, colour = "grey80", linewidth = 0.25) +
    geom_spatvector(data = mcps_congo, fill = NA, color = "black", linewidth = 0.4, linetype = "dashed") +
    t_theme +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 7, face = "bold"),
      legend.text = element_text(size = 6),
      legend.key.height = unit(0.2, "cm"),
      legend.key.width = unit(1.0, "cm"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(
      title = sprintf("Congo Basin (%s)", map_title_suffix),
      subtitle = "Predicted standing mammal biomass"
    )
  
  # Panel B: Amazon Basin Map (Stretched 0-800)
  p_amazon <- ggplot() +
    geom_spatraster(data = r_pred_amazon, aes(fill = pred)) +
    scale_fill_viridis_c(
      option = "plasma",
      name = "Standing Mammal Biomass Index (Amazon)",
      limits = c(0, 800),
      oob = scales::squish,
      na.value = "transparent"
    ) +
    geom_spatvector(data = countries_amazon, fill = NA, colour = "grey80", linewidth = 0.25) +
    geom_spatvector(data = mcps_amazon, fill = NA, color = "black", linewidth = 0.4, linetype = "dashed") +
    t_theme +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 7, face = "bold"),
      legend.text = element_text(size = 6),
      legend.key.height = unit(0.2, "cm"),
      legend.key.width = unit(1.0, "cm"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(
      title = sprintf("Amazon Basin (%s)", map_title_suffix),
      subtitle = "Predicted standing mammal biomass"
    )
  
  # Stack vertically to preserve true relative aspect ratios (ncol = 1)
  fig_combined <- cowplot::plot_grid(
    p_congo, p_amazon,
    ncol = 1,
    rel_heights = c(0.48, 1.0),
    align = "v",
    axis = "lr"
  )
  
  # Save PNG and PDF (PNAS single-column width = 8.7 cm)
  png_file <- sprintf("outputs/%s.png", output_base_name)
  pdf_file <- sprintf("outputs/%s.pdf", output_base_name)
  
  ggsave(filename = png_file, plot = fig_combined, width = 8.7, height = 16.5, units = "cm", dpi = 600, bg = "white")
  ggsave(filename = pdf_file, plot = fig_combined, width = 8.7, height = 16.5, units = "cm", dpi = 600, bg = "white")
  
  # Copy PNG to active brain artifacts folder
  brain_artifacts_dir <- "/home/j/.gemini/antigravity/brain/913e5cea-7c99-4b21-8124-ea8455da8457"
  if (file.exists(brain_artifacts_dir)) {
    file.copy(png_file, file.path(brain_artifacts_dir, sprintf("%s.png", output_base_name)), overwrite = TRUE)
    cat(sprintf("✓ Copied %s.png to brain artifacts folder.\n", output_base_name))
  }
  
  cat(sprintf("✓ Successfully saved map figure to: %s\n\n", png_file))
}


# --- 3. Run Map Projections for 5 km and 20 km Scales ------------------------
generate_predictive_map(5000, "framework2_biomass_predictions_5km", "5 km")
generate_predictive_map(20000, "best_model_biomass_predictions_20km", "20 km")

cat("=== 04: Standing Mammal Biomass Predictive Maps Generation Complete ===\n")
