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
#   - figures/framework2_biomass_predictions_5km.png (and .pdf)
#   - figures/best_model_biomass_predictions_20km.png (and .pdf)
# =============================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(cowplot)
library(dplyr)
library(readr)
library(mgcv)

#' Run Predictive Biomass Mapping
#'
#' Projects the best-fitting Tweedie GLM from Framework 2 across the tropical
#' landscape at designated spatial resolutions, and saves prediction maps to the
#' figures folder.
#'
#' @param scales Numeric vector. Spatial resolutions in meters (default: c(5000, 20000))
#' @param outputs_dir Character. Directory to load models and save RDS (default: "outputs")
#' @param figures_dir Character. Directory to save figures (default: "figures")
#'
#' @return Invisible list of lists of SpatRasters (congo and amazon) of predicted biomass for each scale.
#' @export
run_predictive_biomass_mapping <- function(scales = c(5000, 20000), outputs_dir = "outputs", figures_dir = "figures") {
  cat("=== Generating Standing Mammal Biomass Predictive Maps (5 km & 20 km) ===\n\n")
  
  # --- 1. Load Theme, Data & Vectors --------------------------------------------
  source("code/functions/theme_pnas.R")
  
  dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  geojson_path <- file.path(outputs_dir, "camera_traps_robust_buffered_mcps.geojson")
  if (!file.exists(geojson_path)) {
    stop("MCP GeoJSON missing: ", geojson_path)
  }
  mcps <- terra::vect(geojson_path)
  mcps_congo <- mcps[mcps$region == "Congo", ]
  mcps_amazon <- mcps[mcps$region == "Amazon", ]
  
  # Load country outlines
  countries_v <- terra::vect("data/world-administrative-boundaries")
  
  # Load Saved Best Framework 2 Model
  model_path <- file.path(outputs_dir, "framework2_best_model.RDS")
  if (!file.exists(model_path)) {
    stop("Framework 2 best model RDS file missing: ", model_path)
  }
  m_best <- readRDS(model_path)
  s_best <- summary(m_best)
  
  cat("★ Loaded Best Framework 2 Model from RDS successfully.\n")
  cat(sprintf("Deviance Explained = %.2f%%\n\n", s_best$dev.expl * 100))
  
  # Map lookup table for filenames
  file_suffix_map <- list(
    "5000" = list(name = "framework2_biomass_predictions_5km", title = "5 km"),
    "20000" = list(name = "best_model_biomass_predictions_20km", title = "20 km")
  )
  
  output_rasts <- list()
  
  # --- 2. Functional Mapping Routine (DRY Implementation) ---------------------
  generate_predictive_map <- function(scale_m, output_base_name, map_title_suffix) {
    cat(sprintf("--- Generating map at scale: %d m ---\n", scale_m))
    
    # Load specific scale TIFF stacks
    r_congo_path <- file.path(outputs_dir, "EOdata", sprintf("analysis_stack_%d_Congo.tif", scale_m))
    r_amazon_path <- file.path(outputs_dir, "EOdata", sprintf("analysis_stack_%d_Amazon.tif", scale_m))
    
    # Fallback to synthetic if needed
    if (!file.exists(r_congo_path)) r_congo_path = file.path(outputs_dir, "synthetic_EOdata", sprintf("analysis_stack_%d_Congo.tif", scale_m))
    if (!file.exists(r_amazon_path)) r_amazon_path = file.path(outputs_dir, "synthetic_EOdata", sprintf("analysis_stack_%d_Amazon.tif", scale_m))
    
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
    
    # Save PNG and PDF to the designated figures directory
    png_file <- file.path(figures_dir, sprintf("%s.png", output_base_name))
    pdf_file <- file.path(figures_dir, sprintf("%s.pdf", output_base_name))
    
    ggsave(filename = png_file, plot = fig_combined, width = 8.7, height = 16.5, units = "cm", dpi = 600, bg = "white")
    ggsave(filename = pdf_file, plot = fig_combined, width = 8.7, height = 16.5, units = "cm", dpi = 600, bg = "white")
    
    # Copy PNG to active brain artifacts folder
    brain_artifacts_dir <- "/home/j/.gemini/antigravity/brain/913e5cea-7c99-4b21-8124-ea8455da8457"
    if (file.exists(brain_artifacts_dir)) {
      file.copy(png_file, file.path(brain_artifacts_dir, sprintf("%s.png", output_base_name)), overwrite = TRUE)
      cat(sprintf("✓ Copied %s.png to brain artifacts folder.\n", output_base_name))
    }
    
    cat(sprintf("✓ Successfully saved map figure to: %s\n\n", png_file))
    
    return(list(congo = r_pred_congo, amazon = r_pred_amazon))
  }
  
  # --- 3. Run Map Projections for each scale ---
  for (scale_m in scales) {
    scale_str <- as.character(scale_m)
    meta <- file_suffix_map[[scale_str]]
    if (is.null(meta)) {
      meta <- list(name = sprintf("predicted_biomass_%d_m", scale_m), title = sprintf("%d m", scale_m / 1000))
    }
    
    res_rasts <- generate_predictive_map(scale_m, meta$name, meta$title)
    output_rasts[[scale_str]] <- res_rasts
  }
  
  cat("=== 04: Standing Mammal Biomass Predictive Maps Generation Complete ===\n")
  return(invisible(output_rasts))
}

# --- Execute directly if called from terminal ---
run_predictive_biomass_mapping()
