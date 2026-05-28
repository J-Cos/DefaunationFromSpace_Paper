# =============================================================================
# code/05_Predictive_Biomass_Maps.R
#
# Generates standing mammal biomass predictive mapping for both the 5 km core scale
# and the 20 km peak predictive scale, using the LOBO-selected best model
# (M2.1: UOI Only). Standardizes predictions to the Z-score normalized log-scale MAE
# and renders them over a solid soft-grey country land background polygon layers.
#
# Inputs:
#   - outputs/EOdata/analysis_stack_5000_{Basin}.tif
#   - outputs/EOdata/analysis_stack_20000_{Basin}.tif
#   - outputs/framework2_best_model.RDS
#   - outputs/framework2_covariate_model_selection.csv
#
# Outputs:
#   - figures/figureS5.png (and .pdf) (5 km predictive maps)
#   - figures/figure5.png (and .pdf) (20 km predictive maps)
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
#' Projects the best LOBO-selected Tweedie GLM from Framework 2 across the tropical
#' landscapes at designated spatial resolutions, standardizes values in OOS MAE units,
#' and saves prediction maps to the figures folder.
#'
#' @param scales Numeric vector. Spatial resolutions in meters (default: c(5000, 20000))
#' @param outputs_dir Character. Directory to load models and save RDS (default: "outputs")
#' @param figures_dir Character. Directory to save figures (default: "figures")
#'
#' @return Invisible list of lists of SpatRasters (congo and amazon) of predicted biomass for each scale.
#' @export
run_predictive_biomass_mapping <- function(scales = c(5000, 20000), outputs_dir = "outputs", figures_dir = "figures") {
  cat("=== Generating Z-Score Standardized Standing Mammal Biomass Predictive Maps (5 km & 20 km) ===\n\n")
  
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
  cat(sprintf("Formula: %s\n", paste(deparse(formula(m_best)), collapse = " ")))
  
  # Load Model Selection CSV to extract OOS MAE (log)
  sel_path <- file.path(outputs_dir, "framework2_covariate_model_selection.csv")
  if (!file.exists(sel_path)) {
    stop("Model selection CSV missing: ", sel_path)
  }
  lobo_sel <- read_csv(sel_path, show_col_types = FALSE)
  best_row <- lobo_sel %>% filter(Model == "M2.1: UOI Only")
  if (nrow(best_row) == 0) {
    best_row <- lobo_sel[1, ]
  }
  OOS_MAE_log <- best_row$OOS_MAE_log[1]
  cat(sprintf("OOS MAE (log1p scale) = %.4f\n", OOS_MAE_log))
  
  # Extract calibrated scale data to get full-sample mean of log1p(biomass)
  source("code/functions/calibration_helpers.R")
  joined_data <- extract_scale_data(5000)
  covs_to_check <- c("uoi", "forest_fraction", "B_H_index")
  joined_data <- joined_data %>% filter(complete.cases(joined_data[, covs_to_check]))
  
  mean_log_y_obs <- mean(log1p(joined_data$B_H_index))
  cat(sprintf("Full-sample mean log1p(biomass) = %.4f\n\n", mean_log_y_obs))
  
  # Map lookup table for filenames
  file_suffix_map <- list(
    "5000" = list(name = "figureS5", title = "5 km"),
    "20000" = list(name = "figure5", title = "20 km")
  )
  
  output_rasts <- list()
  
  # --- 2. Functional Mapping Routine -------------------------------------------
  generate_predictive_map <- function(scale_m, output_base_name, map_title_suffix) {
    cat(sprintf("--- Generating map at scale: %d m ---\n", scale_m))
    
    # Load specific scale TIFF stacks
    r_congo_path <- file.path(outputs_dir, "EOdata", sprintf("analysis_stack_%d_Congo.tif", scale_m))
    r_amazon_path <- file.path(outputs_dir, "EOdata", sprintf("analysis_stack_%d_Amazon.tif", scale_m))
    
    # Dynamic aggregation fallback from 5000m real stack if target scale real file is missing
    is_real_file_missing <- !file.exists(file.path(outputs_dir, "EOdata", sprintf("analysis_stack_%d_Congo.tif", scale_m))) ||
                            !file.exists(file.path(outputs_dir, "EOdata", sprintf("analysis_stack_%d_Amazon.tif", scale_m)))
    
    r_congo_5000_path <- file.path(outputs_dir, "EOdata", "analysis_stack_5000_Congo.tif")
    r_amazon_5000_path <- file.path(outputs_dir, "EOdata", "analysis_stack_5000_Amazon.tif")
    
    if (is_real_file_missing && file.exists(r_congo_5000_path) && file.exists(r_amazon_5000_path) && scale_m > 5000) {
      fact <- scale_m / 5000
      cat(sprintf("  ✓ Dynamically aggregating real 5,000m stack by factor of %d to %d m...\n", fact, scale_m))
      r_congo <- terra::aggregate(rast(r_congo_5000_path), fact = fact, fun = "mean", na.rm = TRUE)
      r_amazon <- terra::aggregate(rast(r_amazon_5000_path), fact = fact, fun = "mean", na.rm = TRUE)
      
      # Bypass file loading
      r_congo_path <- "dynamic_aggregated"
      r_amazon_path <- "dynamic_aggregated"
    }
    
    if (r_congo_path != "dynamic_aggregated") {
      # Fallback to synthetic if needed
      if (!file.exists(r_congo_path)) r_congo_path = file.path(outputs_dir, "synthetic_EOdata", sprintf("analysis_stack_%d_Congo.tif", scale_m))
      if (!file.exists(r_amazon_path)) r_amazon_path = file.path(outputs_dir, "synthetic_EOdata", sprintf("analysis_stack_%d_Amazon.tif", scale_m))
      
      if (!file.exists(r_congo_path) || !file.exists(r_amazon_path)) {
        stop(sprintf("GeoTIFF stacks for scale %d m are missing.", scale_m))
      }
      
      r_congo <- rast(r_congo_path)
      r_amazon <- rast(r_amazon_path)
    }
    
    # SE Asia stack loading
    r_seasia_path <- file.path(outputs_dir, "EOdata", sprintf("analysis_stack_%d_SE_Asia.tif", scale_m))
    is_seasia_real_missing <- !file.exists(r_seasia_path)
    r_seasia_5000_path <- file.path(outputs_dir, "EOdata", "analysis_stack_5000_SE_Asia.tif")
    
    if (is_seasia_real_missing && file.exists(r_seasia_5000_path) && scale_m > 5000) {
      fact <- scale_m / 5000
      cat(sprintf("  ✓ Dynamically aggregating real 5,000m SE Asia stack by factor of %d to %d m...\n", fact, scale_m))
      r_seasia <- terra::aggregate(rast(r_seasia_5000_path), fact = fact, fun = "mean", na.rm = TRUE)
    } else if (file.exists(r_seasia_path)) {
      r_seasia <- rast(r_seasia_path)
    } else {
      r_seasia_synth <- file.path(outputs_dir, "synthetic_EOdata", sprintf("analysis_stack_%d_SE_Asia.tif", scale_m))
      if (file.exists(r_seasia_synth)) {
        r_seasia <- rast(r_seasia_synth)
      } else {
        r_seasia <- NULL
      }
    }
    
    aggregate_names <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                         "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
    names(r_congo) <- aggregate_names
    names(r_amazon) <- aggregate_names
    if (!is.null(r_seasia)) {
      names(r_seasia) <- aggregate_names
    }
    
    # Equal dimension boxes: Exactly 50° wide × 30° high (-15 to 15 Latitude)
    ext_amazon_map <- ext(-85, -35, -15, 15)
    ext_congo_map  <- ext(-5, 45, -15, 15)
    ext_seasia_map <- ext(90, 140, -15, 15)
    
    r_congo_cropped <- crop(r_congo, ext_congo_map, extend = TRUE)
    r_amazon_cropped <- crop(r_amazon, ext_amazon_map, extend = TRUE)
    if (!is.null(r_seasia)) {
      r_seasia_cropped <- crop(r_seasia, ext_seasia_map, extend = TRUE)
    }
    
    # Project and crop country borders
    countries_congo <- crop(project(countries_v, crs(r_congo_cropped)), ext_congo_map)
    countries_amazon <- crop(project(countries_v, crs(r_amazon_cropped)), ext_amazon_map)
    if (!is.null(r_seasia)) {
      countries_seasia <- crop(project(countries_v, crs(r_seasia_cropped)), ext_seasia_map)
    }
    
    # Run Predictions on Rasters and transform to Z-score normalized OOS MAE units
    covariate_bands <- c("uoi", "forest_fraction")
    
    # Congo predictions
    congo_cells <- as.data.frame(r_congo_cropped[[covariate_bands]], cells = TRUE, xy = TRUE, na.rm = TRUE)
    congo_cells$pred <- predict(m_best, newdata = congo_cells, type = "response")
    congo_cells$zscore_mae <- (log1p(congo_cells$pred) - mean_log_y_obs) / OOS_MAE_log
    
    r_pred_congo <- rast(r_congo_cropped[["uoi"]])
    names(r_pred_congo) <- "zscore_mae"
    values(r_pred_congo) <- NA
    r_pred_congo[congo_cells$cell] <- as.vector(congo_cells$zscore_mae)
    
    # Amazon predictions
    amazon_cells <- as.data.frame(r_amazon_cropped[[covariate_bands]], cells = TRUE, xy = TRUE, na.rm = TRUE)
    amazon_cells$pred <- predict(m_best, newdata = amazon_cells, type = "response")
    amazon_cells$zscore_mae <- (log1p(amazon_cells$pred) - mean_log_y_obs) / OOS_MAE_log
    
    r_pred_amazon <- rast(r_amazon_cropped[["uoi"]])
    names(r_pred_amazon) <- "zscore_mae"
    values(r_pred_amazon) <- NA
    r_pred_amazon[amazon_cells$cell] <- as.vector(amazon_cells$zscore_mae)
    
    # SE Asia predictions
    r_pred_seasia <- NULL
    if (!is.null(r_seasia)) {
      seasia_cells <- as.data.frame(r_seasia_cropped[[covariate_bands]], cells = TRUE, xy = TRUE, na.rm = TRUE)
      if (nrow(seasia_cells) > 0) {
        seasia_cells$pred <- predict(m_best, newdata = seasia_cells, type = "response")
        seasia_cells$zscore_mae <- (log1p(seasia_cells$pred) - mean_log_y_obs) / OOS_MAE_log
        
        r_pred_seasia <- rast(r_seasia_cropped[["uoi"]])
        names(r_pred_seasia) <- "zscore_mae"
        values(r_pred_seasia) <- NA
        r_pred_seasia[seasia_cells$cell] <- as.vector(seasia_cells$zscore_mae)
      }
    }
    
    # Load and crop protected areas
    pa_congo_file <- file.path(outputs_dir, "WDPA_congo_500km2.gpkg")
    pa_amazon_file <- file.path(outputs_dir, "WDPA_amazon_500km2.gpkg")
    
    pas_congo_cropped <- NULL
    if (file.exists(pa_congo_file)) {
      pas_congo <- vect(pa_congo_file)
      pas_congo <- pas_congo[pas_congo$REP_AREA >= 1000 | pas_congo$GIS_AREA >= 1000, ]
      if (nrow(pas_congo) > 0) {
        pas_congo_proj <- project(pas_congo, crs(r_congo_cropped))
        pas_congo_cropped <- crop(pas_congo_proj, ext_congo_map)
      }
    }
    
    pas_amazon_cropped <- NULL
    if (file.exists(pa_amazon_file)) {
      pas_amazon <- vect(pa_amazon_file)
      pas_amazon <- pas_amazon[pas_amazon$REP_AREA >= 1000 | pas_amazon$GIS_AREA >= 1000, ]
      if (nrow(pas_amazon) > 0) {
        pas_amazon_proj <- project(pas_amazon, crs(r_amazon_cropped))
        pas_amazon_cropped <- crop(pas_amazon_proj, ext_amazon_map)
      }
    }
    
    # Build ggplot Panels
    t_theme <- theme_pnas(base_size = 8)
    
    # Diverging Red-Grey-Blue color scale centered at 0.0 representing mean biomass
    fill_scale <- scale_fill_gradientn(
      colors = c("#b2182b", "#ef8a62", "#f7f7f7", "#67a9cf", "#2166ac"),
      name = "Predicted Biomass Deviation from Global Mean (in units of OOS log-scale MAE)",
      limits = c(-2.5, 2.5),
      breaks = c(-2.0, -1.0, 0, 1.0, 2.0),
      labels = c("-2.0 MAE\n(Low Biomass)", "-1.0 MAE", "0.0\n(Mean Biomass)", "+1.0 MAE", "+2.0 MAE\n(High Biomass)"),
      oob = scales::squish,
      na.value = "transparent",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        barwidth = unit(8.5, "cm"),
        barheight = unit(0.24, "cm")
      )
    )
    
    # Panel A: Amazon Map (with solid soft-grey land borders drawn first)
    p_amazon <- ggplot() +
      geom_spatvector(data = countries_amazon, fill = "#F2F4F4", color = "grey80", linewidth = 0.25) +
      geom_spatraster(data = r_pred_amazon, aes(fill = zscore_mae)) +
      fill_scale
    
    if (!is.null(pas_amazon_cropped) && nrow(pas_amazon_cropped) > 0) {
      p_amazon <- p_amazon + geom_spatvector(data = pas_amazon_cropped, fill = NA, color = "black", linewidth = 0.12)
    }
    
    p_amazon <- p_amazon +
      geom_spatvector(data = mcps_amazon, fill = NA, color = "black", linewidth = 0.4, linetype = "dashed") +
      coord_sf(xlim = c(-85, -35), ylim = c(-15, 15), expand = FALSE) +
      t_theme +
      theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
      ) +
      labs(
        title = sprintf("A. Neotropical Basin (Amazon, %s)", map_title_suffix)
      )
    
    # Panel B: Congo Map (with solid soft-grey land borders drawn first)
    p_congo <- ggplot() +
      geom_spatvector(data = countries_congo, fill = "#F2F4F4", color = "grey80", linewidth = 0.25) +
      geom_spatraster(data = r_pred_congo, aes(fill = zscore_mae)) +
      fill_scale
    
    if (!is.null(pas_congo_cropped) && nrow(pas_congo_cropped) > 0) {
      p_congo <- p_congo + geom_spatvector(data = pas_congo_cropped, fill = NA, color = "black", linewidth = 0.12)
    }
    
    p_congo <- p_congo +
      geom_spatvector(data = mcps_congo, fill = NA, color = "black", linewidth = 0.4, linetype = "dashed") +
      coord_sf(xlim = c(-5, 45), ylim = c(-15, 15), expand = FALSE) +
      t_theme +
      theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
      ) +
      labs(
        title = sprintf("B. Afrotropical Basin (Congo, %s)", map_title_suffix)
      )
      
    # Crop blank SE Asia protected areas
    pa_seasia_file <- file.path(outputs_dir, "WDPA_seasia_500km2.gpkg")
    pas_seasia_cropped <- NULL
    if (file.exists(pa_seasia_file)) {
      pas_seasia <- vect(pa_seasia_file)
      pas_seasia <- pas_seasia[pas_seasia$REP_AREA >= 1000 | pas_seasia$GIS_AREA >= 1000, ]
      if (nrow(pas_seasia) > 0) {
        pas_seasia_proj <- project(pas_seasia, crs(r_seasia_cropped))
        pas_seasia_cropped <- crop(pas_seasia_proj, ext_seasia_map)
      }
    }
    
    # Panel C: SE Asia Map (with solid soft-grey land borders drawn first)
    if (!is.null(r_pred_seasia)) {
      p_seasia <- ggplot() +
        geom_spatvector(data = countries_seasia, fill = "#F2F4F4", color = "grey80", linewidth = 0.25) +
        geom_spatraster(data = r_pred_seasia, aes(fill = zscore_mae)) +
        fill_scale
      
      if (!is.null(pas_seasia_cropped) && nrow(pas_seasia_cropped) > 0) {
        p_seasia <- p_seasia + geom_spatvector(data = pas_seasia_cropped, fill = NA, color = "black", linewidth = 0.12)
      }
      
      mcps_seasia <- mcps[mcps$region == "SE_Asia", ]
      if (nrow(mcps_seasia) > 0) {
        p_seasia <- p_seasia + geom_spatvector(data = mcps_seasia, fill = NA, color = "black", linewidth = 0.4, linetype = "dashed")
      }
      
      p_seasia <- p_seasia +
        coord_sf(xlim = c(90, 140), ylim = c(-15, 15), expand = FALSE) +
        t_theme +
        theme(
          legend.position = "bottom",
          legend.title = element_text(size = 6.5, face = "bold"),
          legend.text = element_text(size = 5.5),
          legend.margin = margin(t = 2, b = 2, unit = "pt"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank()
        ) +
        labs(
          title = sprintf("C. Indo-Malayan Basin (Southeast Asia, %s)", map_title_suffix)
        )
    } else {
      countries_seasia <- crop(project(countries_v, crs(r_congo_cropped)), ext_seasia_map)
      p_seasia <- ggplot() +
        geom_spatvector(data = countries_seasia, fill = "#F2F4F4", color = "grey80", linewidth = 0.25) +
        coord_sf(xlim = c(90, 140), ylim = c(-15, 15), expand = FALSE) +
        t_theme +
        theme(
          legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "#FFFFFF", color = NA)
        ) +
        annotate("text", x = 115, y = 0, label = "Southeast Asia: Stack Data Missing", fontface = "italic", size = 2.4, color = "grey40") +
        labs(
          title = sprintf("C. Indo-Malayan Basin (Southeast Asia, %s)", map_title_suffix)
        )
    }
    
    # Combine maps vertically
    fig_combined <- cowplot::plot_grid(
      p_amazon, p_congo, p_seasia,
      ncol = 1,
      align = "v",
      axis = "lr",
      rel_heights = c(1.0, 1.0, 1.15)
    )
    
    # Save PNG and PDF
    png_file <- file.path(figures_dir, sprintf("%s.png", output_base_name))
    pdf_file <- file.path(figures_dir, sprintf("%s.pdf", output_base_name))
    
    ggsave(filename = png_file, plot = fig_combined, width = 12.0, height = 24.0, units = "cm", dpi = 600, bg = "white")
    ggsave(filename = pdf_file, plot = fig_combined, width = 12.0, height = 24.0, units = "cm", dpi = 600, bg = "white")
    
    # Mirror to active brain artifacts folder
    brain_artifacts_dir <- "/home/j/.gemini/antigravity/brain/8f51df52-4604-48e0-9ce8-1c52d1cb241c"
    if (file.exists(brain_artifacts_dir)) {
      file.copy(png_file, file.path(brain_artifacts_dir, sprintf("%s.png", output_base_name)), overwrite = TRUE)
      file.copy(pdf_file, file.path(brain_artifacts_dir, sprintf("%s.pdf", output_base_name)), overwrite = TRUE)
      cat(sprintf("✓ Copied %s.png and .pdf to brain artifacts folder.\n", output_base_name))
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
  
  cat("=== Standing Mammal Biomass Predictive Maps Generation Complete ===\n")
  return(invisible(output_rasts))
}

# --- Execute directly if called from terminal ---
run_predictive_biomass_mapping()
