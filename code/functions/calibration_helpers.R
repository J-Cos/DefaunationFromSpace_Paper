# =============================================================================
# calibration_helpers.R
#
# Shared helper functions to extract multi-scale pixel data within camera trap
# MCP polygons, clean/filter pixels, calculate weights, compute temporal/spatial
# precision and homogeneity, and fit the Framework 1 (Beta Regression) and
# Framework 2 (Tweedie GLM) models using dynamically exported formulas.
#
# Ensures 100% methodology alignment across downstream pipeline scripts:
#   - code/framework1_integrated_analysis.R
#   - code/framework2_integrated_analysis.R
#   - code/08f_Multiscale_Sensitivity_Analysis.R
#   - code/08g_Representative_Scale_Maps.R
#   - code/08h_25km_Multipanel_Synthesis.R
#   - code/13_Best_Model_Predictive_Map.R
# =============================================================================

library(terra)
library(dplyr)
library(mgcv)
library(readr)
library(jsonlite)

source("code/functions/model_convergence.R")
config_path <- if (file.exists("code/config.json")) "code/config.json" else "config.json"
config <- jsonlite::read_json(config_path)
MIN_TRAP_DAYS <- config$clustering$min_trap_days
CLUSTER_THRESHOLD_KM <- as.numeric(config$clustering$threshold_km)



# --- 1. Extract Raw Pixel-Level Data Within MCP Polygons ---------------------
#' Extract Pixel-Level Data Within MCP Polygons
#'
#' @param scale_m Numeric. Spatial resolution in meters (e.g. 5000, 10000, 25000, etc.)
#' @param mcps SpatVector. Optional loaded MCP vector.
#'
#' @return A data frame of raw extracted pixel values with basin tags.
#' @export
extract_scale_pixels <- function(scale_m, mcps = NULL) {
  if (is.null(mcps)) {
    geojson_path <- "outputs/camera_traps_robust_buffered_mcps.geojson"
    if (!file.exists(geojson_path)) {
      stop("Camera trap buffered MCPs missing. Please run code/visualise_camera_traps.py first.")
    }
    mcps <- terra::vect(geojson_path)
  }
  
  # Compute elephant presence dynamically based on spatial ranges
  if (!file.exists("outputs/elephant_ranges.gpkg")) {
    stop("Critical Error: outputs/elephant_ranges.gpkg is missing! Please run 01_FigureS1_Regional_Bounding_Boxes.R first.")
  }
  
  ele_ranges <- terra::vect("outputs/elephant_ranges.gpkg")
  mcps$elephant_present_strict <- 0
  
  # Strict: only extant (status == "Extant", which corresponds to presence == 1)
  ele_ranges_strict <- ele_ranges[ele_ranges$status == "Extant", ]
  if (nrow(ele_ranges_strict) > 0) {
    intersects_strict <- terra::is.related(mcps, ele_ranges_strict, "intersects")
    mcps$elephant_present_strict <- as.numeric(rowSums(as.matrix(intersects_strict)) > 0)
  }
  
  # Possible: all statuses (Extant, Possibly Extant, Possibly Extinct)
  intersects_possible <- terra::is.related(mcps, ele_ranges, "intersects")
  mcps$elephant_present_possible <- as.numeric(rowSums(as.matrix(intersects_possible)) > 0)
  
  # Retain legacy alias for safety
  mcps$elephant_present <- mcps$elephant_present_possible

  
  load_and_aggregate_if_needed <- function(path, scale_val, basin) {
    if (file.exists(path)) {
      return(terra::rast(path))
    }
    # Dynamic aggregation fallback from 5,000m real stack if target scale real file is missing
    path_5000 <- sprintf("outputs/EOdata/analysis_stack_5000_%s.tif", basin)
    if (file.exists(path_5000) && !is.character(scale_val) && scale_val > 5000) {
      fact <- scale_val / 5000
      message(sprintf("✓ Dynamically aggregating real 5,000m %s stack by factor of %d to %d m", basin, fact, scale_val))
      r_5000 <- terra::rast(path_5000)
      r <- terra::aggregate(r_5000, fact = fact, fun = "mean", na.rm = TRUE)
      return(r)
    }
    stop(sprintf("GeoTIFF analysis stack for %s at scale %s is missing.", basin, as.character(scale_val)))
  }

  if (is.character(scale_m) && scale_m == "native") {
    aggregate_names <- c("uoi", "uoi_sd", "rh98", "gedi_n", "elevation", "slope", "hnd", "precip", "clay", "forest_fraction", "Npp_median")
  } else {
    aggregate_names <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                         "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  }

  basins <- c("Congo", "Amazon", "SE_Asia")
  pixel_list <- list()

  for (b in basins) {
    mcps_b <- mcps[mcps$region == b, ]
    if (nrow(mcps_b) == 0) next

    r_path <- if (is.character(scale_m)) {
      sprintf("outputs/EOdata/analysis_stack_%s_%s.tif", scale_m, b)
    } else {
      sprintf("outputs/EOdata/analysis_stack_%d_%s.tif", scale_m, b)
    }

    r_b <- load_and_aggregate_if_needed(r_path, scale_m, b)
    names(r_b) <- aggregate_names

    ext_b <- terra::extract(r_b, mcps_b, df = TRUE, touches = FALSE)
    all_ids <- 1:nrow(mcps_b)
    extracted_ids <- unique(ext_b$ID)
    empty_ids <- setdiff(all_ids, extracted_ids)

    if (length(empty_ids) > 0) {
      ext_touch <- terra::extract(r_b, mcps_b[empty_ids, ], df = TRUE, touches = TRUE)
      ext_touch$ID <- empty_ids[ext_touch$ID]
      ext_b <- rbind(ext_b %>% filter(ID %in% extracted_ids), ext_touch)
    }

    mcp_df <- as.data.frame(mcps_b)
    mcp_df$ID <- 1:nrow(mcp_df)

    pixel_b <- merge(ext_b, mcp_df, by = "ID") %>%
      filter(!is.na(uoi)) %>%
      select(-ID) %>%
      mutate(basin = b)

    if (is.character(scale_m) && scale_m == "native") {
      pixel_b$frip <- NA
    }

    pixel_list[[b]] <- pixel_b
  }

  pixel_data <- do.call(rbind, pixel_list)
  return(pixel_data)
}

# --- 2. Extract, Clean, and Calibrate Scale-Specific Cluster Data ------------
#' Extract, Clean, and Calibrate Scale-Specific Cluster Data
#'
#' @param scale_m Numeric. Spatial resolution in meters (e.g. 5000, 10000, 25000, etc.)
#' @param mcps SpatVector. Optional loaded MCP vector.
#'
#' @return A data frame of aggregated cluster-level variables with weights and homogeneity.
#' @export
extract_scale_data <- function(scale_m, mcps = NULL) {
  pixel_data <- extract_scale_pixels(scale_m, mcps = mcps)
  
  joined_data <- pixel_data %>%
    group_by(cluster_id, region, basin, elephant_present_strict, elephant_present_possible, elephant_present, trap_days, n_species, B_H_index, M_H_index, B_H_gt50, B_H_gt100, B_H_gt1000, megafauna_fraction, p_keep, w_temp_cluster) %>%
    summarise(
      n_pixels = n(),
      uoi_sd = ifelse(is.na(sd(uoi, na.rm = TRUE)), 0, sd(uoi, na.rm = TRUE)),
      uoi = pmax(pmin(mean(uoi, na.rm = TRUE), 1 - 1e-5), 1e-5),
      elevation = mean(elevation, na.rm = TRUE),
      slope = mean(slope, na.rm = TRUE),
      hnd = mean(hnd, na.rm = TRUE),
      precip = mean(precip, na.rm = TRUE),
      clay = mean(clay, na.rm = TRUE),
      forest_fraction = mean(forest_fraction, na.rm = TRUE),
      frip = mean(frip, na.rm = TRUE),
      rh98 = mean(rh98, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(trap_days >= MIN_TRAP_DAYS)
  
  # Compute standard errors of the mean
  joined_data <- joined_data %>%
    mutate(uoi_se = uoi_sd / sqrt(n_pixels))
  
  # Calculate effective trap days discounted by taxonomic resolution issues
  joined_data <- joined_data %>%
    mutate(trap_days_effective = trap_days * p_keep)

  # Compute precision weights with taxonomic effort discount
  reg_uoi <- median(joined_data$uoi_se[joined_data$uoi_se > 0])
  if (is.na(reg_uoi) || reg_uoi == 0) reg_uoi <- 1e-4
  
  joined_data <- joined_data %>%
    mutate(w_uoi = log10(pmax(trap_days_effective, 1.0)) / (uoi_se + reg_uoi))
  
  # Normalize weights so they sum to N (mean = 1) for proper statistical scale
  joined_data$w_uoi_norm <- joined_data$w_uoi / mean(joined_data$w_uoi)
  basin_levels <- c("Amazon", "Congo")
  if ("SE_Asia" %in% joined_data$basin) {
    basin_levels <- c("Amazon", "Congo", "SE_Asia")
  }
  joined_data$basin <- factor(joined_data$basin, levels = basin_levels)
  joined_data$elephant_present_strict <- factor(ifelse(joined_data$elephant_present_strict == 1, "Present", "Absent"), levels = c("Absent", "Present"))
  joined_data$elephant_present_possible <- factor(ifelse(joined_data$elephant_present_possible == 1, "Present", "Absent"), levels = c("Absent", "Present"))
  joined_data$elephant_present <- factor(ifelse(joined_data$elephant_present == 1, "Present", "Absent"), levels = c("Absent", "Present"))
  
  # Megafauna evolutionary history: Old World (continuous) vs New World (Pleistocene extinction)
  joined_data$megafaunaHistory <- factor(
    ifelse(joined_data$basin %in% c("Congo", "SE_Asia"), "OldWorld", "NewWorld"),
    levels = c("NewWorld", "OldWorld")
  )
  
  # Spatial Homogeneity definition (inverse of standard error, normalized to 0-1)
  raw_homo <- 1 / (joined_data$uoi_se + reg_uoi)
  if (length(raw_homo) > 1 && max(raw_homo) > min(raw_homo)) {
    joined_data$homogeneity <- (raw_homo - min(raw_homo)) / (max(raw_homo) - min(raw_homo))
  } else {
    joined_data$homogeneity <- 1.0
  }
  
  # Calculate Combined Weight (Spatial Precision * Temporal Alignment)
  joined_data <- joined_data %>%
    mutate(w_combined = w_uoi * w_temp_cluster) %>%
    mutate(w_combined_norm = w_combined / mean(w_combined))
  
  return(joined_data)
}

# --- 3. Dynamic Model Fitting Functions with RDS Formula Loading -------------
#' Fit Framework 1 Model (Beta Regression)
#'
#' Fits the best biophysical engineering model mapping GEDI UOI to mammal biomass
#' using a Beta regression with a logit link. Loads the formula dynamically.
#'
#' @param data A data frame returned by extract_scale_data.
#' @param formula_path Path to the RDS formula file.
#'
#' @return A fitted mgcv::gam object.
#' @export
fit_framework1_model <- function(data, formula_path = "outputs/framework1_best_formula.RDS") {
  formula_obj <- if (file.exists(formula_path)) {
    readRDS(formula_path)
  } else {
    uoi ~ B_H_index
  }
  m <- mgcv::gam(formula_obj, family = betar(link = "logit"), weights = w_combined_norm, data = data, method = "REML")
  check_model_convergence(m, "FW1 Calibration Model (REML)")
  m
}

#' Fit Framework 2 Model (Tweedie GLM)
#'
#' Fits the spaceborne detection model mapping mammal biomass to GEDI UOI and elevation
#' using a Tweedie GLM with a log link. Loads the formula dynamically.
#'
#' @param data A data frame returned by extract_scale_data.
#' @param formula_path Path to the RDS formula file.
#'
#' @return A fitted mgcv::gam object.
#' @export
fit_framework2_model <- function(data, formula_path = "outputs/framework2_best_formula.RDS") {
  formula_obj <- if (file.exists(formula_path)) {
    readRDS(formula_path)
  } else {
    B_H_index ~ uoi * basin + elevation
  }
  m <- mgcv::gam(formula_obj, family = tw(), weights = w_combined_norm, data = data, method = "REML")
  check_model_convergence(m, "FW2 Calibration Model (REML)")
  m
}
