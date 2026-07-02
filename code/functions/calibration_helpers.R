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

config_path <- if (file.exists("code/config.json")) "code/config.json" else "config.json"
config <- jsonlite::read_json(config_path)
MIN_TRAP_DAYS <- config$clustering$min_trap_days
CLUSTER_THRESHOLD_KM <- 11.1



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
  
  # Compute elephant presence dynamically based on spatial ranges or historical presence
  mcps$elephant_present_strict <- 0
  mcps$elephant_present_possible <- 0
  if (file.exists("outputs/elephant_ranges.gpkg")) {
    ele_ranges <- terra::vect("outputs/elephant_ranges.gpkg")
    
    # Strict: only extant (status == "Extant", which corresponds to presence == 1)
    ele_ranges_strict <- ele_ranges[ele_ranges$status == "Extant", ]
    if (nrow(ele_ranges_strict) > 0) {
      intersects_strict <- terra::is.related(mcps, ele_ranges_strict, "intersects")
      mcps$elephant_present_strict <- as.numeric(rowSums(as.matrix(intersects_strict)) > 0)
    }
    
    # Possible: all statuses (Extant, Possibly Extant, Possibly Extinct)
    intersects_possible <- terra::is.related(mcps, ele_ranges, "intersects")
    mcps$elephant_present_possible <- as.numeric(rowSums(as.matrix(intersects_possible)) > 0)
  }
  
  # Fallback to continent-level historical presence if GPKG is missing or has no intersection
  if (sum(mcps$elephant_present_strict, na.rm = TRUE) == 0) {
    mcps$elephant_present_strict <- ifelse(mcps$region %in% c("Congo", "SE_Asia"), 1, 0)
  }
  if (sum(mcps$elephant_present_possible, na.rm = TRUE) == 0) {
    mcps$elephant_present_possible <- ifelse(mcps$region %in% c("Congo", "SE_Asia"), 1, 0)
  }
  
  # Retain legacy alias for safety
  mcps$elephant_present <- mcps$elephant_present_possible

  
  mcps_congo <- mcps[mcps$region == "Congo", ]
  mcps_amazon <- mcps[mcps$region == "Amazon", ]
  mcps_seasia <- mcps[mcps$region == "SE_Asia", ]
  
  if (is.character(scale_m)) {
    r_congo_path <- sprintf("outputs/EOdata/analysis_stack_%s_Congo.tif", scale_m)
    r_amazon_path <- sprintf("outputs/EOdata/analysis_stack_%s_Amazon.tif", scale_m)
  } else {
    r_congo_path <- sprintf("outputs/EOdata/analysis_stack_%d_Congo.tif", scale_m)
    r_amazon_path <- sprintf("outputs/EOdata/analysis_stack_%d_Amazon.tif", scale_m)
  }
  
  # Fallback to synthetic data folders if needed
  if (!file.exists(r_congo_path)) {
    if (is.character(scale_m)) {
      r_congo_path <- sprintf("outputs/synthetic_EOdata/analysis_stack_%s_Congo.tif", scale_m)
    } else {
      r_congo_path <- sprintf("outputs/synthetic_EOdata/analysis_stack_%d_Congo.tif", scale_m)
    }
  }
  if (!file.exists(r_amazon_path)) {
    if (is.character(scale_m)) {
      r_amazon_path <- sprintf("outputs/synthetic_EOdata/analysis_stack_%s_Amazon.tif", scale_m)
    } else {
      r_amazon_path <- sprintf("outputs/synthetic_EOdata/analysis_stack_%d_Amazon.tif", scale_m)
    }
  }
  
  if (!file.exists(r_congo_path) || !file.exists(r_amazon_path)) {
    if (is.character(scale_m)) {
      stop(sprintf("GeoTIFF analysis stacks for scale %s are missing.", scale_m))
    } else {
      stop(sprintf("GeoTIFF analysis stacks for scale %d m are missing.", scale_m))
    }
  }
  
  r_congo <- terra::rast(r_congo_path)
  r_amazon <- terra::rast(r_amazon_path)
  
  if (is.character(scale_m) && scale_m == "native") {
    aggregate_names <- c("uoi", "uoi_sd", "rh98", "gedi_n", "elevation", "slope", "hnd", "precip", "clay", "forest_fraction", "Npp_median")
  } else {
    aggregate_names <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                         "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  }
  names(r_congo) <- aggregate_names
  names(r_amazon) <- aggregate_names
  
  # Extract Congo: touches = FALSE by default, touch fallback for empty polygons
  ext_congo <- terra::extract(r_congo, mcps_congo, df = TRUE, touches = FALSE)
  all_congo_ids <- 1:nrow(mcps_congo)
  extracted_congo_ids <- unique(ext_congo$ID)
  empty_congo_ids <- setdiff(all_congo_ids, extracted_congo_ids)
  
  if (length(empty_congo_ids) > 0) {
    ext_congo_touch <- terra::extract(r_congo, mcps_congo[empty_congo_ids, ], df = TRUE, touches = TRUE)
    ext_congo_touch$ID <- empty_congo_ids[ext_congo_touch$ID]
    ext_congo <- rbind(ext_congo %>% filter(ID %in% extracted_congo_ids), ext_congo_touch)
  }
  
  mcp_congo_df <- as.data.frame(mcps_congo)
  mcp_congo_df$ID <- 1:nrow(mcp_congo_df)
  
  if (is.character(scale_m) && scale_m == "native") {
    pixel_congo <- merge(ext_congo, mcp_congo_df, by = "ID") %>%
      filter(!is.na(uoi)) %>%
      select(-ID) %>%
      mutate(basin = "Congo", frip = NA)
  } else {
    pixel_congo <- merge(ext_congo, mcp_congo_df, by = "ID") %>%
      filter(!is.na(uoi)) %>%
      select(-ID) %>%
      mutate(basin = "Congo")
  }
  
  # Extract Amazon: touches = FALSE by default, touch fallback for empty polygons
  ext_amazon <- terra::extract(r_amazon, mcps_amazon, df = TRUE, touches = FALSE)
  all_amazon_ids <- 1:nrow(mcps_amazon)
  extracted_amazon_ids <- unique(ext_amazon$ID)
  empty_amazon_ids <- setdiff(all_amazon_ids, extracted_amazon_ids)
  
  if (length(empty_amazon_ids) > 0) {
    ext_amazon_touch <- terra::extract(r_amazon, mcps_amazon[empty_amazon_ids, ], df = TRUE, touches = TRUE)
    ext_amazon_touch$ID <- empty_amazon_ids[ext_amazon_touch$ID]
    ext_amazon <- rbind(ext_amazon %>% filter(ID %in% extracted_amazon_ids), ext_amazon_touch)
  }
  
  mcp_amazon_df <- as.data.frame(mcps_amazon)
  mcp_amazon_df$ID <- 1:nrow(mcp_amazon_df)
  
  if (is.character(scale_m) && scale_m == "native") {
    pixel_amazon <- merge(ext_amazon, mcp_amazon_df, by = "ID") %>%
      filter(!is.na(uoi)) %>%
      select(-ID) %>%
      mutate(basin = "Amazon", frip = NA)
  } else {
    pixel_amazon <- merge(ext_amazon, mcp_amazon_df, by = "ID") %>%
      filter(!is.na(uoi)) %>%
      select(-ID) %>%
      mutate(basin = "Amazon")
  }
  
  # Extract Southeast Asia if the 5000m file is present and there are SE Asia MCP polygons
  r_seasia_basename <- sprintf("analysis_stack_5000_SE_Asia.tif")
  r_seasia_path <- file.path("outputs", "EOdata", r_seasia_basename)
  
  if (file.exists(r_seasia_path) && nrow(mcps_seasia) > 0) {
    # If other scale requested, dynamically aggregate!
    if (!is.character(scale_m) && scale_m > 5000) {
      fact <- scale_m / 5000
      r_seasia <- terra::aggregate(rast(r_seasia_path), fact = fact, fun = "mean", na.rm = TRUE)
    } else {
      r_seasia <- rast(r_seasia_path)
    }
    
    names(r_seasia) <- aggregate_names
    
    ext_seasia <- terra::extract(r_seasia, mcps_seasia, df = TRUE, touches = FALSE)
    all_seasia_ids <- 1:nrow(mcps_seasia)
    extracted_seasia_ids <- unique(ext_seasia$ID)
    empty_seasia_ids <- setdiff(all_seasia_ids, extracted_seasia_ids)
    
    if (length(empty_seasia_ids) > 0) {
      ext_seasia_touch <- terra::extract(r_seasia, mcps_seasia[empty_seasia_ids, ], df = TRUE, touches = TRUE)
      ext_seasia_touch$ID <- empty_seasia_ids[ext_seasia_touch$ID]
      ext_seasia <- rbind(ext_seasia %>% filter(ID %in% extracted_seasia_ids), ext_seasia_touch)
    }
    
    mcp_seasia_df <- as.data.frame(mcps_seasia)
    mcp_seasia_df$ID <- 1:nrow(mcp_seasia_df)
    
    if (is.character(scale_m) && scale_m == "native") {
      pixel_seasia <- merge(ext_seasia, mcp_seasia_df, by = "ID") %>%
        filter(!is.na(uoi)) %>%
        select(-ID) %>%
        mutate(basin = "SE_Asia", frip = NA)
    } else {
      pixel_seasia <- merge(ext_seasia, mcp_seasia_df, by = "ID") %>%
        filter(!is.na(uoi)) %>%
        select(-ID) %>%
        mutate(basin = "SE_Asia")
    }
    
    pixel_data <- rbind(pixel_congo, pixel_amazon, pixel_seasia)
  } else {
    pixel_data <- rbind(pixel_congo, pixel_amazon)
  }
  
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
      uoi = mean(uoi, na.rm = TRUE),
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
  mgcv::gam(formula_obj, family = betar(link = "logit"), weights = w_combined_norm, data = data, method = "REML")
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
  mgcv::gam(formula_obj, family = tw(), weights = w_combined_norm, data = data, method = "REML")
}
