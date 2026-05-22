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

# --- Internal Helper: Calculate Camera Trap Deployment Temporal Weights ------
calculate_temporal_weights <- function() {
  detections_path <- "outputs/camera_traps_joint_detections.csv"
  if (!file.exists(detections_path)) {
    return(NULL)
  }
  
  det_all <- readr::read_csv(detections_path, show_col_types = FALSE)
  det_all$start_date <- as.Date(det_all$start_date)
  det_all$end_date <- as.Date(det_all$end_date)
  
  # Haversine distance single-linkage clustering (threshold = 11.1 km)
  haversine_dist <- function(lon1, lat1, lon2, lat2) {
    r <- 6371.0
    rad <- pi / 180
    dlon <- (lon2 - lon1) * rad
    dlat <- (lat2 - lat1) * rad
    lat1 <- lat1 * rad
    lat2 <- lat2 * rad
    a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
    c <- 2 * asin(sqrt(a))
    return(r * c)
  }
  
  coords_df <- det_all %>% 
    select(region, longitude, latitude) %>% 
    distinct() %>% 
    mutate(cluster_id_geo = "")
  
  for (reg in unique(coords_df$region)) {
    sub_indices <- which(coords_df$region == reg)
    sub <- coords_df[sub_indices, ]
    n <- nrow(sub)
    if (n == 0) next
    if (n > 1) {
      dist_mat <- matrix(0, nrow=n, ncol=n)
      for (i in 1:n) {
        for (j in 1:n) {
          dist_mat[i,j] <- haversine_dist(sub$longitude[i], sub$latitude[i], sub$longitude[j], sub$latitude[j])
        }
      }
      hc <- hclust(as.dist(dist_mat), method="single")
      labels <- cutree(hc, h=11.1)
    } else {
      labels <- 1
    }
    coords_df$cluster_id_geo[sub_indices] <- paste0(reg, "_", sprintf("%02d", labels))
  }
  
  det_all <- det_all %>%
    left_join(coords_df, by = c("region", "longitude", "latitude"))
  
  deployments <- det_all %>%
    select(region, cluster_id_geo, deployment_id, start_date, end_date, trap_days) %>%
    distinct()
  
  gedi_start <- as.Date("2019-04-17")
  
  deployments <- deployments %>%
    mutate(
      years_before_gedi = as.numeric(gedi_start - start_date) / 365.25,
      w_temp = case_when(
        years_before_gedi <= 1.0  ~ 1.0,
        years_before_gedi <= 6.0  ~ 0.5,
        years_before_gedi <= 11.0 ~ 0.25,
        TRUE                      ~ 0.1
      )
    )
  
  cluster_temp_metrics <- deployments %>%
    group_by(cluster_id_geo) %>%
    summarise(
      w_temp_cluster = sum(trap_days * w_temp) / sum(trap_days),
      .groups = "drop"
    )
  
  return(cluster_temp_metrics)
}

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
  
  mcps_congo <- mcps[mcps$region == "Congo", ]
  mcps_amazon <- mcps[mcps$region == "Amazon", ]
  
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
      filter(!is.na(uoi) & !is.na(frip)) %>%
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
      filter(!is.na(uoi) & !is.na(frip)) %>%
      select(-ID) %>%
      mutate(basin = "Amazon")
  }
  
  pixel_data <- rbind(pixel_congo, pixel_amazon)
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
    group_by(cluster_id, region, basin, trap_days, n_species, B_H_index, M_H_index, B_H_gt50, B_H_gt100, megafauna_fraction) %>%
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
      .groups = "drop"
    ) %>%
    filter(trap_days >= 10)
  
  # Compute standard errors of the mean
  joined_data <- joined_data %>%
    mutate(uoi_se = uoi_sd / sqrt(n_pixels))
  
  # Compute precision weights
  reg_uoi <- median(joined_data$uoi_se[joined_data$uoi_se > 0])
  if (is.na(reg_uoi) || reg_uoi == 0) reg_uoi <- 1e-4
  
  joined_data <- joined_data %>%
    mutate(w_uoi = log10(trap_days) / (uoi_se + reg_uoi))
  
  # Normalize weights so they sum to N (mean = 1) for proper statistical scale
  joined_data$w_uoi_norm <- joined_data$w_uoi / mean(joined_data$w_uoi)
  joined_data$basin <- factor(joined_data$basin, levels = c("Amazon", "Congo"))
  
  # Spatial Homogeneity definition (inverse of standard error, normalized to 0-1)
  raw_homo <- 1 / (joined_data$uoi_se + reg_uoi)
  if (length(raw_homo) > 1 && max(raw_homo) > min(raw_homo)) {
    joined_data$homogeneity <- (raw_homo - min(raw_homo)) / (max(raw_homo) - min(raw_homo))
  } else {
    joined_data$homogeneity <- 1.0
  }
  
  # Join temporal weights back to spatial dataset
  cluster_temp <- calculate_temporal_weights()
  if (!is.null(cluster_temp)) {
    joined_data <- joined_data %>%
      left_join(cluster_temp, by = c("cluster_id" = "cluster_id_geo"))
    
    # Fallback for NAs
    joined_data$w_temp_cluster[is.na(joined_data$w_temp_cluster)] <- median(joined_data$w_temp_cluster, na.rm = TRUE)
  } else {
    joined_data$w_temp_cluster <- 1.0 # Default fallback
  }
  
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
  mgcv::gam(formula_obj, family = betar(link = "logit"), weights = w_uoi_norm, data = data, method = "REML")
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
    B_H_index ~ uoi:basin + elevation
  }
  mgcv::gam(formula_obj, family = tw(), weights = w_uoi_norm, data = data, method = "REML")
}
