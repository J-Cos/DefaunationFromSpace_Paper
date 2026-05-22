# =============================================================================
# 08g_Representative_Scale_Maps.R
#
# Formally compares and maps predicted mammal biomass index and spatial uncertainty
# across three representative spatial scales:
#   1. MODIS Native Scale (~1 km, from analysis_stack_native_{basin}.tif)
#   2. 10 km GEDI Aggregate Scale (from analysis_stack_10000_{basin}.tif)
#   3. 25 km GEDI Aggregate Scale (from analysis_stack_25000_{basin}.tif)
#
# Generates two premium, publication-quality PNAS figures:
#   - outputs/congo_ct_gee_scatter_comparison_multi_res.png (1x3 scatter panels)
#   - outputs/congo_ct_gee_predictions_multi_res.png (3x4 landscape prediction master map)
# =============================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(cowplot)
library(dplyr)
library(readr)
library(mgcv)

cat("=== Starting Representative-Scale Biomass Modeling & Spatial Mapping ===\n\n")

# --- 1. Load Theme and GIS Polygons ------------------------------------------
source("code/functions/theme_pnas.R")

# Ensure outputs directory exists
dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

# Load robust cluster buffered MCPs
geojson_path <- "outputs/camera_traps_robust_buffered_mcps.geojson"
if (!file.exists(geojson_path)) {
  stop("Camera trap buffered MCPs GeoJSON missing: ", geojson_path, 
       "\nPlease run python3 code/visualise_camera_traps.py first to generate it.")
}
mcps <- terra::vect(geojson_path)
mcps_congo <- mcps[as.vector(values(mcps)$region) == "Congo", ]
mcps_amazon <- mcps[as.vector(values(mcps)$region) == "Amazon", ]
cat(sprintf("✓ Loaded %d Congo and %d Amazon study polygons.\n", nrow(mcps_congo), nrow(mcps_amazon)))

# Establish consistent cropped geographic study extents with 1.0 degree padding
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

# Set base font family from theme_pnas
base_family <- if (requireNamespace("showtext", quietly = TRUE)) "Roboto Condensed" else "sans"

# Config structure for representative scales
scale_configs <- list(
  list(label = "Native (~1 km)", suffix = "native", is_native = TRUE, color = "#1B5E20"),
  list(label = "10 km GEDI",     suffix = "10000",  is_native = FALSE, color = "#4527A0"),
  list(label = "25 km GEDI",     suffix = "25000",  is_native = FALSE, color = "#0D47A1")
)

scatter_panels <- list()
map_panels <- list()

# --- 2. Loop Through Scales --------------------------------------------------
for (cfg in scale_configs) {
  cat(sprintf("\nProcessing scale: %s...\n", cfg$label))
  
  r_congo_path <- sprintf("outputs/EOdata/analysis_stack_%s_Congo.tif", cfg$suffix)
  r_amazon_path <- sprintf("outputs/EOdata/analysis_stack_%s_Amazon.tif", cfg$suffix)
  
  if (!file.exists(r_congo_path) || !file.exists(r_amazon_path)) {
    stop(sprintf("GEOTIFFs missing for scale: %s\n  Congo: %s\n  Amazon: %s", cfg$label, r_congo_path, r_amazon_path))
  }
  
  r_congo <- rast(r_congo_path)
  r_amazon <- rast(r_amazon_path)
  
  # Conditionally assign names (Native stack has 11 layers, aggregate has 12 layers)
  if (cfg$is_native) {
    native_names <- c("uoi", "uoi_sd", "rh98", "gedi_n", "elevation", "slope", 
                      "hnd", "precip", "clay", "forest_fraction", "Npp_median")
    names(r_congo) <- native_names
    names(r_amazon) <- native_names
  } else {
    aggregate_names <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                         "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
    names(r_congo) <- aggregate_names
    names(r_amazon) <- aggregate_names
  }
  
  # --- 2.1 Extract Pixels & Merge With Camera Trap Cluster Metadata ----------
  cat("  Extracting pixel values inside MCP polygons...\n")
  extracted_congo <- terra::extract(r_congo, mcps_congo, df = TRUE, touches = TRUE)
  mcp_congo_df <- as.data.frame(mcps_congo)
  mcp_congo_df$ID <- 1:nrow(mcp_congo_df)
  pixel_congo <- merge(extracted_congo, mcp_congo_df, by = "ID") %>%
    filter(!is.na(uoi)) %>%
    select(-ID) %>%
    mutate(basin = "Congo")
  
  extracted_amazon <- terra::extract(r_amazon, mcps_amazon, df = TRUE, touches = TRUE)
  mcp_amazon_df <- as.data.frame(mcps_amazon)
  mcp_amazon_df$ID <- 1:nrow(mcp_amazon_df)
  pixel_amazon <- merge(extracted_amazon, mcp_amazon_df, by = "ID") %>%
    filter(!is.na(uoi)) %>%
    select(-ID) %>%
    mutate(basin = "Amazon")
  
  pixel_data <- rbind(pixel_congo, pixel_amazon)
  cat(sprintf("  Extracted %d Congo and %d Amazon pixels.\n", nrow(pixel_congo), nrow(pixel_amazon)))
  
  # --- 2.2 Aggregate to Polygon Means and Calculate Weights ------------------
  joined_data <- pixel_data %>%
    group_by(cluster_id, region, basin, trap_days, n_species, B_H_index, M_H_index, B_H_gt50, B_H_gt100, megafauna_fraction) %>%
    summarise(
      n_pixels = n(),
      uoi_sd = ifelse(is.na(sd(uoi, na.rm = TRUE)), 0, sd(uoi, na.rm = TRUE)),
      uoi = mean(uoi, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(trap_days >= 10)
  
  # Compute UOI standard error of the mean
  joined_data <- joined_data %>%
    mutate(uoi_se = uoi_sd / sqrt(n_pixels))
  
  # Calculate regularizing constant
  reg_uoi <- median(joined_data$uoi_se[joined_data$uoi_se > 0], na.rm = TRUE)
  if (is.na(reg_uoi) || reg_uoi == 0) reg_uoi <- 1e-4
  
  # Dual-weight precision framework
  joined_data <- joined_data %>%
    mutate(w_uoi = log10(trap_days) / (uoi_se + reg_uoi))
  
  # Weight normalization (mean = 1)
  joined_data$w_uoi_norm <- joined_data$w_uoi / mean(joined_data$w_uoi)
  
  # --- 2.3 Fit Weight-Normalized Tweedie GLM ----------------------------------
  tw_model <- gam(
    B_H_index ~ uoi,
    family = tw(),
    weights = w_uoi_norm,
    data = joined_data,
    method = "REML"
  )
  
  cat(sprintf("  Model fitted successfully for %s.\n", cfg$label))
  s_tw <- summary(tw_model)
  
  # --- 2.4 Build Scatter Plot Panel -------------------------------------------
  uoi_seq <- seq(min(joined_data$uoi) - 0.002, max(joined_data$uoi) + 0.002, length.out = 200)
  newdata <- data.frame(uoi = uoi_seq)
  pred <- predict(tw_model, newdata = newdata, type = "link", se.fit = TRUE)
  
  pred_df <- data.frame(
    uoi = uoi_seq,
    fit = exp(pred$fit),
    lower = exp(pred$fit - 1.96 * pred$se.fit),
    upper = exp(pred$fit + 1.96 * pred$se.fit)
  )
  
  p_val <- s_tw$p.pv["uoi"]
  dev_expl <- s_tw$dev.expl * 100
  p_text <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)
  
  stat_label <- sprintf(
    "%s\nDev. expl. = %.1f%%\n%s (N = %d)",
    cfg$label, dev_expl, p_text, nrow(joined_data)
  )
  
  p_scatter <- ggplot() +
    # 95% Confidence ribbon
    geom_ribbon(
      data = pred_df,
      aes(x = uoi, ymin = lower, ymax = upper),
      fill = cfg$color,
      alpha = 0.12
    ) +
    # Regression line
    geom_line(
      data = pred_df,
      aes(x = uoi, y = fit),
      color = cfg$color,
      linewidth = 0.7
    ) +
    # Empirical data points
    geom_point(
      data = joined_data,
      aes(x = uoi, y = B_H_index, size = trap_days,
          fill = uoi_se, shape = basin),
      color = "black",
      stroke = 0.3,
      alpha = 0.85
    ) +
    scale_shape_manual(
      values = c("Congo" = 21, "Amazon" = 24),
      name = "Basin"
    ) +
    scale_size_continuous(
      range = c(1.5, 6.0),
      guide = "none"
    ) +
    scale_fill_viridis_c(
      option = "magma",
      name = "UOI Spatial SE",
      guide = guide_colorbar(
        title.position = "top",
        barwidth = unit(1.5, "cm"),
        barheight = unit(0.12, "cm")
      )
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = stat_label,
      hjust = 1.05, vjust = 1.2,
      family = base_family,
      size = 2.0,
      fontface = "bold",
      color = "grey20"
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    theme_pnas(base_size = 7) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.size = unit(0.25, "cm"),
      legend.title = element_text(size = 5.5, face = "bold"),
      legend.text = element_text(size = 5.0),
      legend.margin = margin(t = -5, r = 0, b = 0, l = 0, unit = "pt")
    ) +
    labs(
      x = "GEDI Understory Openness Index (UOI)",
      y = "Total Biomass Index"
    )
  
  scatter_panels[[cfg$suffix]] <- p_scatter
  
  # --- 2.5 Run Geographic Spatial Predictions --------------------------------
  cat("  Running pixel-level spatial predictions for Congo...\n")
  r_congo_cropped <- crop(r_congo, study_extent_congo)
  congo_cells <- as.data.frame(r_congo_cropped[["uoi"]], cells = TRUE, xy = TRUE, na.rm = TRUE)
  names(congo_cells)[names(congo_cells) == "uoi"] <- "uoi"
  
  pred_congo <- predict(tw_model, newdata = congo_cells, type = "link", se.fit = TRUE)
  congo_cells$pred <- exp(pred_congo$fit)
  congo_cells$lower <- exp(pred_congo$fit - 1.96 * pred_congo$se.fit)
  congo_cells$upper <- exp(pred_congo$fit + 1.96 * pred_congo$se.fit)
  congo_cells$ci_range <- congo_cells$upper - congo_cells$lower
  
  r_pred_congo <- rast(r_congo_cropped[["uoi"]])
  names(r_pred_congo) <- "pred"
  values(r_pred_congo) <- NA
  r_pred_congo[congo_cells$cell] <- as.vector(congo_cells$pred)
  
  r_unc_congo <- rast(r_congo_cropped[["uoi"]])
  names(r_unc_congo) <- "ci_range"
  values(r_unc_congo) <- NA
  r_unc_congo[congo_cells$cell] <- as.vector(congo_cells$ci_range)
  
  cat("  Running pixel-level spatial predictions for Amazon...\n")
  r_amazon_cropped <- crop(r_amazon, study_extent_amazon)
  amazon_cells <- as.data.frame(r_amazon_cropped[["uoi"]], cells = TRUE, xy = TRUE, na.rm = TRUE)
  names(amazon_cells)[names(amazon_cells) == "uoi"] <- "uoi"
  
  pred_amazon <- predict(tw_model, newdata = amazon_cells, type = "link", se.fit = TRUE)
  amazon_cells$pred <- exp(pred_amazon$fit)
  amazon_cells$lower <- exp(pred_amazon$fit - 1.96 * pred_amazon$se.fit)
  amazon_cells$upper <- exp(pred_amazon$fit + 1.96 * pred_amazon$se.fit)
  amazon_cells$ci_range <- amazon_cells$upper - amazon_cells$lower
  
  r_pred_amazon <- rast(r_amazon_cropped[["uoi"]])
  names(r_pred_amazon) <- "pred"
  values(r_pred_amazon) <- NA
  r_pred_amazon[amazon_cells$cell] <- as.vector(amazon_cells$pred)
  
  r_unc_amazon <- rast(r_amazon_cropped[["uoi"]])
  names(r_unc_amazon) <- "ci_range"
  values(r_unc_amazon) <- NA
  r_unc_amazon[amazon_cells$cell] <- as.vector(amazon_cells$ci_range)
  
  # --- 2.6 Build Map Panels --------------------------------------------------
  make_map_panel <- function(r_data, mcps_vector, palette_option, title, legend_title) {
    ggplot() +
      geom_spatraster(data = r_data) +
      scale_fill_viridis_c(
        option = palette_option,
        name = legend_title,
        limits = c(0, 5000),
        oob = scales::squish,
        na.value = "transparent",
        guide = guide_colorbar(
          title.position = "top",
          barwidth = unit(1.0, "cm"),
          barheight = unit(0.08, "cm")
        )
      ) +
      geom_spatvector(data = mcps_vector, fill = NA, color = "black", linewidth = 0.2, alpha = 0.8) +
      theme_pnas(base_size = 5) +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 4.0, face = "bold"),
        legend.text = element_text(size = 3.5),
        legend.key.height = unit(0.15, "cm"),
        legend.key.width = unit(0.06, "cm"),
        legend.margin = margin(l = -2, r = 0, t = 0, b = 0, unit = "pt"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 5.5, face = "bold", margin = margin(b = 1)),
        plot.subtitle = element_blank()
      ) +
      labs(title = title)
  }
  
  # Label prefix assignment based on resolution index
  r_idx <- ifelse(cfg$suffix == "native", "Native", ifelse(cfg$suffix == "10000", "10km", "25km"))
  
  p_congo_pred <- make_map_panel(r_pred_congo, mcps_congo, "inferno", sprintf("Congo Predicted Biomass (%s)", r_idx), "Biomass")
  p_congo_unc  <- make_map_panel(r_unc_congo,  mcps_congo, "mako",    sprintf("Congo 95%% CI Uncertainty (%s)", r_idx), "CI Range")
  p_amazon_pred <- make_map_panel(r_pred_amazon, mcps_amazon, "inferno", sprintf("Amazon Predicted Biomass (%s)", r_idx), "Biomass")
  p_amazon_unc  <- make_map_panel(r_unc_amazon,  mcps_amazon, "mako",    sprintf("Amazon 95%% CI Uncertainty (%s)", r_idx), "CI Range")
  
  map_panels[[sprintf("%s_congo_pred", cfg$suffix)]] <- p_congo_pred
  map_panels[[sprintf("%s_congo_unc",  cfg$suffix)]] <- p_congo_unc
  map_panels[[sprintf("%s_amazon_pred", cfg$suffix)]] <- p_amazon_pred
  map_panels[[sprintf("%s_amazon_unc",  cfg$suffix)]] <- p_amazon_unc
}

# --- 3. Assemble and Export Comparative Scatter Figures -----------------------
cat("\nAssembling comparative scatter plot (1x3 panel)...\n")
fig_scatter <- cowplot::plot_grid(
  scatter_panels[["native"]]  + labs(title = "A. Native Resolution (~1 km)"),
  scatter_panels[["10000"]]  + labs(title = "B. 10 km Aggregate GEDI"),
  scatter_panels[["25000"]]  + labs(title = "C. 25 km Aggregate GEDI"),
  ncol = 3,
  align = "vh"
)

save_pnas(
  plot = fig_scatter,
  filename = "outputs/congo_ct_gee_scatter_comparison_multi_res.png",
  type = "double",
  height_cm = 6.5
)
cat("✓ Saved premium 1x3 scatter plot to outputs/congo_ct_gee_scatter_comparison_multi_res.png\n")

# --- 4. Assemble and Export 3x4 Landscape Master Map -------------------------
cat("Assembling comparative 3x4 spatial prediction map...\n")

# Lay out the grid row by row
fig_master_maps <- cowplot::plot_grid(
  # Row 1: Native
  map_panels[["native_congo_pred"]]  + labs(title = "A. Congo Expected Biomass (~1 km)"),
  map_panels[["native_congo_unc"]]   + labs(title = "B. Congo 95% CI Uncertainty (~1 km)"),
  map_panels[["native_amazon_pred"]]  + labs(title = "C. Amazon Expected Biomass (~1 km)"),
  map_panels[["native_amazon_unc"]]   + labs(title = "D. Amazon 95% CI Uncertainty (~1 km)"),
  
  # Row 2: 10km
  map_panels[["10000_congo_pred"]]   + labs(title = "E. Congo Expected Biomass (10 km)"),
  map_panels[["10000_congo_unc"]]    + labs(title = "F. Congo 95% CI Uncertainty (10 km)"),
  map_panels[["10000_amazon_pred"]]   + labs(title = "G. Amazon Expected Biomass (10 km)"),
  map_panels[["10000_amazon_unc"]]    + labs(title = "H. Amazon 95% CI Uncertainty (10 km)"),
  
  # Row 3: 25km
  map_panels[["25000_congo_pred"]]   + labs(title = "I. Congo Expected Biomass (25 km)"),
  map_panels[["25000_congo_unc"]]    + labs(title = "J. Congo 95% CI Uncertainty (25 km)"),
  map_panels[["25000_amazon_pred"]]   + labs(title = "K. Amazon Expected Biomass (25 km)"),
  map_panels[["25000_amazon_unc"]]    + labs(title = "L. Amazon 95% CI Uncertainty (25 km)"),
  
  ncol = 4,
  align = "hv",
  axis = "tblr"
)

save_pnas(
  plot = fig_master_maps,
  filename = "outputs/congo_ct_gee_predictions_multi_res.png",
  type = "double",
  height_cm = 21.0
)
cat("✓ Saved premium 3x4 prediction map to outputs/congo_ct_gee_predictions_multi_res.png\n")

cat("\n=== Representative-Scale Spatial Mapping Analysis Complete! ===\n")
