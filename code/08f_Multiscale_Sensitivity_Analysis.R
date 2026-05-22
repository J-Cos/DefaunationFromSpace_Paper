# =============================================================================
# 08f_Multiscale_Sensitivity_Analysis.R
#
# Loops through all 20 remote sensing resolution stacks (5 km to 100 km) under
# outputs/EOdata/, extracts buffered camera trap MCP averages, fits dual-weighted
# linear regressions and weight-normalized Tweedie GLMs at each scale, and
# exports scale-dependency results.
#
# Generates a premium 4-panel (2x2) PNAS-style publication scaling figure:
#   outputs/congo_ct_gee_multiscale_scaling_plots.png
# =============================================================================

library(terra)
library(dplyr)
library(readr)
library(ggplot2)
library(scales)
library(cowplot)
library(mgcv)

cat("=== Starting Multi-Scale Spatial Sensitivity Analysis (5 km - 100 km) ===\n\n")

# --- 1. Load Theme, Polygons and Scales ---------------------------------------
source("code/functions/theme_pnas.R")

# Ensure outputs directory exists
dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

# Load mcps geojson
geojson_path <- "outputs/camera_traps_robust_buffered_mcps.geojson"
if (!file.exists(geojson_path)) {
  stop("Camera trap buffered MCPs missing. Please run code/visualise_camera_traps.py first.")
}
mcps <- terra::vect(geojson_path)
mcps_congo <- mcps[mcps$region == "Congo", ]
mcps_amazon <- mcps[mcps$region == "Amazon", ]

scales_seq <- seq(5000, 100000, by = 5000)
scale_results <- list()

# --- 2. Loop Through All Spatial Resolutions ----------------------------------
for (s in scales_seq) {
  cat(sprintf("Analyzing scale: %d meters (%d km)...\n", s, s / 1000))
  
  r_congo_path <- sprintf("outputs/EOdata/analysis_stack_%d_Congo.tif", s)
  r_amazon_path <- sprintf("outputs/EOdata/analysis_stack_%d_Amazon.tif", s)
  
  if (!file.exists(r_congo_path) || !file.exists(r_amazon_path)) {
    cat(sprintf("  Skipping scale %d (TIFFs missing)\n", s))
    next
  }
  
  r_congo <- rast(r_congo_path)
  r_amazon <- rast(r_amazon_path)
  
  names(r_congo) <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                      "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  names(r_amazon) <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                       "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  
  # Extract Congo pixels using touches = TRUE to prevent NA at coarse resolutions
  extracted_congo <- terra::extract(r_congo, mcps_congo, df = TRUE, touches = TRUE)
  mcp_congo_df <- as.data.frame(mcps_congo)
  mcp_congo_df$ID <- 1:nrow(mcp_congo_df)
  pixel_congo <- merge(extracted_congo, mcp_congo_df, by = "ID") %>%
    filter(!is.na(uoi) & !is.na(frip)) %>%
    select(-ID) %>%
    mutate(basin = "Congo")
  
  # Extract Amazon pixels
  extracted_amazon <- terra::extract(r_amazon, mcps_amazon, df = TRUE, touches = TRUE)
  mcp_amazon_df <- as.data.frame(mcps_amazon)
  mcp_amazon_df$ID <- 1:nrow(mcp_amazon_df)
  pixel_amazon <- merge(extracted_amazon, mcp_amazon_df, by = "ID") %>%
    filter(!is.na(uoi) & !is.na(frip)) %>%
    select(-ID) %>%
    mutate(basin = "Amazon")
  
  # Combine
  pixel_data <- rbind(pixel_congo, pixel_amazon)
  
  # Compute cluster-level polygon means
  joined_data <- pixel_data %>%
    group_by(cluster_id, region, basin, trap_days, n_species, B_H_index, M_H_index, B_H_gt50, B_H_gt100, megafauna_fraction) %>%
    summarise(
      n_pixels = n(),
      uoi_sd = ifelse(is.na(sd(uoi, na.rm = TRUE)), 0, sd(uoi, na.rm = TRUE)),
      frip_sd = ifelse(is.na(sd(frip, na.rm = TRUE)), 0, sd(frip, na.rm = TRUE)),
      uoi = mean(uoi, na.rm = TRUE),
      frip = mean(frip, na.rm = TRUE),
      frip_mk_tau = mean(frip_mk_tau, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(trap_days >= 10)
  
  if (nrow(joined_data) < 5) {
    cat(sprintf("  Scale %d has too few clusters with valid overlap, skipping...\n", s))
    next
  }
  
  # Compute standard errors
  joined_data <- joined_data %>%
    mutate(
      uoi_se = uoi_sd / sqrt(n_pixels),
      frip_se = frip_sd / sqrt(n_pixels)
    )
  
  reg_uoi <- median(joined_data$uoi_se[joined_data$uoi_se > 0])
  if (is.na(reg_uoi) || reg_uoi == 0) reg_uoi <- 1e-4
  
  reg_frip <- median(joined_data$frip_se[joined_data$frip_se > 0])
  if (is.na(reg_frip) || reg_frip == 0) reg_frip <- 1e-4
  
  # Re-compute composite dual weights
  joined_data <- joined_data %>%
    mutate(
      w_uoi = log10(trap_days) / (uoi_se + reg_uoi),
      w_frip = log10(trap_days) / (frip_se + reg_frip)
    )
  
  # Weighted Linear Regressions
  # 1. Total Biomass vs GEDI UOI
  lm_uoi <- lm(B_H_index ~ uoi, data = joined_data, weights = w_uoi)
  sum_lm_uoi <- summary(lm_uoi)
  
  r_uoi <- cov.wt(data.frame(joined_data$B_H_index, joined_data$uoi), wt = joined_data$w_uoi, cor = TRUE)$cor[1, 2]
  p_uoi <- sum_lm_uoi$coefficients["uoi", "Pr(>|t|)"]
  slope_uoi <- sum_lm_uoi$coefficients["uoi", "Estimate"]
  se_uoi <- sum_lm_uoi$coefficients["uoi", "Std. Error"]
  r2_uoi <- sum_lm_uoi$r.squared
  
  # 2. Total Biomass vs GEE FRIP
  lm_frip <- lm(B_H_index ~ frip, data = joined_data, weights = w_frip)
  sum_lm_frip <- summary(lm_frip)
  
  r_frip <- cov.wt(data.frame(joined_data$B_H_index, joined_data$frip), wt = joined_data$w_frip, cor = TRUE)$cor[1, 2]
  p_frip <- sum_lm_frip$coefficients["frip", "Pr(>|t|)"]
  
  # 3. Fit weight-normalized Tweedie GLM (using mgcv::gam with Tweedie family)
  joined_data$w_uoi_norm <- joined_data$w_uoi / mean(joined_data$w_uoi)
  tw_fit <- gam(B_H_index ~ uoi, data = joined_data, weights = w_uoi_norm, family = Tweedie(p = 1.5, link = "log"))
  sum_tw <- summary(tw_fit)
  
  p_tw_slope <- sum_tw$p.pv["uoi"]
  dev_expl_tw <- sum_tw$dev.expl
  
  scale_results[[as.character(s)]] <- data.frame(
    scale_m = s,
    scale_km = s / 1000,
    n_clusters = nrow(joined_data),
    uoi_r = r_uoi,
    uoi_slope = slope_uoi,
    uoi_se = se_uoi,
    uoi_p = p_uoi,
    uoi_r2 = r2_uoi,
    frip_r = r_frip,
    frip_p = p_frip,
    tweedie_p = p_tw_slope,
    tweedie_dev_expl = dev_expl_tw
  )
}

# --- 3. Save Scale Dependency CSV ---------------------------------------------
results_df <- do.call(rbind, scale_results)
write_csv(results_df, "outputs/multiscale_spatial_sensitivity_results.csv")
cat("\n✓ Multi-scale spatial results saved to outputs/multiscale_spatial_sensitivity_results.csv\n\n")

# --- 4. Build ggplot Scaling Panels -------------------------------------------
cat("Building premium 4-panel spatial scaling plots...\n")

# Retrieve font
base_family <- if (requireNamespace("showtext", quietly = TRUE)) "Helvetica Neue" else "sans"

# Panel A: GEDI UOI Correlation vs. Scale
p_a <- ggplot(results_df, aes(x = scale_km, y = uoi_r)) +
  geom_line(color = "#1B5E20", linewidth = 0.8) +
  geom_point(color = "#2E7D32", size = 1.8, shape = 21, fill = "white", stroke = 1.0) +
  scale_x_continuous(breaks = seq(0, 100, by = 20), limits = c(5, 100)) +
  scale_y_continuous(limits = c(0.2, 0.65), breaks = seq(0.2, 0.6, by = 0.1)) +
  labs(
    x = "Spatial Aggregation Scale (km)",
    y = "Weighted Pearson Correlation (r)",
    title = "A. GEDI UOI vs. Total Biomass Index Correlation",
    subtitle = "Relationship strength across spatial resolutions"
  ) +
  theme_pnas(base_size = 8)

# Panel B: GEE FRIP Correlation vs. Scale
p_b <- ggplot(results_df, aes(x = scale_km, y = frip_r)) +
  geom_line(color = "#BF360C", linewidth = 0.8) +
  geom_point(color = "#FF5722", size = 1.8, shape = 24, fill = "white", stroke = 1.0) +
  scale_x_continuous(breaks = seq(0, 100, by = 20), limits = c(5, 100)) +
  scale_y_continuous(limits = c(-0.1, 0.3), breaks = seq(-0.1, 0.3, by = 0.1)) +
  geom_hline(yintercept = 0.0, linetype = "dotted", color = "grey50", linewidth = 0.3) +
  labs(
    x = "Spatial Aggregation Scale (km)",
    y = "Weighted Pearson Correlation (r)",
    title = "B. GEE FRIP vs. Total Biomass Index Correlation",
    subtitle = "Decoupled relationship strength across spatial resolutions"
  ) +
  theme_pnas(base_size = 8)

# Panel C: Statistical Significance (p-values) vs. Scale
p_c <- ggplot(results_df) +
  geom_line(aes(x = scale_km, y = uoi_p, color = "Weighted Linear"), linewidth = 0.8) +
  geom_point(aes(x = scale_km, y = uoi_p, color = "Weighted Linear"), size = 1.5, shape = 21, fill = "white", stroke = 0.8) +
  geom_line(aes(x = scale_km, y = tweedie_p, color = "Tweedie GLM"), linewidth = 0.8) +
  geom_point(aes(x = scale_km, y = tweedie_p, color = "Tweedie GLM"), size = 1.5, shape = 22, fill = "white", stroke = 0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#D32F2F", linewidth = 0.4) +
  geom_hline(yintercept = 0.01, linetype = "dashed", color = "#388E3C", linewidth = 0.4) +
  annotate("text", x = 90, y = 0.065, label = "p = 0.05", color = "#D32F2F", size = 2.0, fontface = "bold") +
  annotate("text", x = 90, y = 0.013, label = "p = 0.01", color = "#388E3C", size = 2.0, fontface = "bold") +
  scale_x_continuous(breaks = seq(0, 100, by = 20), limits = c(5, 100)) +
  scale_y_log10(
    labels = trans_format("log10", math_format(10^.x)),
    limits = c(0.0001, 1.0),
    breaks = c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5, 1.0)
  ) +
  scale_color_manual(
    name = "Model Form",
    values = c("Weighted Linear" = "#43A047", "Tweedie GLM" = "#8E24AA")
  ) +
  labs(
    x = "Spatial Aggregation Scale (km)",
    y = "Slope Significance (p-value, Log Scale)",
    title = "C. UOI Slope Statistical Significance",
    subtitle = "Scale dependency of p-values (alpha threshold benchmarks)"
  ) +
  theme_pnas(base_size = 8) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.size = unit(0.35, "cm"),
    legend.margin = margin(t = -5, r = 0, b = 0, l = 0, unit = "pt")
  )

# Panel D: Tweedie Deviance Explained (%) vs. Scale
p_d <- ggplot(results_df, aes(x = scale_km, y = tweedie_dev_expl * 100)) +
  geom_line(color = "#8E24AA", linewidth = 0.8) +
  geom_point(color = "#BA68C8", size = 1.8, shape = 22, fill = "white", stroke = 1.0) +
  scale_x_continuous(breaks = seq(0, 100, by = 20), limits = c(5, 100)) +
  scale_y_continuous(limits = c(10, 55), breaks = seq(10, 50, by = 10)) +
  labs(
    x = "Spatial Aggregation Scale (km)",
    y = "Deviance Explained (%)",
    title = "D. Tweedie GLM Model Performance",
    subtitle = "Deviance explained as a function of aggregation resolution"
  ) +
  theme_pnas(base_size = 8)

# --- 5. Assemble and Export PNAS Figure 10 ------------------------------------
fig_scaling <- cowplot::plot_grid(
  p_a, p_b,
  p_c, p_d,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

save_pnas(
  plot = fig_scaling,
  filename = "outputs/congo_ct_gee_multiscale_scaling_plots.png",
  type = "double",
  height_cm = 16.0
)

cat("✓ Saved premium 4-panel figure: outputs/congo_ct_gee_multiscale_scaling_plots.png\n")
cat("=== Multi-Scale Spatial Sensitivity Analysis Complete ===\n")
