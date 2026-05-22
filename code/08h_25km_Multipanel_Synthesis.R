# =============================================================================
# 08h_25km_Multipanel_Synthesis.R
#
# Generates a premium, publication-quality multipanel figure for the 25km scale:
#   - Panel A: Scatterplot of empirical data, fitted weight-normalized Sigmoid model,
#              and robust bootstrapped extrapolation up to UOI = 1.0.
#   - Panel B: Boxplots of GEDI UOI pixels inside each buffered MCP polygon with
#              individual pixel jitter overlaid, ordered by cluster biomass index
#              and colored by continent.
#   - Panel C: Spatial maps for each continent:
#              - GEDI UOI (25km) with highly informative viridis palette (0.91 - 0.97)
#              - Predicted Biomass Index (25km) from Sigmoid model with winsorized
#                inferno palette (0 - 2000)
#
# Saves figure to:
#   - outputs/congo_ct_gee_25km_multipanel_synthesis.png
# Copies to brain artifacts folder:
#   - /home/j/.gemini/antigravity/brain/913e5cea-7c99-4b21-8124-ea8455da8457/
# =============================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(cowplot)
library(dplyr)
library(readr)
library(mgcv)
library(scales)

cat("=== Starting 25km Multipanel Synthesis Figure Generation ===\n\n")

# --- 1. Load Theme, Colors, and GIS Polygons ---------------------------------
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

base_family <- if (requireNamespace("showtext", quietly = TRUE)) "Roboto Condensed" else "sans"

# --- 2. Load 25 km Spatial Stacks --------------------------------------------
r_congo_path <- "outputs/EOdata/analysis_stack_25000_Congo.tif"
r_amazon_path <- "outputs/EOdata/analysis_stack_25000_Amazon.tif"

if (!file.exists(r_congo_path) || !file.exists(r_amazon_path)) {
  stop("25km GEE geotiffs missing. Please run GEE aggregation first.")
}

r_congo <- rast(r_congo_path)
r_amazon <- rast(r_amazon_path)

aggregate_names <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                     "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
names(r_congo) <- aggregate_names
names(r_amazon) <- aggregate_names

# --- 3. Extract Pixel Values and Merge Metadata -----------------------------
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

# Aggregate to Polygon Means for modeling
joined_data <- pixel_data %>%
  group_by(cluster_id, region, basin, trap_days, n_species, B_H_index, M_H_index, B_H_gt50, B_H_gt100, megafauna_fraction) %>%
  summarise(
    n_pixels = n(),
    uoi_sd = ifelse(is.na(sd(uoi, na.rm = TRUE)), 0, sd(uoi, na.rm = TRUE)),
    uoi = mean(uoi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(trap_days >= 10)

# Compute standard error of the mean
joined_data <- joined_data %>%
  mutate(uoi_se = uoi_sd / sqrt(n_pixels))

# Regularizing constant
reg_uoi <- median(joined_data$uoi_se[joined_data$uoi_se > 0], na.rm = TRUE)
if (is.na(reg_uoi) || reg_uoi == 0) reg_uoi <- 1e-4

# Dual-weight precision framework
joined_data <- joined_data %>%
  mutate(w_uoi = log10(trap_days) / (uoi_se + reg_uoi))

# Weight normalization (mean = 1)
joined_data$w_uoi_norm <- joined_data$w_uoi / mean(joined_data$w_uoi)

# --- 4. Fit Weight-Normalized Sigmoid Model (25 km) --------------------------
# Define Sigmoid (3-parameter logistic) function
sigmoid <- function(x, L, k, x0) {
  L / (1 + exp(-k * (x - x0)))
}

# Weighted Least Squares Objective Function
wls_obj <- function(params, x_val, y_val, w_val) {
  L <- params[1]
  k <- params[2]
  x0 <- params[3]
  pred <- sigmoid(x_val, L, k, x0)
  sum(w_val * (y_val - pred)^2)
}

init_params <- c(L = 4000, k = 100, x0 = 0.94)
res_opt <- optim(
  init_params, wls_obj,
  x_val = joined_data$uoi,
  y_val = joined_data$B_H_index,
  w_val = joined_data$w_uoi_norm,
  method = "L-BFGS-B",
  lower = c(100, 10, 0.85),
  upper = c(10000, 500, 1.0)
)

L_fit <- res_opt$par[1]
k_fit <- res_opt$par[2]
x0_fit <- res_opt$par[3]

# Compute R-squared and Weighted R-squared
mean_y <- mean(joined_data$B_H_index)
w_mean_y <- sum(joined_data$w_uoi_norm * joined_data$B_H_index) / sum(joined_data$w_uoi_norm)

tss <- sum((joined_data$B_H_index - mean_y)^2)
rss <- sum((joined_data$B_H_index - sigmoid(joined_data$uoi, L_fit, k_fit, x0_fit))^2)
r2_val <- 1 - rss/tss

wtss <- sum(joined_data$w_uoi_norm * (joined_data$B_H_index - w_mean_y)^2)
wrss <- sum(joined_data$w_uoi_norm * (joined_data$B_H_index - sigmoid(joined_data$uoi, L_fit, k_fit, x0_fit))^2)
wr2_val <- 1 - wrss/wtss

cat(sprintf("✓ Sigmoid Model (25 km) fitted: Weighted R2 = %.2f%%, Unweighted R2 = %.2f%%\n", wr2_val * 100, r2_val * 100))
cat(sprintf("  Asymptote (L) = %.2f, Growth Rate (k) = %.2f, Midpoint (x0) = %.4f\n", L_fit, k_fit, x0_fit))

# --- 5. PANEL A: Scatter Plot and Model Extrapolation ------------------------
cat("  Generating Panel A (Scatterplot with sigmoid extrapolation to UOI=1.0)...\n")

# Bootstrap 500 times for robust non-linear confidence intervals
set.seed(42)
n_boot <- 500
uoi_seq <- seq(0.91, 1.00, length.out = 300)
boot_preds <- matrix(NA, nrow = n_boot, ncol = length(uoi_seq))

cat("  Bootstrapping sigmoid curve for 95% confidence ribbon...\n")
for (b in 1:n_boot) {
  boot_idx <- sample(1:nrow(joined_data), replace = TRUE)
  x_b <- joined_data$uoi[boot_idx]
  y_b <- joined_data$B_H_index[boot_idx]
  w_b <- joined_data$w_uoi_norm[boot_idx]
  w_b_norm <- w_b / mean(w_b)
  
  res_b <- tryCatch({
    optim(init_params, wls_obj, x_val = x_b, y_val = y_b, w_val = w_b_norm, method = "L-BFGS-B",
          lower = c(100, 10, 0.85), upper = c(10000, 500, 1.0))
  }, error = function(e) NULL)
  
  if (!is.null(res_b) && res_b$convergence == 0) {
    boot_preds[b, ] <- sigmoid(uoi_seq, res_b$par[1], res_b$par[2], res_b$par[3])
  } else {
    boot_preds[b, ] <- sigmoid(uoi_seq, L_fit, k_fit, x0_fit)
  }
}

pred_df <- data.frame(
  uoi = uoi_seq,
  fit = sigmoid(uoi_seq, L_fit, k_fit, x0_fit),
  lower = apply(boot_preds, 2, function(col) quantile(col, 0.025, na.rm = TRUE)),
  upper = apply(boot_preds, 2, function(col) quantile(col, 0.975, na.rm = TRUE))
)

stat_label <- sprintf(
  "Weighted Sigmoid (25 km)\nWeighted R² = %.1f%%\nN = %d clusters",
  wr2_val * 100, nrow(joined_data)
)

p_scatter <- ggplot() +
  # 95% Confidence ribbon (extended to 1.0)
  geom_ribbon(
    data = pred_df,
    aes(x = uoi, ymin = lower, ymax = upper),
    fill = "#0D47A1",
    alpha = 0.12
  ) +
  # Regression line
  geom_line(
    data = pred_df,
    aes(x = uoi, y = fit),
    color = "#0D47A1",
    linewidth = 0.8
  ) +
  # Empirical data points
  geom_point(
    data = joined_data,
    aes(x = uoi, y = B_H_index, size = trap_days,
        fill = uoi_se, shape = basin),
    color = "black",
    stroke = 0.4,
    alpha = 0.85
  ) +
  scale_shape_manual(
    values = c("Congo" = 21, "Amazon" = 24),
    name = "Continent"
  ) +
  scale_size_continuous(
    range = c(1.5, 6.0),
    name = "Effort (Days)"
  ) +
  scale_fill_viridis_c(
    option = "magma",
    name = "UOI Spatial SE",
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(1.8, "cm"),
      barheight = unit(0.12, "cm")
    )
  ) +
  annotate(
    "text",
    x = 0.912, y = 7000,
    label = stat_label,
    hjust = 0, vjust = 1,
    family = base_family,
    size = 2.4,
    fontface = "bold",
    color = "grey20"
  ) +
  scale_x_continuous(limits = c(0.91, 1.00), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 8500), expand = expansion(mult = c(0, 0.05))) +
  theme_pnas(base_size = 7.5) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 6.0, face = "bold"),
    legend.text = element_text(size = 5.5),
    legend.key.size = unit(0.28, "cm"),
    legend.margin = margin(t = -5, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  labs(
    title = "A. Sigmoid Model Fit & Extrapolation (25 km)",
    x = "GEDI Understory Openness Index (UOI)",
    y = "Total Mammal Biomass Index"
  )

# --- 6. PANEL B: Boxplots with Jittered Pixels -------------------------------
cat("  Generating Panel B (Boxplots with jittered pixels)...\n")

pixel_data_filtered <- pixel_data %>%
  filter(trap_days >= 10)

# Sort cluster_id by actual standing biomass index (B_H_index)
pixel_data_filtered$cluster_id <- reorder(pixel_data_filtered$cluster_id, pixel_data_filtered$B_H_index)

p_box <- ggplot(pixel_data_filtered, aes(x = cluster_id, y = uoi)) +
  # Semi-transparent boxplots colored by continent
  geom_boxplot(
    aes(fill = basin, color = basin),
    alpha = 0.22,
    linewidth = 0.3,
    outlier.shape = NA
  ) +
  # Jittered individual pixel points
  geom_jitter(
    aes(color = basin),
    width = 0.15,
    height = 0,
    size = 0.7,
    alpha = 0.65
  ) +
  scale_fill_manual(
    values = pal_basin,
    name = "Continent"
  ) +
  scale_color_manual(
    values = pal_basin,
    name = "Continent"
  ) +
  theme_pnas(base_size = 7.5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.0, family = base_family),
    legend.position = "bottom",
    legend.title = element_text(size = 6.0, face = "bold"),
    legend.text = element_text(size = 5.5),
    legend.margin = margin(t = -5, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  labs(
    title = "B. Within-Cluster Pixel Heterogeneity (25 km)",
    x = "Camera Trap Cluster (Ordered by Mammal Biomass Index)",
    y = "Pixel Understory Openness (UOI)"
  )

# --- 7. PANEL C: Spatial Maps of UOI & Predicted Biomass ----------------------
cat("  Generating Panel C Spatial Maps (25 km)...\n")

# --- 7.1 Congo Predictions and Mapping
cat("    Processing Congo spatial predictions...\n")
r_congo_cropped <- crop(r_congo, study_extent_congo)
congo_cells <- as.data.frame(r_congo_cropped[["uoi"]], cells = TRUE, xy = TRUE, na.rm = TRUE)
names(congo_cells)[names(congo_cells) == "uoi"] <- "uoi"

congo_cells$pred <- sigmoid(congo_cells$uoi, L_fit, k_fit, x0_fit)

r_pred_congo <- rast(r_congo_cropped[["uoi"]])
names(r_pred_congo) <- "pred"
values(r_pred_congo) <- NA
r_pred_congo[congo_cells$cell] <- as.vector(congo_cells$pred)

# --- 7.2 Amazon Predictions and Mapping
cat("    Processing Amazon spatial predictions...\n")
r_amazon_cropped <- crop(r_amazon, study_extent_amazon)
amazon_cells <- as.data.frame(r_amazon_cropped[["uoi"]], cells = TRUE, xy = TRUE, na.rm = TRUE)
names(amazon_cells)[names(amazon_cells) == "uoi"] <- "uoi"

amazon_cells$pred <- sigmoid(amazon_cells$uoi, L_fit, k_fit, x0_fit)

r_pred_amazon <- rast(r_amazon_cropped[["uoi"]])
names(r_pred_amazon) <- "pred"
values(r_pred_amazon) <- NA
r_pred_amazon[amazon_cells$cell] <- as.vector(amazon_cells$pred)

# --- 7.3 Map Panel Helper
make_map_panel <- function(r_data, mcps_vector, palette_option, title, legend_title, limits = NULL, winsorize = FALSE) {
  p <- ggplot() +
    geom_spatraster(data = r_data)
  
  if (winsorize) {
    p <- p + scale_fill_viridis_c(
      option = palette_option,
      name = legend_title,
      limits = limits,
      oob = scales::squish,
      na.value = "transparent",
      guide = guide_colorbar(
        title.position = "top",
        barwidth = unit(1.2, "cm"),
        barheight = unit(0.08, "cm")
      )
    )
  } else {
    p <- p + scale_fill_viridis_c(
      option = palette_option,
      name = legend_title,
      limits = limits,
      na.value = "transparent",
      guide = guide_colorbar(
        title.position = "top",
        barwidth = unit(1.2, "cm"),
        barheight = unit(0.08, "cm")
      )
    )
  }
  
  p <- p +
    geom_spatvector(data = mcps_vector, fill = NA, color = "black", linewidth = 0.25, alpha = 0.8) +
    theme_pnas(base_size = 5.5) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 4.5, face = "bold"),
      legend.text = element_text(size = 4.0),
      legend.key.height = unit(0.18, "cm"),
      legend.key.width = unit(0.08, "cm"),
      legend.margin = margin(l = -2, r = 0, t = 0, b = 0, unit = "pt"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 6.0, face = "bold", margin = margin(b = 2)),
      plot.subtitle = element_blank()
    ) +
    labs(title = title)
  
  return(p)
}

# Generate Map Sub-panels (UOI using viridis, Predictions using inferno)
# Optimized dynamic ranges based on empirical raster percentiles to maximize contrast:
# GEDI UOI is bounded between 0.91 and 0.97.
# Predicted Biomass spans 0 to 3311, with 99% of pixels below 1400.
uoi_lims <- c(0.91, 0.97)
pred_lims <- c(0, 2000)

p_congo_uoi  <- make_map_panel(r_congo_cropped[["uoi"]], mcps_congo, "viridis", "C.1 Congo GEDI UOI (25 km)", "UOI", limits = uoi_lims, winsorize = TRUE)
p_congo_pred <- make_map_panel(r_pred_congo,             mcps_congo, "inferno", "C.2 Congo Predicted Biomass (25 km)", "Biomass", limits = pred_lims, winsorize = TRUE)
p_amazon_uoi  <- make_map_panel(r_amazon_cropped[["uoi"]], mcps_amazon, "viridis", "C.3 Amazon GEDI UOI (25 km)", "UOI", limits = uoi_lims, winsorize = TRUE)
p_amazon_pred <- make_map_panel(r_pred_amazon,             mcps_amazon, "inferno", "C.4 Amazon Predicted Biomass (25 km)", "Biomass", limits = pred_lims, winsorize = TRUE)

# Lay out Panel C as a 2x2 grid
fig_maps <- cowplot::plot_grid(
  p_congo_uoi, p_congo_pred,
  p_amazon_uoi, p_amazon_pred,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

# --- 8. Assemble and Save Master Figure --------------------------------------
cat("\nAssembling final 25km Multipanel Master Figure...\n")

# Row 1: A (Scatterplot) and B (Boxplots) side-by-side (giving B more width)
fig_top_row <- cowplot::plot_grid(
  p_scatter,
  p_box,
  ncol = 2,
  rel_widths = c(1.0, 1.45),
  align = "h",
  axis = "tb"
)

# Combine Row 1 (Scatter & Box) and Row 2 (2x2 Maps)
fig_master <- cowplot::plot_grid(
  fig_top_row,
  fig_maps,
  ncol = 1,
  rel_heights = c(1.0, 1.7)
)

# Save using the standard PNAS exporter
save_pnas(
  plot = fig_master,
  filename = "outputs/congo_ct_gee_25km_multipanel_synthesis.png",
  type = "double",
  height_cm = 19.5
)
cat("✓ Saved premium master figure to outputs/congo_ct_gee_25km_multipanel_synthesis.png\n")

# Copy to the active brain artifacts folder
brain_artifacts_dir <- "/home/j/.gemini/antigravity/brain/913e5cea-7c99-4b21-8124-ea8455da8457"
file.copy("outputs/congo_ct_gee_25km_multipanel_synthesis.png", 
          file.path(brain_artifacts_dir, "congo_ct_gee_25km_multipanel_synthesis.png"), 
          overwrite = TRUE)
cat("✓ Copied master figure to brain artifacts folder.\n")

cat("\n=== 25km Multipanel Synthesis Complete! ===\n")
