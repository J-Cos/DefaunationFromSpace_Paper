# =============================================================================
# 08_Congo_Camera_Trap_Comparison.R
#
# Analyzes and visualizes GEE 5000m Congo remote sensing data (GEDI, FRIP) 
# against robust camera trap metrics for the Congo Basin (clusters >= 10 trap-days).
#
# Incorporates spatial polygon extraction:
#   - Loads the Python-generated robust cluster MCPs (+ 5.5 km buffer) GeoJSON.
#   - Extracts all raster pixels falling inside the polygons to visualize
#     spatial heterogeneity with publication-quality PNAS-style box plots.
#   - Calculates polygon-wide average remote sensing signatures for robust regressions.
#
# Saves high-resolution, publication-quality PNAS-style figures:
#   1. outputs/congo_ct_gee_spatial_maps.png (2 panels: spatial maps with MCPs)
#   2. outputs/congo_ct_gee_spatial_distributions.png (2 panels: pixel box plots)
#   3. outputs/congo_ct_gee_statistical_plots.png (5 panels: scatter plots with effort weighting)
# =============================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(cowplot)
library(dplyr)
library(readr)
library(stringr)

cat("=== Starting Congo Camera Trap & GEE Joint Analysis (MCP Buffers) ===\n\n")

# --- 1. Load Theme & Configurations -----------------------------------------
source("code/functions/theme_pnas.R")

# Ensure outputs directory exists
dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

# --- 2. Load Camera Trap Datasets -------------------------------------------
# Load the robust cluster buffered MCPs
geojson_path <- "outputs/camera_traps_robust_buffered_mcps.geojson"
if (!file.exists(geojson_path)) {
  stop("Python robust cluster buffered MCPs GeoJSON missing: ", geojson_path, 
       "\nPlease run python3 code/visualise_camera_traps.py first to generate it.")
}

mcps <- terra::vect(geojson_path)
mcps_congo <- mcps[mcps$region == "Congo", ]
mcps_amazon <- mcps[mcps$region == "Amazon", ]
cat(sprintf("✓ Loaded %d Congo and %d Amazon buffered MCP polygons from GeoJSON.\n", nrow(mcps_congo), nrow(mcps_amazon)))

# Load individual camera trap locations to plot on maps
ct_detections_path <- "outputs/camera_traps_joint_detections.csv"
if (!file.exists(ct_detections_path)) {
  stop("Camera trap detections file missing: ", ct_detections_path)
}

ct_all <- read_csv(ct_detections_path, show_col_types = FALSE)
ct_congo <- ct_all %>%
  filter(region == "Congo") %>%
  distinct(project_id, deployment_id, longitude, latitude, trap_days)
ct_amazon <- ct_all %>%
  filter(region == "Amazon") %>%
  distinct(project_id, deployment_id, longitude, latitude, trap_days)

cat(sprintf("✓ Loaded %d Congo and %d Amazon camera trap coordinates.\n", nrow(ct_congo), nrow(ct_amazon)))

# --- 3. Load GEE 5000m Congo & Amazon Raster Stacks -------------------------
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

cat("✓ Loaded Congo and Amazon GEE raster stacks.\n")

# --- 4. Crop Rasters & Extract Pixels (MCP Polygon-Level) -------------------
# Congo spatial extent dynamically calculated from MCP vector layers
study_extent_congo <- ext(mcps_congo)
study_extent_congo <- ext(
  xmin(study_extent_congo) - 1.0,
  xmax(study_extent_congo) + 1.0,
  ymin(study_extent_congo) - 1.0,
  ymax(study_extent_congo) + 1.0
)
r_congo_cropped <- terra::crop(r_congo, study_extent_congo)

# Amazon spatial extent dynamically calculated from MCP vector layers
study_extent_amazon <- ext(mcps_amazon)
study_extent_amazon <- ext(
  xmin(study_extent_amazon) - 1.0,
  xmax(study_extent_amazon) + 1.0,
  ymin(study_extent_amazon) - 1.0,
  ymax(study_extent_amazon) + 1.0
)
r_amazon_cropped <- terra::crop(r_amazon, study_extent_amazon)

# Extract pixel values for Congo
cat("Extracting raster pixel distributions inside Congo buffered MCP polygons...\n")
extracted_congo <- terra::extract(r_congo_cropped, mcps_congo, df = TRUE)
mcp_congo_df <- as.data.frame(mcps_congo)
mcp_congo_df$ID <- 1:nrow(mcp_congo_df)
pixel_congo <- merge(extracted_congo, mcp_congo_df, by = "ID") %>%
  filter(!is.na(uoi) & !is.na(frip)) %>%
  select(-ID) %>%
  mutate(basin = "Congo")

# Extract pixel values for Amazon
cat("Extracting raster pixel distributions inside Amazon buffered MCP polygons...\n")
extracted_amazon <- terra::extract(r_amazon_cropped, mcps_amazon, df = TRUE)
mcp_amazon_df <- as.data.frame(mcps_amazon)
mcp_amazon_df$ID <- 1:nrow(mcp_amazon_df)
pixel_amazon <- merge(extracted_amazon, mcp_amazon_df, by = "ID") %>%
  filter(!is.na(uoi) & !is.na(frip)) %>%
  select(-ID) %>%
  mutate(basin = "Amazon")

# Combine pixel datasets
pixel_data <- rbind(pixel_congo, pixel_amazon)
cat(sprintf("✓ Extracted %d Congo and %d Amazon pixels.\n", nrow(pixel_congo), nrow(pixel_amazon)))

# Compute cluster-level polygon means, standard deviations, and pixel counts
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
  )

# Apply minimum trap-day threshold: exclude under-sampled clusters (< 10 trap-days)
# Below ~10 days, community biomass estimates are ecologically unreliable and
# log10 scaling produces steep, noisy weight differences.
n_before <- nrow(joined_data)
joined_data <- joined_data %>% filter(trap_days >= 10)
cat(sprintf("  Trap-day threshold (>= 10): %d → %d clusters retained.\n", n_before, nrow(joined_data)))

# Compute Standard Error of the Mean for each predictor's polygon average:
#   SE = SD / sqrt(N_pixels)
# SE is the statistically appropriate measure of uncertainty in the polygon-wide
# mean, naturally rewarding larger polygons (more pixels) and penalising high
# spatial heterogeneity without the blow-up risk of using raw 1/SD.
joined_data <- joined_data %>%
  mutate(
    uoi_se = uoi_sd / sqrt(n_pixels),
    frip_se = frip_sd / sqrt(n_pixels)
  )

# Regularization constants (median of non-zero SEs) to prevent division by zero
reg_uoi <- median(joined_data$uoi_se[joined_data$uoi_se > 0])
reg_frip <- median(joined_data$frip_se[joined_data$frip_se > 0])

# Compute dual weights using log-effort and Standard Error:
#   w = log10(trap_days) / (SE + reg)
# Log10 captures diminishing returns of sampling effort (100 days = 2x weight
# of 10 days, not 10x). SE in the denominator properly accounts for both
# spatial heterogeneity and polygon size (number of raster pixels).
joined_data <- joined_data %>%
  mutate(
    w_uoi = log10(trap_days) / (uoi_se + reg_uoi),
    w_frip = log10(trap_days) / (frip_se + reg_frip)
  )

# Save the polygon-joined data for records or user review
write_csv(joined_data, "outputs/congo_camera_trap_gee_joined.csv")
cat("✓ Created outputs/congo_camera_trap_gee_joined.csv using polygon extractions.\n")

# --- 5. Extract Regional Sample (1000 pixels) for Panel A --------------------
cat("Extracting a random sample of 1000 regional pixels for UOI vs. FRIP (Panel A)...\n")
set.seed(42) # For scientific reproducibility
all_congo_pixels <- as.data.frame(r_congo_cropped[[c("uoi", "frip")]], na.rm = TRUE) %>% mutate(basin = "Congo")
all_amazon_pixels <- as.data.frame(r_amazon_cropped[[c("uoi", "frip")]], na.rm = TRUE) %>% mutate(basin = "Amazon")

# Sample 500 from each to get a representative 1000 pixel sample
sample_congo <- all_congo_pixels[sample(1:nrow(all_congo_pixels), 500), ]
sample_amazon <- all_amazon_pixels[sample(1:nrow(all_amazon_pixels), 500), ]
pixel_sample <- rbind(sample_congo, sample_amazon)

# --- 6. Generate Figure 1: Spatial Map Comparison (4 panels) ----------------
cat("\nGenerating Figure 1: Spatial Maps with Buffered MCPs (Congo & Amazon)...\n")

# Panel A: Congo GEDI UOI Map
p_map_a <- ggplot() +
  geom_spatraster(data = r_congo_cropped[["uoi"]]) +
  scale_fill_viridis_c(
    option = "plasma", 
    name = "GEDI UOI\n(Congo)",
    na.value = "transparent"
  ) +
  geom_spatvector(
    data = mcps_congo,
    fill = NA, color = "black", linewidth = 0.5, alpha = 0.8
  ) +
  geom_point(
    data = ct_congo,
    aes(x = longitude, y = latitude),
    color = "#FF1744", size = 0.4, alpha = 0.7
  ) +
  theme_pnas(base_size = 8) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6.0, face = "bold"),
    legend.text = element_text(size = 5.5),
    plot.title = element_text(size = 9, face = "bold"),
    axis.title = element_blank()
  ) +
  labs(
    title = "A. Congo Basin Understory Openness Index (UOI)",
    subtitle = "Buffered MCPs & deployments overlaid on UOI"
  )

# Panel B: Congo FRIP Map
p_map_b <- ggplot() +
  geom_spatraster(data = r_congo_cropped[["frip"]]) +
  scale_fill_viridis_c(
    option = "viridis", 
    name = "FRIP\n(Congo)",
    na.value = "transparent"
  ) +
  geom_spatvector(
    data = mcps_congo,
    fill = NA, color = "black", linewidth = 0.5, alpha = 0.8
  ) +
  geom_point(
    data = ct_congo,
    aes(x = longitude, y = latitude),
    color = "#FF1744", size = 0.4, alpha = 0.7
  ) +
  theme_pnas(base_size = 8) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6.0, face = "bold"),
    legend.text = element_text(size = 5.5),
    plot.title = element_text(size = 9, face = "bold"),
    axis.title = element_blank()
  ) +
  labs(
    title = "B. Congo Basin Flooding Role in Productivity (FRIP)",
    subtitle = "Buffered MCPs & deployments overlaid on FRIP"
  )

# Panel C: Amazon GEDI UOI Map
p_map_c <- ggplot() +
  geom_spatraster(data = r_amazon_cropped[["uoi"]]) +
  scale_fill_viridis_c(
    option = "plasma", 
    name = "GEDI UOI\n(Amazon)",
    na.value = "transparent"
  ) +
  geom_spatvector(
    data = mcps_amazon,
    fill = NA, color = "black", linewidth = 0.5, alpha = 0.8
  ) +
  geom_point(
    data = ct_amazon,
    aes(x = longitude, y = latitude),
    color = "#FF1744", size = 0.4, alpha = 0.7
  ) +
  theme_pnas(base_size = 8) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6.0, face = "bold"),
    legend.text = element_text(size = 5.5),
    plot.title = element_text(size = 9, face = "bold"),
    axis.title = element_blank()
  ) +
  labs(
    title = "C. Amazon Basin Understory Openness Index (UOI)",
    subtitle = "Buffered MCPs & deployments overlaid on UOI"
  )

# Panel D: Amazon FRIP Map
p_map_d <- ggplot() +
  geom_spatraster(data = r_amazon_cropped[["frip"]]) +
  scale_fill_viridis_c(
    option = "viridis", 
    name = "FRIP\n(Amazon)",
    na.value = "transparent"
  ) +
  geom_spatvector(
    data = mcps_amazon,
    fill = NA, color = "black", linewidth = 0.5, alpha = 0.8
  ) +
  geom_point(
    data = ct_amazon,
    aes(x = longitude, y = latitude),
    color = "#FF1744", size = 0.4, alpha = 0.7
  ) +
  theme_pnas(base_size = 8) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6.0, face = "bold"),
    legend.text = element_text(size = 5.5),
    plot.title = element_text(size = 9, face = "bold"),
    axis.title = element_blank()
  ) +
  labs(
    title = "D. Amazon Basin Flooding Role in Productivity (FRIP)",
    subtitle = "Buffered MCPs & deployments overlaid on FRIP"
  )

# Combine using cowplot (Double column layout)
fig_spatial <- cowplot::plot_grid(
  p_map_a, p_map_b,
  p_map_c, p_map_d,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

# Save spatial map comparison
save_pnas(
  plot = fig_spatial, 
  filename = "outputs/congo_ct_gee_spatial_maps.png", 
  type = "double", 
  height_cm = 16.0
)
cat("✓ Saved Figure 1: outputs/congo_ct_gee_spatial_maps.png\n")

# --- 7. Generate Figure 2: Spatial Heterogeneity Box Plots (2 panels) --------
cat("Generating Figure 2: Pixel Distribution Box Plots...\n")

# Order the factor levels of cluster_id by Large Fauna Biomass Index (B_H_gt50) ascending
pixel_data <- pixel_data %>%
  mutate(cluster_id = reorder(factor(cluster_id), B_H_gt50))

# Panel A: Box plots of GEDI UOI pixels in each robust cluster, faceted by Basin
p_box_a <- ggplot(pixel_data, aes(x = cluster_id, y = uoi, fill = B_H_gt50)) +
  geom_boxplot(color = "black", outlier.size = 0.3, outlier.alpha = 0.4, linewidth = 0.3) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "Large Fauna Biomass\nIndex (B_H_gt50)"
  ) +
  facet_wrap(~basin, scales = "free_x", ncol = 2) +
  theme_pnas(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 4.5),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    plot.title = element_text(size = 9, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 6.0, face = "bold"),
    legend.text = element_text(size = 5.5)
  ) +
  labs(
    title = "A. GEDI Understory Openness Index (UOI) Pixel Distribution",
    x = "Camera Trap Cluster (ordered by Large Fauna Biomass, B_H_gt50)",
    y = "GEDI UOI Pixel Values"
  )

# Panel B: Box plots of GEE FRIP pixels in each robust cluster, faceted by Basin
p_box_b <- ggplot(pixel_data, aes(x = cluster_id, y = frip, fill = B_H_gt50)) +
  geom_boxplot(color = "black", outlier.size = 0.3, outlier.alpha = 0.4, linewidth = 0.3) +
  scale_fill_viridis_c(
    option = "viridis",
    name = "Large Fauna Biomass\nIndex (B_H_gt50)"
  ) +
  facet_wrap(~basin, scales = "free_x", ncol = 2) +
  theme_pnas(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 4.5),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    plot.title = element_text(size = 9, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 6.0, face = "bold"),
    legend.text = element_text(size = 5.5)
  ) +
  labs(
    title = "B. Flooding Role in Productivity (FRIP) Pixel Distribution",
    x = "Camera Trap Cluster (ordered by Large Fauna Biomass, B_H_gt50)",
    y = "Flooding Role in Productivity (FRIP)"
  )

# Combine using cowplot (Double column layout)
fig_boxplots <- cowplot::plot_grid(
  p_box_a, p_box_b,
  ncol = 1,
  align = "v",
  axis = "lr"
)

# Save the box plots
save_pnas(
  plot = fig_boxplots, 
  filename = "outputs/congo_ct_gee_spatial_distributions.png", 
  type = "double", 
  height_cm = 16.0
)
cat("✓ Saved Figure 2 (Spatial Distributions): outputs/congo_ct_gee_spatial_distributions.png\n")

# --- 8. Generate Figure 3: Focused Statistical Relationships (5 panels) ------
cat("Generating Figure 3: Scatter Plots (Effort-Weighted)...\n")

# ── Panel A: Regional UOI vs. FRIP (Unweighted, 1000 random pixels) ────────
fit_a <- lm(frip ~ uoi, data = pixel_sample)
summary_a <- summary(fit_a)
r_val_a <- cor(pixel_sample$uoi, pixel_sample$frip)
p_val_a <- summary_a$coefficients[2, 4]
p_text_a <- if (p_val_a < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val_a)
stat_label_a <- sprintf("Pearson r = %.2f\n%s (N = 1000)", r_val_a, p_text_a)

p_stat_a <- ggplot(pixel_sample, aes(x = uoi, y = frip)) +
  geom_smooth(
    method = "lm", 
    color = "#6A1B9A", 
    fill = "#E1BEE7", 
    alpha = 0.2, 
    linewidth = 0.6,
    na.rm = TRUE
  ) +
  geom_point(
    aes(fill = basin), 
    color = "black", 
    shape = 21, 
    stroke = 0.15,
    size = 0.8,
    alpha = 0.4
  ) +
  scale_fill_manual(
    values = c("Congo" = "#1E88E5", "Amazon" = "#D81B60"),
    name = "Basin"
  ) +
  annotate(
    "text", 
    x = Inf, y = Inf, 
    label = stat_label_a, 
    hjust = 1.1, vjust = 1.2, 
    family = base_family,
    size = 2.4, 
    fontface = "bold",
    color = "grey20"
  ) +
  theme_pnas(base_size = 8) +
  labs(
    title = "A. Regional Understory Openness vs. Flooding Role in Productivity (FRIP)",
    x = "GEDI Understory Openness Index (UOI)",
    y = "Flooding Role in Productivity (FRIP)"
  )

# ── Helper for weighted relation plots (Panels B, C, D, E) ──────────────────
# Weights use log10(trap_days) / (SE + reg) — see weighting block above.
make_weighted_relation_plot <- function(df, x_col, y_col, se_col, weight_col, x_label, y_label, title, fill_scale, line_color) {
  # Calculate weighted linear model
  w_fit <- lm(as.formula(sprintf("%s ~ %s", y_col, x_col)), weights = df[[weight_col]], data = df)
  summary_fit <- summary(w_fit)
  
  # Weighted Pearson correlation r = sqrt(R2) * sign(slope)
  r_val <- sqrt(summary_fit$r.squared) * sign(coef(w_fit)[2])
  
  # P-value for the weighted slope
  p_val <- summary_fit$coefficients[2, 4]
  
  # Format p-value elegantly
  p_text <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)
  stat_label <- sprintf("Log-weighted r = %.2f\n%s (N = %d)", r_val, p_text, nrow(df))
  
  ggplot(df, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    geom_smooth(
      method = "lm", 
      aes(weight = !!sym(weight_col)),
      color = line_color, 
      fill = line_color, 
      alpha = 0.12, 
      linewidth = 0.6,
      na.rm = TRUE
    ) +
    geom_point(
      aes(size = trap_days, fill = !!sym(se_col), shape = basin), 
      color = "black", 
      stroke = 0.3,
      alpha = 0.85
    ) +
    scale_shape_manual(
      values = c("Congo" = 21, "Amazon" = 24),
      name = "Basin"
    ) +
    annotate(
      "text", 
      x = Inf, y = Inf, 
      label = stat_label, 
      hjust = 1.1, vjust = 1.2, 
      family = base_family,
      size = 2.4, 
      fontface = "bold",
      color = "grey20"
    ) +
    scale_size_continuous(
      range = c(1.5, 6.0), 
      guide = "none" # Hide size legend to keep plots compact
    ) +
    fill_scale +
    theme_pnas(base_size = 8) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 5.5, face = "bold"),
      legend.text = element_text(size = 5.0),
      legend.key.height = unit(0.25, "cm"),
      legend.key.width = unit(0.12, "cm"),
      legend.margin = margin(l = -2, r = 0, t = 0, b = 0, unit = "pt")
    ) +
    labs(
      title = title,
      x = x_label,
      y = y_label
    )
}

# ── Generate 4 individual weighted panels ───────────────────────────────────
# Panel B: GEDI UOI vs. Large Fauna Biomass Index (B_H_gt50)
p_stat_b <- make_weighted_relation_plot(
  df = joined_data,
  x_col = "uoi",
  y_col = "B_H_gt50",
  se_col = "uoi_se",
  weight_col = "w_uoi",
  x_label = "GEDI Understory Openness Index (UOI)",
  y_label = "Large Fauna Biomass Index (B_H_gt50)",
  title = "B. Large Fauna Biomass vs. Understory Openness",
  fill_scale = scale_fill_viridis_c(option = "magma", name = "UOI Spatial SE"),
  line_color = "#4527A0"
)

# Panel C: GEE FRIP vs. Large Fauna Biomass Index (B_H_gt50)
p_stat_c <- make_weighted_relation_plot(
  df = joined_data,
  x_col = "frip",
  y_col = "B_H_gt50",
  se_col = "frip_se",
  weight_col = "w_frip",
  x_label = "Flooding Role in Productivity (FRIP)",
  y_label = "Large Fauna Biomass Index (B_H_gt50)",
  title = "C. Large Fauna Biomass vs. Flooding Role in Productivity (FRIP)",
  fill_scale = scale_fill_viridis_c(option = "mako", name = "FRIP Spatial SE"),
  line_color = "#C62828"
)

# Panel D: GEDI UOI vs. Megafauna Biomass Index (B_H_gt100)
p_stat_d <- make_weighted_relation_plot(
  df = joined_data,
  x_col = "uoi",
  y_col = "B_H_gt100",
  se_col = "uoi_se",
  weight_col = "w_uoi",
  x_label = "GEDI Understory Openness Index (UOI)",
  y_label = "Megafauna Biomass Index (B_H_gt100)",
  title = "D. Megafauna Biomass vs. Understory Openness",
  fill_scale = scale_fill_viridis_c(option = "magma", name = "UOI Spatial SE"),
  line_color = "#1565C0"
)

# Panel E: GEE FRIP vs. Megafauna Biomass Index (B_H_gt100)
p_stat_e <- make_weighted_relation_plot(
  df = joined_data,
  x_col = "frip",
  y_col = "B_H_gt100",
  se_col = "frip_se",
  weight_col = "w_frip",
  x_label = "Flooding Role in Productivity (FRIP)",
  y_label = "Megafauna Biomass Index (B_H_gt100)",
  title = "E. Megafauna Biomass vs. Flooding Role in Productivity (FRIP)",
  fill_scale = scale_fill_viridis_c(option = "mako", name = "FRIP Spatial SE"),
  line_color = "#2E7D32"
)

# ── Assemble final layout using cowplot ─────────────────────────────────────
# Combine bottom 4 panels into a 2x2 grid
fig_bottom <- cowplot::plot_grid(
  p_stat_b, p_stat_c, p_stat_d, p_stat_e,
  ncol = 2,
  align = "vh",
  labels = NULL
)

# Combine Panel A (spanning full top row) and the bottom 2x2 grid
fig_stats <- cowplot::plot_grid(
  p_stat_a, fig_bottom,
  ncol = 1,
  rel_heights = c(1.0, 2.0), # Panel A spans top, bottom is 2x2
  align = "v"
)

# Save scatter plot relationship matrix
save_pnas(
  plot = fig_stats,
  filename = "outputs/congo_ct_gee_statistical_plots.png",
  type = "double",
  height_cm = 18.0
)
cat("✓ Saved Figure 3: outputs/congo_ct_gee_statistical_plots.png\n")

# ── Generate Figure 4: Focused Biomass & UOI Relationships (4 panels, 2x2) ──
cat("\nGenerating Figure 4: Biomass & UOI Relationships (2x2)...\n")

p_biomass_a <- make_weighted_relation_plot(
  df = joined_data,
  x_col = "uoi",
  y_col = "B_H_index",
  se_col = "uoi_se",
  weight_col = "w_uoi",
  x_label = "GEDI Understory Openness Index (UOI)",
  y_label = "Total Biomass Index (B_H_index)",
  title = "A. Total Biomass vs. Understory Openness",
  fill_scale = scale_fill_viridis_c(option = "magma", name = "UOI Spatial SE"),
  line_color = "#311B92"
)

p_biomass_b <- make_weighted_relation_plot(
  df = joined_data,
  x_col = "uoi",
  y_col = "B_H_gt50",
  se_col = "uoi_se",
  weight_col = "w_uoi",
  x_label = "GEDI Understory Openness Index (UOI)",
  y_label = "Large Fauna Biomass Index (>50kg, B_H_gt50)",
  title = "B. Large Fauna Biomass vs. Understory Openness",
  fill_scale = scale_fill_viridis_c(option = "magma", name = "UOI Spatial SE"),
  line_color = "#4527A0"
)

p_biomass_c <- make_weighted_relation_plot(
  df = joined_data,
  x_col = "uoi",
  y_col = "B_H_gt100",
  se_col = "uoi_se",
  weight_col = "w_uoi",
  x_label = "GEDI Understory Openness Index (UOI)",
  y_label = "Megafauna Biomass Index (>100kg, B_H_gt100)",
  title = "C. Megafauna Biomass vs. Understory Openness",
  fill_scale = scale_fill_viridis_c(option = "magma", name = "UOI Spatial SE"),
  line_color = "#1565C0"
)

p_biomass_d <- make_weighted_relation_plot(
  df = joined_data,
  x_col = "uoi",
  y_col = "megafauna_fraction",
  se_col = "uoi_se",
  weight_col = "w_uoi",
  x_label = "GEDI Understory Openness Index (UOI)",
  y_label = "Megafauna Biomass Fraction (%)",
  title = "D. Megafauna Fraction vs. Understory Openness",
  fill_scale = scale_fill_viridis_c(option = "magma", name = "UOI Spatial SE"),
  line_color = "#2E7D32"
)

fig_biomass_uoi <- cowplot::plot_grid(
  p_biomass_a, p_biomass_b, p_biomass_c, p_biomass_d,
  ncol = 2,
  align = "vh",
  labels = NULL
)

save_pnas(
  plot = fig_biomass_uoi,
  filename = "outputs/congo_ct_gee_biomass_uoi_relationships.png",
  type = "double",
  height_cm = 16.0
)
cat("✓ Saved Figure 4: outputs/congo_ct_gee_biomass_uoi_relationships.png\n")

cat("\n=== Congo Camera Trap & GEE Analysis Completed Successfully! ===\n")
