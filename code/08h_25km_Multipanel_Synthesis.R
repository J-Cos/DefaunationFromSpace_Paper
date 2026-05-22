# =============================================================================
# 08h_25km_Multipanel_Synthesis.R
#
# Generates a premium, publication-quality multipanel figure for the 25km scale:
#   - Panel A: Scatterplot of empirical data, fitted weight-normalized Tweedie GLM,
#              and extrapolation up to UOI = 1.0.
#   - Panel B: Boxplots of GEDI UOI pixels inside each buffered MCP polygon with
#              individual pixel jitter overlaid, ordered by cluster biomass index
#              and colored by continent.
#   - Panel C: Spatial maps for each continent:
#              - GEDI UOI (25km)
#              - Predicted Biomass Index (25km) from Tweedie GLM
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

# --- 4. Fit Weight-Normalized Tweedie GLM (25 km) ----------------------------
tw_model <- gam(
  B_H_index ~ uoi,
  family = tw(),
  weights = w_uoi_norm,
  data = joined_data,
  method = "REML"
)

s_tw <- summary(tw_model)
p_val <- s_tw$p.pv["uoi"]
dev_expl <- s_tw$dev.expl * 100
p_text <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.4f", p_val)

cat(sprintf("✓ Tweedie GLM (25 km) fitted: DevExpl = %.2f%%, %s\n", dev_expl, p_text))

# --- 5. PANEL A: Scatter Plot and Model Extrapolation ------------------------
cat("  Generating Panel A (Scatterplot with extrapolation to UOI=1.0)...\n")

# Extend the sequence up to UOI = 1.00
uoi_seq <- seq(0.91, 1.00, length.out = 300)
pred_df <- data.frame(uoi = uoi_seq)
pred <- predict(tw_model, newdata = pred_df, type = "link", se.fit = TRUE)

pred_df$fit <- exp(pred$fit)
pred_df$lower <- exp(pred$fit - 1.96 * pred$se.fit)
pred_df$upper <- exp(pred$fit + 1.96 * pred$se.fit)

stat_label <- sprintf(
  "Tweedie GLM (25 km)\nDeviance expl. = %.1f%%\n%s\nN = %d clusters",
  dev_expl, p_text, nrow(joined_data)
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
    title = "A. Tweedie Model Fit & Extrapolation (25 km)",
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

pred_congo <- predict(tw_model, newdata = congo_cells, type = "link", se.fit = TRUE)
congo_cells$pred <- exp(pred_congo$fit)

r_pred_congo <- rast(r_congo_cropped[["uoi"]])
names(r_pred_congo) <- "pred"
values(r_pred_congo) <- NA
r_pred_congo[congo_cells$cell] <- as.vector(congo_cells$pred)

# --- 7.2 Amazon Predictions and Mapping
cat("    Processing Amazon spatial predictions...\n")
r_amazon_cropped <- crop(r_amazon, study_extent_amazon)
amazon_cells <- as.data.frame(r_amazon_cropped[["uoi"]], cells = TRUE, xy = TRUE, na.rm = TRUE)
names(amazon_cells)[names(amazon_cells) == "uoi"] <- "uoi"

pred_amazon <- predict(tw_model, newdata = amazon_cells, type = "link", se.fit = TRUE)
amazon_cells$pred <- exp(pred_amazon$fit)

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

# Generate Map Sub-panels (UOI using viridis, Predictions using inferno clamped at 5000)
p_congo_uoi  <- make_map_panel(r_congo_cropped[["uoi"]], mcps_congo, "viridis", "C.1 Congo GEDI UOI (25 km)", "UOI", limits = c(0.91, 1.00))
p_congo_pred <- make_map_panel(r_pred_congo,             mcps_congo, "inferno", "C.2 Congo Predicted Biomass (25 km)", "Biomass", limits = c(0, 5000), winsorize = TRUE)
p_amazon_uoi  <- make_map_panel(r_amazon_cropped[["uoi"]], mcps_amazon, "viridis", "C.3 Amazon GEDI UOI (25 km)", "UOI", limits = c(0.91, 1.00))
p_amazon_pred <- make_map_panel(r_pred_amazon,             mcps_amazon, "inferno", "C.4 Amazon Predicted Biomass (25 km)", "Biomass", limits = c(0, 5000), winsorize = TRUE)

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
