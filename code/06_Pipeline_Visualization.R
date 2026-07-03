# =============================================================================
# code/06_Pipeline_Visualization.R
#
# Generates two premium 6-panel PNAS-style publication-grade figures summarizing
# the GEDI raw data pipeline (Figure 1) and the Camera Trap data pipeline (Figure 2).
#
# Figure 1: GEDI Raw Data Pipeline & Shot-Sparsity Distribution (6 Panels)
#   (A) Congo Basin regional map of GEDI UOI (5km) with cluster points
#   (B) Amazon Basin regional map of GEDI UOI (5km) with cluster points
#   (C) Zoomed-in inset at raw 500m scale showing GEDI shot density sparsity (gedi_n)
#   (D) Zoomed-in inset at raw 500m scale showing GEDI understory openness (uoi)
#   (E) Shot-noise uncertainty reduction: GEDI UOI Standard Error vs. Shot Count
#   (F) UOI Distribution Scaling Law: raw 500m vs aggregated 5km vs aggregated 20km
#
# Figure 2: Camera Trap Wildlife Ingestion, Clustering, and Temporal Calibration (6 Panels)
#   (A) Species rank-abundance curves for both basins
#   (B) Spatial clustering map of Congo deployments colored by mathematical cluster ID
#   (C) Spatial clustering map of Amazon deployments colored by mathematical cluster ID
#   (D) Temporal calibration timeline showing cluster survey spans vs. GEDI mission launch
#   (E) Proportional horizontal bar plot of taxonomic order composition
#   (F) Vertebrate body mass probability density distribution (Congo vs. Amazon)
# =============================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(cowplot)
library(patchwork)
library(dplyr)
library(readr)
library(sf)
library(scales)
library(tidyr)
library(jsonlite)

config_path <- if (file.exists("code/config.json")) "code/config.json" else "config.json"
config <- jsonlite::read_json(config_path)
brackets <- unlist(config$temporal_decay$brackets)
weights <- unlist(config$temporal_decay$weights)

cat("============================================================\n")
cat("=== Generating 6-Panel Manuscript Figures 1 and 2 (Pipeline Summary) ===\n")
cat("============================================================\n\n")

# --- Source standard styling and helpers -------------------------------------
source("code/functions/theme_pnas.R")
source("code/functions/calibration_helpers.R")

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# Color palettes (defined in theme_pnas.R: pal_region, pal_taxon, etc.)

# =============================================================================
# FIGURE 1: GEDI RAW DATA PIPELINE & SHOT-SPARSITY DISTRIBUTION (6 PANELS)
# =============================================================================
cat("--- Drafting Figure 1: GEDI Raw Data Pipeline (6 Panels) ---\n")

# Helper to dynamically assign names based on layer counts
assign_multiscale_names <- function(r) {
  num_layers <- nlyr(r)
  if (num_layers == 12) {
    names(r) <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                  "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  } else if (num_layers == 11) {
    names(r) <- c("frip", "frip_mk_tau", "uoi", "rh98", "gedi_n",
                  "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  }
  return(r)
}

assign_native_names <- function(r) {
  num_layers <- nlyr(r)
  if (num_layers == 11) {
    names(r) <- c("uoi", "uoi_sd", "rh98", "gedi_n", "elevation", "slope", "hnd",
                  "precip", "clay", "forest_fraction", "Npp_median")
  } else if (num_layers == 10) {
    names(r) <- c("uoi", "rh98", "gedi_n", "elevation", "slope", "hnd",
                  "precip", "clay", "forest_fraction", "Npp_median")
  }
  return(r)
}

# Load 5km aggregated maps for all three basins
r_congo_5km  <- assign_multiscale_names(rast("outputs/EOdata/analysis_stack_5000_Congo.tif"))
r_amazon_5km <- assign_multiscale_names(rast("outputs/EOdata/analysis_stack_5000_Amazon.tif"))
r_seasia_5km <- assign_multiscale_names(rast("outputs/EOdata/analysis_stack_5000_SE_Asia.tif"))

# Load country outlines and crop for background
countries_v <- vect("data/world-administrative-boundaries")

# Load camera trap cluster coordinates for overlay
mcps <- vect("outputs/camera_traps_robust_buffered_mcps.geojson")
mcps_congo  <- mcps[mcps$region == "Congo", ]
mcps_amazon  <- mcps[mcps$region == "Amazon", ]
mcps_seasia  <- mcps[mcps$region == "SE_Asia", ]

# Set up study bounding boxes: EXACTLY 50° Wide × 30° High (-15 to 15 Lat) to match Figure 5 and S1
ext_a <- ext(-85, -35, -15, 15)
ext_c <- ext(-5, 45, -15, 15)
ext_s <- ext(90, 140, -15, 15)

r_c_crop <- crop(r_congo_5km, ext_c)
r_a_crop <- crop(r_amazon_5km, ext_a)
r_s_crop <- crop(r_seasia_5km, ext_s)

# Scale GEDI shot density to raw shot count (number of GEDI footprints in 5km cell)
r_c_crop$gedi_n <- r_c_crop$gedi_n * 40000
r_a_crop$gedi_n <- r_a_crop$gedi_n * 40000
r_s_crop$gedi_n <- r_s_crop$gedi_n * 40000

countries_c <- crop(countries_v, ext_c)
countries_a <- crop(countries_v, ext_a)
countries_s <- crop(countries_v, ext_s)

# --- UOI Fill Scale ---
uoi_fill_scale <- scale_fill_viridis_c(
  option = "viridis",
  name = "GEDI Understory Openness (UOI)",
  limits = c(0.92, 0.975),
  oob = scales::squish,
  na.value = "transparent"
)

# --- Shot Count Fill Scale ---
n_fill_scale <- scale_fill_viridis_c(
  option = "mako",
  name = "GEDI Shot Count",
  trans = "log10",
  limits = c(10, 6000),
  breaks = c(10, 100, 1000, 6000),
  labels = c("10", "100", "1K", "6K"),
  oob = scales::squish,
  na.value = "transparent"
)

t_theme <- theme_pnas(base_size = 7.5)

# Define premium map theme for Figure 1
map_theme_fig1 <- theme_pnas(base_size = 7.5) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 5.0),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "white", linewidth = 0.2), # Soft white grid lines
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#EBF5FB", color = NA), # Soft blue oceans
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(4, 4, 4, 4, "pt"),
    panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5)
  )

# --- Figure 1 Panels A, C, E: UOI Maps ---
p1_a <- ggplot() +
  geom_spatvector(data = countries_a, fill = "white", colour = NA) +
  geom_spatraster(data = r_a_crop, aes(fill = uoi)) +
  uoi_fill_scale +
  geom_spatvector(data = countries_a, fill = NA, colour = "grey80", linewidth = 0.2) +
  coord_sf(xlim = c(-85, -35), ylim = c(-15, 15), expand = FALSE) +
  map_theme_fig1 +
  labs(title = "A")

p1_b <- ggplot() +
  geom_spatvector(data = countries_c, fill = "white", colour = NA) +
  geom_spatraster(data = r_c_crop, aes(fill = uoi)) +
  uoi_fill_scale +
  geom_spatvector(data = countries_c, fill = NA, colour = "grey80", linewidth = 0.2) +
  coord_sf(xlim = c(-5, 45), ylim = c(-15, 15), expand = FALSE) +
  map_theme_fig1 +
  labs(title = "C")

p1_c <- ggplot() +
  geom_spatvector(data = countries_s, fill = "white", colour = NA) +
  geom_spatraster(data = r_s_crop, aes(fill = uoi)) +
  uoi_fill_scale +
  geom_spatvector(data = countries_s, fill = NA, colour = "grey80", linewidth = 0.2) +
  coord_sf(xlim = c(90, 140), ylim = c(-15, 15), expand = FALSE) +
  map_theme_fig1 +
  labs(title = "E")

# --- Figure 1 Panels B, D, F: Shot Count Maps ---
p1_d <- ggplot() +
  geom_spatvector(data = countries_a, fill = "white", colour = NA) +
  geom_spatraster(data = r_a_crop, aes(fill = gedi_n)) +
  n_fill_scale +
  geom_spatvector(data = countries_a, fill = NA, colour = "grey80", linewidth = 0.2) +
  coord_sf(xlim = c(-85, -35), ylim = c(-15, 15), expand = FALSE) +
  map_theme_fig1 +
  labs(title = "B")

p1_e <- ggplot() +
  geom_spatvector(data = countries_c, fill = "white", colour = NA) +
  geom_spatraster(data = r_c_crop, aes(fill = gedi_n)) +
  n_fill_scale +
  geom_spatvector(data = countries_c, fill = NA, colour = "grey80", linewidth = 0.2) +
  coord_sf(xlim = c(-5, 45), ylim = c(-15, 15), expand = FALSE) +
  map_theme_fig1 +
  labs(title = "D")

p1_f <- ggplot() +
  geom_spatvector(data = countries_s, fill = "white", colour = NA) +
  geom_spatraster(data = r_s_crop, aes(fill = gedi_n)) +
  n_fill_scale +
  geom_spatvector(data = countries_s, fill = NA, colour = "grey80", linewidth = 0.2) +
  coord_sf(xlim = c(90, 140), ylim = c(-15, 15), expand = FALSE) +
  map_theme_fig1 +
  labs(title = "F")

# --- Figure 1 Panel G: Uncertainty Decay vs. GEDI Shot Density ---
df_c_pts <- as.data.frame(r_c_crop[[c("uoi_sd", "gedi_n")]], na.rm = TRUE)
df_c_pts$basin <- "Congo"
df_a_pts <- as.data.frame(r_a_crop[[c("uoi_sd", "gedi_n")]], na.rm = TRUE)
df_a_pts$basin <- "Amazon"
df_s_pts <- as.data.frame(r_s_crop[[c("uoi_sd", "gedi_n")]], na.rm = TRUE)
df_s_pts$basin <- "SE_Asia"
df_pts <- rbind(df_c_pts, df_a_pts, df_s_pts)

# Downsample for clean plotting density
set.seed(42)
df_pts <- df_pts %>%
  group_by(basin) %>%
  slice_sample(n = 2000) %>%
  ungroup()

p1_g <- ggplot(df_pts, aes(x = gedi_n, y = uoi_sd, color = basin)) +
  geom_point(alpha = 0.25, size = 0.5, stroke = 0) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE, linewidth = 0.75) +
  scale_color_manual(values = pal_region, labels = c("Amazon" = "Amazon", "Congo" = "Congo", "SE_Asia" = "SE Asia"), name = "Region") +
  scale_x_continuous(trans = "log10", labels = comma_format()) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), limits = c(0, 0.008), oob = scales::squish) +
  labs(
    title = "G",
    x = "GEDI Shot Count (log scale)",
    y = "GEDI UOI Standard Error (SE)"
  ) +
  t_theme +
  theme(
    legend.position = "none",
    aspect.ratio = 0.6,
    plot.title = element_text(size = 8.5, face = "bold", hjust = 0),
    panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5),
    plot.margin = margin(4, 4, 4, 4, "pt")
  )

# --- Figure 1 Panel H: UOI Distribution divided by Continent ---
df_c_uoi <- data.frame(uoi = as.data.frame(r_c_crop[["uoi"]], na.rm = TRUE)$uoi, basin = "Congo")
df_a_uoi <- data.frame(uoi = as.data.frame(r_a_crop[["uoi"]], na.rm = TRUE)$uoi, basin = "Amazon")
df_s_uoi <- data.frame(uoi = as.data.frame(r_s_crop[["uoi"]], na.rm = TRUE)$uoi, basin = "SE_Asia")
df_uoi_dist <- rbind(df_c_uoi, df_a_uoi, df_s_uoi)

p1_h <- ggplot(df_uoi_dist, aes(x = uoi, fill = basin, color = basin)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  scale_fill_manual(values = pal_region, name = "Region") +
  scale_color_manual(values = pal_region, name = "Region") +
  scale_x_continuous(limits = c(0.92, 0.975), breaks = seq(0.92, 0.97, by = 0.01)) +
  labs(
    title = "H",
    x = "GEDI Understory Openness Index (UOI)",
    y = "Probability Density"
  ) +
  t_theme +
  theme(
    legend.position = "none",
    aspect.ratio = 0.6,
    plot.title = element_text(size = 8.5, face = "bold", hjust = 0),
    panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5),
    plot.margin = margin(4, 4, 4, 4, "pt")
  )

# --- Extract legends as grobs ---
get_robust_legend <- function(plot) {
  g <- ggplotGrob(plot)
  grob_names <- sapply(g$grobs, function(x) x$name)
  idx <- grep("guide-box", grob_names)
  if (length(idx) > 0) {
    for (i in idx) {
      if (!inherits(g$grobs[[i]], "zeroGrob")) {
        return(g$grobs[[i]])
      }
    }
  }
  return(NULL)
}

# UOI colorbar legend
p_uoi_leg <- ggplot() +
  geom_spatraster(data = r_c_crop, aes(fill = uoi)) +
  scale_fill_viridis_c(
    option = "viridis",
    name = "GEDI Understory Openness Index (UOI)",
    limits = c(0.92, 0.975), oob = scales::squish, na.value = "transparent",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
      label.position = "bottom", barwidth = unit(5, "cm"), barheight = unit(0.12, "cm"))
  ) + t_theme + theme(legend.position = "bottom",
    legend.title = element_text(size = 5.5, face = "bold"),
    legend.text = element_text(size = 4.5))
uoi_legend <- get_robust_legend(p_uoi_leg)

# Shot Count colorbar legend
p_n_leg <- ggplot() +
  geom_spatraster(data = r_c_crop, aes(fill = gedi_n)) +
  scale_fill_viridis_c(
    option = "mako",
    name = "GEDI Shot Count (5 km cell, log scale)",
    trans = "log10", limits = c(10, 6000),
    breaks = c(10, 100, 1000, 6000), labels = c("10", "100", "1K", "6K"),
    oob = scales::squish, na.value = "transparent",
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
      label.position = "bottom", barwidth = unit(5, "cm"), barheight = unit(0.12, "cm"))
  ) + t_theme + theme(legend.position = "bottom",
    legend.title = element_text(size = 5.5, face = "bold"),
    legend.text = element_text(size = 4.5))
n_legend <- get_robust_legend(p_n_leg)

# Region legend
p_reg_leg <- ggplot(df_pts, aes(x = gedi_n, y = uoi_sd, color = basin)) +
  geom_point() +
  scale_color_manual(values = pal_region,
    labels = c("Amazon" = "Amazon", "Congo" = "Congo", "SE_Asia" = "SE Asia"),
    name = "Region:") +
  theme_pnas(base_size = 7.5) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
    legend.title = element_text(size = 6.0, face = "bold"),
    legend.text = element_text(size = 5.5),
    legend.key.height = unit(0.25, "cm"), legend.key.width = unit(0.3, "cm"))
region_legend <- get_robust_legend(p_reg_leg)

# --- Patchwork assembly ---
# Flat grid with wrap_elements() for legends (these don't affect panel alignment).
# Heights: maps (1), maps (1), maps (1), colorbars (0.12), G/H (1)
fig1_grid <- p1_a + p1_d +
             p1_b + p1_e +
             p1_c + p1_f +
             wrap_elements(full = uoi_legend) + wrap_elements(full = n_legend) +
             p1_g + p1_h +
             plot_layout(ncol = 2, heights = c(1, 1, 1, 0.12, 1))

# Stack the single centered region legend below the grid
fig1_final <- fig1_grid / wrap_elements(full = region_legend) +
  plot_layout(heights = c(1, 0.025))

fig1_png <- "figures/figure1_gedi_pipeline.png"
ggsave(filename = fig1_png, plot = fig1_final, width = 20.0, height = 27.5, units = "cm", dpi = 600, bg = "white")
cat("✓ 8-Panel Figure 1 successfully saved to:", fig1_png, "\n\n")

# Save PDF version
fig1_pdf <- sub("\\.png$", ".pdf", fig1_png)
ggsave(filename = fig1_pdf, plot = fig1_final, width = 20.0, height = 27.5, units = "cm", dpi = 600, bg = "white")


# =============================================================================
# FIGURE 2: CAMERA TRAP WILDLIFE INGESTION, CLUSTERING, AND TEMPORAL CALIBRATION
# =============================================================================
cat("--- Drafting Figure 2: Camera Trap Wildlife Ingestion & Calibration (7 Panels) ---\n")

# Load robust, GEDI-overlapping detections with pre-assigned cluster_id
det_path <- "outputs/camera_traps_robust_detections.csv"
if (!file.exists(det_path)) {
  stop("Camera trap robust detections missing: ", det_path)
}
det_all <- read_csv(det_path, show_col_types = FALSE)

# --- Panel A: Individual size distribution curve ---
bin_width <- 0.25
df_binned_mass <- det_all %>%
  filter(!is.na(body_mass_kg) & body_mass_kg > 0) %>%
  mutate(log_mass = log10(body_mass_kg)) %>%
  mutate(bin_center = round(log_mass / bin_width) * bin_width) %>%
  group_by(region, bin_center) %>%
  summarise(n_ind = sum(n_detections), .groups = "drop") %>%
  mutate(body_mass_kg = 10^bin_center)

p2_a <- ggplot(df_binned_mass, aes(x = body_mass_kg, y = n_ind, color = region)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2, alpha = 0.8) +
  scale_color_manual(values = pal_region, labels = c("Amazon" = "Amazon", "Congo" = "Congo", "SE_Asia" = "SE Asia"), name = NULL) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  labs(
    title = "D",
    x = "Mammal Body Mass (kg, log scale)",
    y = "Total Detections (log scale)"
  ) +
  t_theme +
  theme(
    aspect.ratio = 0.45,
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.text = element_text(size = 5.5),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.key = element_blank(),
    legend.key.size = unit(0.2, "cm"),
    plot.margin = margin(4, 4, 4, 4, "pt")
  )

get_clustered_deployments <- function(basin_name) {
  unique_deps <- det_all %>%
    filter(region == basin_name) %>%
    select(longitude, latitude, deployment_id, trap_days, cluster_id) %>%
    distinct()
  
  # Calculate centroids and sum trap days per cluster using pre-assigned cluster_id
  centroids <- unique_deps %>%
    group_by(cluster_id) %>%
    summarise(
      lon = mean(longitude),
      lat = mean(latitude),
      total_trap_days = sum(trap_days, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(list(deps = unique_deps, centroids = centroids))
}

# Load pre-clustered deployments for each basin
congo_c_info <- get_clustered_deployments("Congo")
amazon_c_info <- get_clustered_deployments("Amazon")
seasia_c_info <- get_clustered_deployments("SE_Asia")

# Define study bounding boxes: Exactly 50° Wide × 30° High (-15 to 15 Lat)
ext_amazon_map <- ext(-85, -35, -15, 15)
ext_congo_map  <- ext(-5, 45, -15, 15)
ext_seasia_map <- ext(90, 140, -15, 15)

# Crop country boundaries to each region's bounding box
countries_amazon_map <- crop(countries_v, ext_amazon_map)
countries_congo_map  <- crop(countries_v, ext_congo_map)
countries_seasia_map <- crop(countries_v, ext_seasia_map)

# Define clean PNAS-style map theme with soft light blue ocean fill
map_p_theme <- t_theme +
  theme(
    panel.background = element_rect(fill = "#EBF5FB", color = NA), # Soft light blue ocean fill
    panel.grid.major = element_line(color = "white", linewidth = 0.2), # Soft white grid lines
    axis.text = element_text(size = 5.0),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(4, 4, 4, 4, "pt"),
    panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5)
  )

# --- Panel B: Amazon Spatial Clustering Map ---
p2_b <- ggplot() +
  geom_spatvector(data = countries_amazon_map, fill = "#F4F6F7", colour = "grey70", linewidth = 0.2) +
  geom_point(data = amazon_c_info$centroids, aes(x = lon, y = lat, size = total_trap_days), 
             shape = 21, color = "black", fill = pal_region[["Amazon"]], stroke = 0.4, alpha = 0.5) +
  scale_size_continuous(name = "Trap Days", range = c(1.0, 5.0), breaks = c(100, 500, 1000, 2000, 5000), limits = c(10, 100000)) +
  coord_sf(xlim = c(-85, -35), ylim = c(-15, 15), expand = FALSE) +
  labs(title = "A") +
  map_p_theme +
  theme(legend.position = "none")
 
# --- Panel C: Congo Spatial Clustering Map ---
p2_c <- ggplot() +
  geom_spatvector(data = countries_congo_map, fill = "#F4F6F7", colour = "grey70", linewidth = 0.2) +
  geom_point(data = congo_c_info$centroids, aes(x = lon, y = lat, size = total_trap_days), 
             shape = 21, color = "black", fill = pal_region[["Congo"]], stroke = 0.4, alpha = 0.5) +
  scale_size_continuous(name = "Trap Days", range = c(1.0, 5.0), breaks = c(100, 500, 1000, 2000, 5000), limits = c(10, 100000)) +
  coord_sf(xlim = c(-5, 45), ylim = c(-15, 15), expand = FALSE) +
  labs(title = "B") +
  map_p_theme +
  theme(legend.position = "none")
 
# --- Panel D: Southeast Asia Spatial Clustering Map ---
p2_g <- ggplot() +
  geom_spatvector(data = countries_seasia_map, fill = "#F4F6F7", colour = "grey70", linewidth = 0.2) +
  geom_point(data = seasia_c_info$centroids, aes(x = lon, y = lat, size = total_trap_days), 
             shape = 21, color = "black", fill = pal_region[["SE_Asia"]], stroke = 0.4, alpha = 0.5) +
  scale_size_continuous(name = "Trap Days", range = c(1.0, 5.0), breaks = c(100, 500, 1000, 2000, 5000), limits = c(10, 100000)) +
  coord_sf(xlim = c(90, 140), ylim = c(-15, 15), expand = FALSE) +
  labs(title = "C") +
  map_p_theme +
  theme(
    legend.position = c(0.02, 0.02),
    legend.justification = c(0, 0),
    legend.title = element_text(size = 5.0, face = "bold"),
    legend.text = element_text(size = 4.5),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.key = element_blank(),
    legend.key.size = unit(0.2, "cm")
  ) +
  guides(size = guide_legend(override.aes = list(fill = "grey70", alpha = 0.7)))

# --- Panel E: Temporal Calibration Timeline ---
gedi_start <- as.Date("2019-04-17")

cluster_temporal_spans <- det_all %>%
  select(region, cluster_id_geo = cluster_id, deployment_id, start_date, end_date, trap_days) %>%
  distinct() %>%
  mutate(
    start_date = as.Date(start_date),
    end_date = as.Date(end_date),
    years_before_gedi = as.numeric(gedi_start - start_date) / 365.25,
    w_temp = case_when(
      years_before_gedi <= brackets[1] ~ weights[1],
      years_before_gedi <= brackets[2] ~ weights[2],
      years_before_gedi <= brackets[3] ~ weights[3],
      TRUE                             ~ weights[4]
    )
  )

# Group by cluster to get aggregate spans
cluster_spans <- cluster_temporal_spans %>%
  group_by(region, cluster_id_geo) %>%
  summarise(
    start = min(start_date),
    end = max(end_date),
    w_temp_cluster = sum(trap_days * w_temp) / sum(trap_days),
    .groups = "drop"
  ) %>%
  arrange(start) %>%
  mutate(row_idx = row_number())

p2_d <- ggplot(cluster_spans) +
  geom_vline(xintercept = gedi_start, linetype = "solid", color = "black", linewidth = 0.75) +
  annotate("text", x = gedi_start - 180, y = 90, label = "GEDI Launch\n(April 2019)", color = "black", size = 2.0, fontface = "bold", hjust = 1) +
  geom_segment(aes(x = start, xend = end, y = row_idx, yend = row_idx, color = region), linewidth = 0.6, alpha = 0.8) +
  scale_color_manual(values = pal_region, name = "Region") +
  scale_x_date(date_breaks = "4 years", date_labels = "%Y", limits = c(as.Date("2003-01-01"), as.Date("2024-12-31"))) +
  labs(
    title = "E",
    x = "Survey Era",
    y = "Camera Trap Clusters"
  ) +
  t_theme +
  theme(aspect.ratio = 0.45,
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(4, 4, 4, 4, "pt"))

# --- Panel G: Taxonomic order composition horizontal bar plot ---
# ORDER_COLOURS now provided by pal_taxon from theme_pnas.R

df_tax <- read_csv("outputs/camera_traps_robust_taxon_biomass.csv", show_col_types = FALSE) %>%
  filter(!is.na(order) & !is.na(body_mass_kg) & body_mass_kg > 0) %>%
  group_by(region, order) %>%
  summarise(total_biomass = sum(biomass_contrib, na.rm = TRUE), .groups = "drop")



top_orders <- c("Cetartiodactyla", "Proboscidea", "Carnivora", "Rodentia", "Primates", "Cingulata", "Perissodactyla")
df_tax <- df_tax %>%
  mutate(order_clean = ifelse(order %in% top_orders, order, "Other")) %>%
  group_by(region, order_clean) %>%
  summarise(total_biomass = sum(total_biomass), .groups = "drop")

df_tax <- df_tax %>%
  group_by(region) %>%
  mutate(prop = total_biomass / sum(total_biomass) * 100) %>%
  ungroup()

df_tax$order_clean <- factor(df_tax$order_clean, levels = rev(c(top_orders, "Other")))

p2_e <- ggplot(df_tax, aes(x = prop, y = region, fill = order_clean)) +
  geom_bar(stat = "identity", width = 0.85, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = pal_taxon, name = "Taxonomic Order:") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(labels = c("Amazon" = "Amazon", "Congo" = "Congo", "SE_Asia" = "SE Asia"), expand = c(0, 0.3)) +
  labs(
    title = "G",
    x = expression(Proportion~of~Total~Biomass~(B[H]*" %")),
    y = NULL
  ) +
  t_theme +
  theme(aspect.ratio = 0.45,
        legend.position = "bottom",
        legend.title = element_text(size = 5.0, face = "bold"),
        legend.text = element_text(size = 4.5),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(t = -2, b = -2, unit = "pt"),
        legend.spacing.x = unit(0.1, "cm"),
        legend.spacing.y = unit(0.05, "cm"),
        plot.margin = margin(4, 4, 4, 4, "pt")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# --- Panel F: Vertebrate Biomass Index distribution across clusters ---
df_clusters <- read_csv("outputs/camera_traps_cluster_level_metrics_robust.csv", show_col_types = FALSE)

bin_width <- 0.4
df_binned_biomass <- df_clusters %>%
  filter(!is.na(B_H_index)) %>%
  mutate(log_biomass = log10(B_H_index + 1)) %>%
  mutate(bin_center = round(log_biomass / bin_width) * bin_width) %>%
  group_by(region, bin_center) %>%
  summarise(n_clusters = n(), .groups = "drop") %>%
  complete(region, bin_center = seq(0, max(bin_center, na.rm = TRUE), by = bin_width), fill = list(n_clusters = 0)) %>%
  mutate(B_H_index_val = 10^bin_center - 1)

p2_f <- ggplot(df_binned_biomass, aes(x = B_H_index_val, y = n_clusters, color = region)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2, alpha = 0.8) +
  scale_color_manual(values = pal_region, labels = c("Amazon" = "Amazon", "Congo" = "Congo", "SE_Asia" = "SE Asia"), name = NULL) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 10000), labels = c("0", "10", "100", "1K", "10K")) +
  labs(
    title = "F",
    x = expression(Cluster~Biomass~Index~(B[H]*", log1p scale")),
    y = "Number of Clusters"
  ) +
  t_theme +
  theme(
    aspect.ratio = 0.45,
    legend.position = "none",
    plot.margin = margin(4, 4, 4, 4, "pt")
  )

# Assemble Figure 2 using patchwork for perfect panel alignment
col1 <- p2_b / p2_c / p2_g + plot_layout(heights = c(1, 1, 1))
col2 <- p2_a / p2_d / p2_f / p2_e + plot_layout(heights = c(0.9, 1.1, 0.9, 1.1))

fig2_final <- (col1 | col2) + plot_layout(widths = c(1.2, 1.0))

fig2_png <- "figures/figure2_camera_trap_pipeline.png"
ggsave(filename = fig2_png, plot = fig2_final, width = 22.0, height = 23.0, units = "cm", dpi = 600, bg = "white")
cat("✓ 7-Panel Figure 2 successfully saved to:", fig2_png, "\n\n")

fig2_pdf <- "figures/figure2_camera_trap_pipeline.pdf"
ggsave(filename = fig2_pdf, plot = fig2_final, width = 22.0, height = 23.0, units = "cm", dpi = 600, bg = "white")
cat("✓ 7-Panel Figure 2 PDF successfully saved to:", fig2_pdf, "\n\n")




cat("============================================================\n")
cat("=== Manuscript Figures 1 & 2 Completed Successfully ===\n")
cat("============================================================\n")
