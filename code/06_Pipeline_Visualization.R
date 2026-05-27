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
library(dplyr)
library(readr)
library(sf)
library(scales)
library(tidyr)

cat("============================================================\n")
cat("=== Generating 6-Panel Manuscript Figures 1 and 2 (Pipeline Summary) ===\n")
cat("============================================================\n\n")

# --- Source standard styling and helpers -------------------------------------
source("code/functions/theme_pnas.R")
source("code/functions/calibration_helpers.R")

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# Color palettes
pal_basin <- c("Amazon" = "#E65100", "Congo" = "#1B5E20", "SE_Asia" = "#0D47A1")

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

# Load 5km aggregated maps
r_congo_5km <- assign_multiscale_names(rast("outputs/EOdata/analysis_stack_5000_Congo.tif"))
r_amazon_5km <- assign_multiscale_names(rast("outputs/EOdata/analysis_stack_5000_Amazon.tif"))

# Load country outlines and crop for background
countries_v <- vect("data/world-administrative-boundaries")

# Load camera trap cluster coordinates for overlay
mcps <- vect("outputs/camera_traps_robust_buffered_mcps.geojson")
mcps_congo <- mcps[mcps$region == "Congo", ]
mcps_amazon <- mcps[mcps$region == "Amazon", ]

# Set up study bounding boxes with padding
ext_c <- ext(mcps_congo)
ext_c <- ext(xmin(ext_c)-1.0, xmax(ext_c)+1.0, ymin(ext_c)-1.0, ymax(ext_c)+1.0)

ext_a <- ext(mcps_amazon)
ext_a <- ext(xmin(ext_a)-1.0, xmax(ext_a)+1.0, ymin(ext_a)-1.0, ymax(ext_a)+1.0)

r_c_crop <- crop(r_congo_5km, ext_c)
r_a_crop <- crop(r_amazon_5km, ext_a)

# Scale GEDI shot density to raw shot count (number of GEDI footprints in 5km cell)
r_c_crop$gedi_n <- r_c_crop$gedi_n * 40000
r_a_crop$gedi_n <- r_a_crop$gedi_n * 40000

countries_c <- crop(project(countries_v, crs(r_c_crop)), ext_c)
countries_a <- crop(project(countries_v, crs(r_a_crop)), ext_a)

# --- Figure 1 Panels A & B: Regional 5km UOI maps ---
t_theme <- theme_pnas(base_size = 7.5)

p1_a <- ggplot() +
  geom_spatraster(data = r_c_crop, aes(fill = uoi)) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "GEDI Understory Openness Index (UOI)",
    limits = c(0.92, 0.975),
    oob = scales::squish,
    na.value = "transparent"
  ) +
  geom_spatvector(data = countries_c, fill = NA, colour = "grey80", linewidth = 0.25) +
  geom_spatvector(data = mcps_congo, fill = NA, color = "white", linewidth = 0.4, linetype = "solid") +
  t_theme +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(2, 2, 2, 2, "pt")
  ) +
  labs(title = "A. Congo Basin GEDI Openness (5 km)")

p1_b <- ggplot() +
  geom_spatraster(data = r_a_crop, aes(fill = uoi)) +
  scale_fill_viridis_c(
    option = "plasma",
    name = "GEDI Understory Openness Index (UOI)",
    limits = c(0.92, 0.975),
    oob = scales::squish,
    na.value = "transparent",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      barwidth = unit(5.5, "cm"),
      barheight = unit(0.12, "cm")
    )
  ) +
  geom_spatvector(data = countries_a, fill = NA, colour = "grey80", linewidth = 0.25) +
  geom_spatvector(data = mcps_amazon, fill = NA, color = "white", linewidth = 0.4, linetype = "solid") +
  t_theme +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 5.5, face = "bold"),
    legend.text = element_text(size = 4.5),
    legend.margin = margin(t = -2, b = 2, unit = "pt"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(2, 2, 2, 2, "pt")
  ) +
  labs(title = "B. Amazon Basin GEDI Openness (5 km)")

# --- Figure 1 Panels C & D: Regional 5km GEDI Shot Density maps ---
p1_c <- ggplot() +
  geom_spatraster(data = r_c_crop, aes(fill = gedi_n)) +
  scale_fill_viridis_c(
    option = "mako",
    name = "GEDI Shot Count (log scale)",
    trans = "log10",
    limits = c(10, 6000),
    breaks = c(10, 100, 1000, 6000),
    labels = c("10", "100", "1,000", "6,000"),
    oob = scales::squish,
    na.value = "transparent"
  ) +
  geom_spatvector(data = countries_c, fill = NA, colour = "grey80", linewidth = 0.25) +
  geom_spatvector(data = mcps_congo, fill = NA, color = "white", linewidth = 0.4, linetype = "solid") +
  t_theme +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(2, 2, 2, 2, "pt")
  ) +
  labs(title = "C. Congo Basin GEDI Shot Count (5 km)")

p1_d <- ggplot() +
  geom_spatraster(data = r_a_crop, aes(fill = gedi_n)) +
  scale_fill_viridis_c(
    option = "mako",
    name = "GEDI Shot Count (log scale)",
    trans = "log10",
    limits = c(10, 6000),
    breaks = c(10, 100, 1000, 6000),
    labels = c("10", "100", "1,000", "6,000"),
    oob = scales::squish,
    na.value = "transparent",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      barwidth = unit(5.5, "cm"),
      barheight = unit(0.12, "cm")
    )
  ) +
  geom_spatvector(data = countries_a, fill = NA, colour = "grey80", linewidth = 0.25) +
  geom_spatvector(data = mcps_amazon, fill = NA, color = "white", linewidth = 0.4, linetype = "solid") +
  t_theme +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 5.5, face = "bold"),
    legend.text = element_text(size = 4.5),
    legend.margin = margin(t = -2, b = 2, unit = "pt"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(2, 2, 2, 2, "pt")
  ) +
  labs(title = "D. Amazon Basin GEDI Shot Count (5 km)")

# --- Figure 1 Panel E: Uncertainty Decay vs. GEDI Shot Density ---
df_c_pts <- as.data.frame(r_c_crop[[c("uoi_sd", "gedi_n")]], na.rm = TRUE)
df_c_pts$basin <- "Congo"
df_a_pts <- as.data.frame(r_a_crop[[c("uoi_sd", "gedi_n")]], na.rm = TRUE)
df_a_pts$basin <- "Amazon"
df_pts <- rbind(df_c_pts, df_a_pts)

# Downsample for clean plotting density
set.seed(42)
df_pts <- df_pts %>%
  group_by(basin) %>%
  slice_sample(n = 2000) %>%
  ungroup()

p1_e <- ggplot(df_pts, aes(x = gedi_n, y = uoi_sd, color = basin)) +
  geom_point(alpha = 0.25, size = 0.5, stroke = 0) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE, linewidth = 0.75) +
  scale_color_manual(values = pal_basin, name = "Basin") +
  scale_x_continuous(trans = "log10", labels = comma_format()) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), limits = c(0, 0.008), oob = scales::squish) +
  labs(
    title = "E. Uncertainty Decay vs. GEDI Shot Count",
    x = "GEDI Shot Count (log scale)",
    y = "GEDI UOI Standard Error (SE)"
  ) +
  t_theme +
  theme(
    legend.position = c(0.76, 0.78),
    legend.title = element_text(size = 5.5, face = "bold"),
    legend.text = element_text(size = 5.0),
    legend.key.height = unit(0.12, "cm"),
    legend.key.width = unit(0.3, "cm"),
    plot.margin = margin(4, 4, 4, 4, "pt")
  )

# --- Figure 1 Panel F: UOI Distribution divided by Continent ---
df_c_uoi <- data.frame(uoi = as.data.frame(r_c_crop[["uoi"]], na.rm = TRUE)$uoi, basin = "Congo")
df_a_uoi <- data.frame(uoi = as.data.frame(r_a_crop[["uoi"]], na.rm = TRUE)$uoi, basin = "Amazon")
df_uoi_dist <- rbind(df_c_uoi, df_a_uoi)

p1_f <- ggplot(df_uoi_dist, aes(x = uoi, fill = basin, color = basin)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  scale_fill_manual(values = pal_basin, name = "Basin") +
  scale_color_manual(values = pal_basin, name = "Basin") +
  scale_x_continuous(limits = c(0.92, 0.975), breaks = seq(0.92, 0.97, by = 0.01)) +
  labs(
    title = "F. UOI Distribution by Continent",
    x = "GEDI Understory Openness Index (UOI)",
    y = "Probability Density"
  ) +
  t_theme +
  theme(
    legend.position = "none",
    plot.margin = margin(4, 4, 4, 4, "pt")
  )

# Assemble Figure 1 (6 Panels arranged in a flat perfectly aligned 3x2 Grid)
fig1_final <- cowplot::plot_grid(
  p1_a, p1_b,
  p1_c, p1_d,
  p1_e, p1_f,
  ncol = 2,
  align = "hv",
  axis = "tblr",
  rel_heights = c(1.0, 1.0, 1.0)
)

fig1_png <- "figures/figure1_gedi_pipeline.png"
ggsave(filename = fig1_png, plot = fig1_final, width = 17.8, height = 22.0, units = "cm", dpi = 600, bg = "white")
cat("✓ 6-Panel Figure 1 successfully saved to:", fig1_png, "\n\n")


# =============================================================================
# FIGURE 2: CAMERA TRAP WILDLIFE INGESTION, CLUSTERING, AND TEMPORAL CALIBRATION
# =============================================================================
cat("--- Drafting Figure 2: Camera Trap Wildlife Ingestion & Calibration (7 Panels) ---\n")

# Load raw detections
det_path <- "outputs/camera_traps_joint_detections.csv"
if (!file.exists(det_path)) {
  stop("Camera trap detections missing: ", det_path)
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
  scale_color_manual(values = pal_basin, name = "Basin") +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  labs(
    title = "A. Vertebrate Individual Size Distribution",
    x = "Mammal Body Mass (kg, log scale)",
    y = "Total Detections (log scale)"
  ) +
  t_theme +
  theme(legend.position = c(0.78, 0.82),
        legend.title = element_text(size = 6.0, face = "bold"),
        legend.text = element_text(size = 5.5),
        plot.margin = margin(2, 4, 2, 4, "pt"))

# --- Dynamic Haversine Clustering Maps Helper ---
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

get_clustered_deployments <- function(basin_name, buffer_val = 0.1) {
  unique_deps <- det_all %>%
    filter(region == basin_name) %>%
    select(longitude, latitude, deployment_id, trap_days) %>%
    distinct()
  
  n_deps <- nrow(unique_deps)
  if (n_deps > 1) {
    dist_mat <- matrix(0, nrow = n_deps, ncol = n_deps)
    for (i in 1:n_deps) {
      for (j in 1:n_deps) {
        dist_mat[i,j] <- haversine_dist(unique_deps$longitude[i], unique_deps$latitude[i], unique_deps$longitude[j], unique_deps$latitude[j])
      }
    }
    hc <- hclust(as.dist(dist_mat), method = "single")
    unique_deps$cluster_num <- cutree(hc, h = CLUSTER_THRESHOLD_KM)
  } else {
    unique_deps$cluster_num <- 1
  }
  
  # Calculate centroids and sum trap days per cluster
  centroids <- unique_deps %>%
    group_by(cluster_num) %>%
    summarise(
      lon = mean(longitude),
      lat = mean(latitude),
      total_trap_days = sum(trap_days, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(list(deps = unique_deps, centroids = centroids))
}

# Run Haversine Single-Linkage Clustering for all three basins
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
    plot.margin = margin(2, 2, 2, 2, "pt"),
    panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5)
  )

# --- Panel B: Amazon Spatial Clustering Map ---
p2_b <- ggplot() +
  geom_spatvector(data = countries_amazon_map, fill = "#F4F6F7", colour = "grey70", linewidth = 0.2) +
  geom_point(data = amazon_c_info$centroids, aes(x = lon, y = lat, size = total_trap_days), 
             shape = 21, color = "black", fill = "#E65100", stroke = 0.4, alpha = 0.5) +
  scale_size_continuous(name = "Trap Days", range = c(1.0, 5.0), breaks = c(100, 500, 1000, 2000, 5000), limits = c(10, 100000)) +
  coord_sf(xlim = c(-85, -35), ylim = c(-15, 15), expand = FALSE) +
  labs(title = sprintf("B. Amazon Clusters (n = %d)", nrow(amazon_c_info$centroids))) +
  map_p_theme +
  theme(legend.position = "none")

# --- Panel C: Congo Spatial Clustering Map ---
p2_c <- ggplot() +
  geom_spatvector(data = countries_congo_map, fill = "#F4F6F7", colour = "grey70", linewidth = 0.2) +
  geom_point(data = congo_c_info$centroids, aes(x = lon, y = lat, size = total_trap_days), 
             shape = 21, color = "black", fill = "#1B5E20", stroke = 0.4, alpha = 0.5) +
  scale_size_continuous(name = "Trap Days", range = c(1.0, 5.0), breaks = c(100, 500, 1000, 2000, 5000), limits = c(10, 100000)) +
  coord_sf(xlim = c(-5, 45), ylim = c(-15, 15), expand = FALSE) +
  labs(title = sprintf("C. Congo Clusters (n = %d)", nrow(congo_c_info$centroids))) +
  map_p_theme +
  theme(legend.position = "none")

# --- Panel D: Southeast Asia Spatial Clustering Map ---
p2_g <- ggplot() +
  geom_spatvector(data = countries_seasia_map, fill = "#F4F6F7", colour = "grey70", linewidth = 0.2) +
  geom_point(data = seasia_c_info$centroids, aes(x = lon, y = lat, size = total_trap_days), 
             shape = 21, color = "black", fill = "#0D47A1", stroke = 0.4, alpha = 0.5) +
  scale_size_continuous(name = "Trap Days", range = c(1.0, 5.0), breaks = c(100, 500, 1000, 2000, 5000), limits = c(10, 100000)) +
  coord_sf(xlim = c(90, 140), ylim = c(-15, 15), expand = FALSE) +
  labs(title = sprintf("D. SE Asia Clusters (n = %d)", nrow(seasia_c_info$centroids))) +
  map_p_theme +
  theme(
    legend.position = c(0.18, 0.22),
    legend.title = element_text(size = 5.0, face = "bold"),
    legend.text = element_text(size = 4.5),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.key = element_blank(),
    legend.key.size = unit(0.2, "cm")
  )

# --- Panel E: Temporal Calibration Timeline ---
gedi_start <- as.Date("2019-04-17")

cluster_temporal_spans <- det_all %>%
  select(region, cluster_id_geo = project_id, deployment_id, start_date, end_date, trap_days) %>%
  distinct() %>%
  mutate(
    start_date = as.Date(start_date),
    end_date = as.Date(end_date),
    years_before_gedi = as.numeric(gedi_start - start_date) / 365.25,
    w_temp = case_when(
      years_before_gedi <= 1.0  ~ 1.0,
      years_before_gedi <= 6.0  ~ 0.5,
      years_before_gedi <= 11.0 ~ 0.25,
      TRUE                      ~ 0.1
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
  geom_vline(xintercept = gedi_start, linetype = "solid", color = "#D32F2F", linewidth = 0.75) +
  annotate("text", x = gedi_start + 180, y = 15, label = "GEDI Launch\n(April 2019)", color = "#D32F2F", size = 2.0, fontface = "bold", hjust = 0) +
  geom_segment(aes(x = start, xend = end, y = row_idx, yend = row_idx, color = region, alpha = w_temp_cluster), linewidth = 1.5) +
  scale_color_manual(values = pal_basin, name = "Basin") +
  scale_alpha_continuous(name = "Temporal Weight", range = c(0.35, 1.0), breaks = c(0.1, 0.25, 0.5, 1.0)) +
  scale_x_date(date_breaks = "4 years", date_labels = "%Y", limits = c(as.Date("2003-01-01"), as.Date("2024-12-31"))) +
  labs(
    title = "E. Survey Temporal Calibration Span",
    x = "Survey Era",
    y = "Spatial Clusters (Ranked by Start Date)"
  ) +
  t_theme +
  theme(legend.position = "bottom",
        legend.box = "vertical", # Stack legends vertically to prevent overflow
        legend.title = element_text(size = 5.0, face = "bold"),
        legend.text = element_text(size = 4.5),
        legend.key.height = unit(0.10, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.margin = margin(t = -2, b = -2, unit = "pt"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(2, 4, 2, 4, "pt"))

# --- Panel F: Taxonomic order composition horizontal bar plot ---
ORDER_COLOURS <- c(
  "Cetartiodactyla"  = "#1B7837",  # Green
  "Proboscidea"      = "#2166AC",  # Blue
  "Carnivora"        = "#E08214",  # Orange
  "Rodentia"         = "#6A3D9A",  # Purple
  "Primates"         = "#C51B7D",  # Pink
  "Cingulata"        = "#35978F",  # Teal
  "Perissodactyla"   = "#4D4D4D",  # Charcoal
  "Other"            = "#878787"   # Grey
)

df_tax <- det_all %>%
  filter(!is.na(order) & !is.na(body_mass_kg) & body_mass_kg > 0) %>%
  group_by(region, order) %>%
  summarise(total_biomass = sum(n_detections * body_mass_kg), .groups = "drop")

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
  geom_bar(stat = "identity", width = 0.55, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = ORDER_COLOURS, name = "Taxonomic Order:") +
  labs(
    title = "F. Vertebrate Biomass Composition",
    x = "Proportion of Total Biomass (%)",
    y = "Basin/Continent"
  ) +
  t_theme +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 4.5, face = "bold"),
        legend.text = element_text(size = 4.0),
        legend.key.width = unit(0.12, "cm"),
        legend.key.height = unit(0.12, "cm"),
        legend.margin = margin(t = -2, b = -2, unit = "pt"),
        legend.spacing.x = unit(0.08, "cm"),
        legend.spacing.y = unit(0.05, "cm"),
        plot.margin = margin(2, 4, 2, 4, "pt")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# --- Panel G: Vertebrate Biomass Index distribution across clusters ---
df_clusters <- read_csv("outputs/camera_traps_cluster_level_metrics.csv", show_col_types = FALSE) %>%
  mutate(region = ifelse(region == "SE_Asia", "SE Asia", region))

p2_f <- ggplot(df_clusters, aes(x = B_H_index, fill = region)) +
  geom_histogram(bins = 10, color = "black", linewidth = 0.25, show.legend = FALSE) +
  scale_fill_manual(values = c("Amazon" = "#E65100", "Congo" = "#1B5E20", "SE Asia" = "#0D47A1")) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 100000, 10000000), labels = c("0", "10", "100", "1K", "100K", "10M")) +
  facet_wrap(~region, ncol = 1) +
  labs(
    title = "G. Standing Biomass Across Spatial Clusters",
    x = "Cluster Biomass Index (log1p scale)",
    y = "Number of Clusters"
  ) +
  t_theme +
  theme(
    strip.text = element_text(size = 5.0, face = "bold", margin = margin(1, 1, 1, 1)),
    plot.margin = margin(2, 4, 8, 4, "pt")
  )

# Assemble Figure 2: Column of 3 maps on the left, other 4 plots on the right
column_maps <- cowplot::plot_grid(
  p2_b, p2_c, p2_g,
  ncol = 1,
  align = "v"
)

column_right <- cowplot::plot_grid(
  p2_a,
  p2_d,
  p2_e,
  p2_f,
  ncol = 1,
  align = "v",
  rel_heights = c(0.9, 1.1, 1.1, 0.9)
)

fig2_final <- cowplot::plot_grid(
  column_maps,
  column_right,
  ncol = 2,
  rel_widths = c(1.2, 1.0)
)

fig2_png <- "figures/figure2_camera_trap_pipeline.png"
ggsave(filename = fig2_png, plot = fig2_final, width = 17.8, height = 22.0, units = "cm", dpi = 600, bg = "white")
cat("✓ 7-Panel Figure 2 successfully saved to:", fig2_png, "\n\n")

# Copy figures to active brain artifacts folder
brain_artifacts_dir <- "/home/j/.gemini/antigravity/brain/913e5cea-7c99-4b21-8124-ea8455da8457"
if (file.exists(brain_artifacts_dir)) {
  file.copy(fig2_png, file.path(brain_artifacts_dir, "figure2_camera_trap_pipeline.png"), overwrite = TRUE)
  cat("✓ Copied 7-Panel Figure 2 to brain artifacts folder.\n")
}

cat("============================================================\n")
cat("=== Manuscript Figures 1 & 2 Completed Successfully ===\n")
cat("============================================================\n")
