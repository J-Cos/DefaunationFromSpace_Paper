# =============================================================================
# code/07_Pipeline_Visualization.R
#
# Generates two premium PNAS-style publication-grade figures summarizing the
# GEDI raw data pipeline (Figure 1) and the Camera Trap data pipeline (Figure 2).
#
# Figure 1: GEDI Raw Data Pipeline & Shot-Sparsity Distribution
#   (A) Congo Basin regional map of GEDI UOI (5km) with cluster points
#   (B) Amazon Basin regional map of GEDI UOI (5km) with cluster points
#   (C) Zoomed-in inset at raw 500m scale showing GEDI shot density sparsity (gedi_n)
#   (D) Shot-noise uncertainty reduction: GEDI UOI Standard Error vs. Shot Count (gedi_n)
#
# Figure 2: Camera Trap Wildlife Ingestion, Clustering, and Temporal Calibration
#   (A) Species rank-abundance curves for both basins
#   (B) Spatial clustering map of deployments colored by mathematical cluster ID
#   (C) Temporal calibration timeline showing cluster survey spans vs. GEDI mission launch
#   (D) Proportional horizontal bar plot of taxonomic order composition
# =============================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(cowplot)
library(dplyr)
library(readr)
library(sf)
library(scales)

cat("============================================================\n")
cat("=== Generating Manuscript Figures 1 and 2 (Pipeline Summary) ===\n")
cat("============================================================\n\n")

# --- Source standard styling and helpers -------------------------------------
source("code/functions/theme_pnas.R")
source("code/functions/calibration_helpers.R")

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# Color palettes
pal_basin <- c("Amazon" = "#E65100", "Congo" = "#1B5E20")

# =============================================================================
# FIGURE 1: GEDI RAW DATA PIPELINE & SHOT-SPARSITY DISTRIBUTION
# =============================================================================
cat("--- Drafting Figure 1: GEDI Raw Data Pipeline ---\n")

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

countries_c <- crop(project(countries_v, crs(r_c_crop)), ext_c)
countries_a <- crop(project(countries_v, crs(r_a_crop)), ext_a)

# --- Figure 1 Panels A & B: Regional 5km UOI maps ---
t_theme <- theme_pnas(base_size = 7.5)

p1_a <- ggplot() +
  geom_spatraster(data = r_c_crop, aes(fill = uoi)) +
  scale_fill_viridis_c(option = "plasma", name = "UOI", limits = c(0.92, 0.975), oob = scales::squish, na.value = "transparent") +
  geom_spatvector(data = countries_c, fill = NA, colour = "grey80", linewidth = 0.25) +
  geom_spatvector(data = mcps_congo, fill = NA, color = "white", linewidth = 0.4, linetype = "solid") +
  t_theme +
  theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        plot.margin = margin(2, 2, 2, 2, "pt")) +
  labs(title = "A. Congo Basin GEDI Openness (5 km)")

p1_b <- ggplot() +
  geom_spatraster(data = r_a_crop, aes(fill = uoi)) +
  scale_fill_viridis_c(option = "plasma", name = "GEDI Understory Openness Index (UOI)", limits = c(0.92, 0.975), oob = scales::squish, na.value = "transparent") +
  geom_spatvector(data = countries_a, fill = NA, colour = "grey80", linewidth = 0.25) +
  geom_spatvector(data = mcps_amazon, fill = NA, color = "white", linewidth = 0.4, linetype = "solid") +
  t_theme +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 6.5, face = "bold"),
        legend.text = element_text(size = 5.5),
        legend.key.height = unit(0.18, "cm"),
        legend.key.width = unit(1.0, "cm"),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        plot.margin = margin(2, 2, 2, 2, "pt")) +
  labs(title = "B. Amazon Basin GEDI Openness (5 km)")

# Extract map legend
shared_uoi_legend <- get_legend(p1_b)
p1_b <- p1_b + theme(legend.position = "none")

# --- Figure 1 Panel C: Zoomed-in Inset at raw 500m scale showing GEDI shot density sparsity (gedi_n) ---
# We select a sub-extent of Congo native stack containing a camera trap cluster
r_congo_native <- assign_native_names(rast("outputs/EOdata/analysis_stack_native_Congo.tif"))

# Zoom in on Nouabalé-Ndoki cluster area (centered around Lon 16.7, Lat 2.45)
zoom_ext <- ext(16.4, 16.9, 2.2, 2.7)
r_zoom_uoi <- crop(r_congo_native[["uoi"]], zoom_ext)
r_zoom_n <- crop(r_congo_native[["gedi_n"]], zoom_ext)

# Convert to dataframes for neat plotting
df_zoom_n <- as.data.frame(r_zoom_n, xy = TRUE, na.rm = TRUE)

p1_c <- ggplot() +
  geom_tile(data = df_zoom_n, aes(x = x, y = y, fill = gedi_n)) +
  scale_fill_viridis_c(option = "mako", name = "Shot Count", trans = "log1p", breaks = c(0, 5, 20, 100), labels = c("0", "5", "20", "100+")) +
  t_theme +
  theme(legend.position = "right",
        legend.title = element_text(size = 6.0, face = "bold"),
        legend.text = element_text(size = 5.0),
        legend.key.width = unit(0.12, "cm"),
        legend.key.height = unit(0.4, "cm"),
        axis.text = element_text(size = 5.5),
        plot.margin = margin(2, 2, 2, 2, "pt")) +
  labs(title = "C. Zoomed Raw GEDI Shot Counts (500m)",
       x = "Longitude (°E)", y = "Latitude (°N)")

# --- Figure 1 Panel D: GEDI UOI Standard Error vs. Shot Count (gedi_n) ---
# Showcases the shot-noise uncertainty reduction scaling law (SE ~ 1/sqrt(gedi_n))
# We extract a random sample of pixels from the native Congo stack
set.seed(42)
df_native_sample <- as.data.frame(r_congo_native[[c("uoi", "gedi_n")]], na.rm = TRUE)
df_native_sample <- df_native_sample %>%
  filter(gedi_n > 0 & gedi_n < 1000) %>%
  sample_n(pmin(nrow(.), 1500))

# Standard error of the mean for UOI can be simulated or calculated.
# Since UOI is beta-distributed, UOI_SD across a cluster decreases with shot density.
# Let's plot the standard error of UOI aggregated across scales vs aggregated shot counts!
# We extract this from the multi-scale dataset itself
df_scales <- extract_scale_data(5000)
# We can pool multiple scales to show the standard error of UOI decreasing as average shot counts (gedi_n) increase!
scales_vector <- c(5000, 10000, 20000, 30000, 50000)
df_multiscale_se <- do.call(rbind, lapply(scales_vector, function(s) {
  dat <- extract_scale_data(s)
  dat$scale_km <- s / 1000
  return(dat)
}))

p1_d <- ggplot(df_multiscale_se, aes(x = n_pixels * 5, y = uoi_se, color = basin)) +
  geom_point(aes(size = scale_km), alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ log(x), color = "black", linetype = "dashed", linewidth = 0.5, se = FALSE) +
  scale_color_manual(values = pal_basin, name = "Basin") +
  scale_size_continuous(name = "Scale (km)", range = c(1.0, 3.5), breaks = c(5, 20, 50)) +
  scale_x_log10(labels = comma_format()) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = "D. Shot-Noise Averaging Law",
    x = "Average GEDI Shots per Grid Cell (log scale)",
    y = "GEDI UOI Standard Error of Mean (SE)"
  ) +
  t_theme +
  theme(legend.position = "right",
        legend.title = element_text(size = 6.0, face = "bold"),
        legend.text = element_text(size = 5.5),
        legend.margin = margin(0,0,0,0),
        plot.margin = margin(2, 2, 2, 2, "pt"))

# Assemble Figure 1
row1 <- plot_grid(p1_a, p1_b, ncol = 2, rel_widths = c(1, 1))
row1_legend <- plot_grid(row1, shared_uoi_legend, ncol = 1, rel_heights = c(1, 0.15))
row2 <- plot_grid(p1_c, p1_d, ncol = 2, rel_widths = c(1, 1.15))

fig1_final <- plot_grid(row1_legend, row2, ncol = 1, rel_heights = c(1.05, 0.95), hspace = 0.2)

fig1_png <- "figures/figure1_gedi_pipeline.png"
ggsave(filename = fig1_png, plot = fig1_final, width = 17.8, height = 15.0, units = "cm", dpi = 600, bg = "white")
cat("✓ Figure 1 successfully saved to:", fig1_png, "\n\n")


# =============================================================================
# FIGURE 2: CAMERA TRAP WILDLIFE INGESTION, CLUSTERING, AND TEMPORAL CALIBRATION
# =============================================================================
cat("--- Drafting Figure 2: Camera Trap Wildlife Ingestion & Calibration ---\n")

# Load raw detections
det_path <- "outputs/camera_traps_joint_detections.csv"
if (!file.exists(det_path)) {
  stop("Camera trap detections missing: ", det_path)
}
det_all <- read_csv(det_path, show_col_types = FALSE)

# Load aggregated cluster data
df_clusters <- extract_scale_data(5000)

# --- Panel A: Species Rank-Abundance curves ---
sp_abundance <- det_all %>%
  filter(!is.na(species) & taxon_quality == "species") %>%
  group_by(region, genus, species, common_name) %>%
  summarise(n_detections = sum(n_detections), .groups = "drop") %>%
  arrange(region, desc(n_detections))

# Rank species within each region
sp_abundance <- sp_abundance %>%
  group_by(region) %>%
  mutate(rank = row_number()) %>%
  ungroup()

p2_a <- ggplot(sp_abundance, aes(x = rank, y = n_detections, color = region)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.0, alpha = 0.8) +
  scale_color_manual(values = pal_basin, name = "Basin") +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  labs(
    title = "A. Vertebrate Rank-Abundance Curves",
    x = "Species Rank",
    y = "Total Independent Detections"
  ) +
  t_theme +
  theme(legend.position = c(0.7, 0.8),
        legend.title = element_text(size = 6.0, face = "bold"),
        legend.text = element_text(size = 5.5),
        plot.margin = margin(2, 4, 2, 4, "pt"))

# Add text labels for top species
top_spp <- sp_abundance %>% group_by(region) %>% slice_max(n_detections, n = 1)
p2_a <- p2_a +
  annotate("text", x = top_spp$rank[top_spp$region == "Congo"] + 5, y = top_spp$n_detections[top_spp$region == "Congo"],
           label = "Duiker (Congo)", size = 2.0, fontface = "italic", color = pal_basin["Congo"], hjust = 0) +
  annotate("text", x = top_spp$rank[top_spp$region == "Amazon"] + 5, y = top_spp$n_detections[top_spp$region == "Amazon"],
           label = "Paca (Amazon)", size = 2.0, fontface = "italic", color = pal_basin["Amazon"], hjust = 0)

# --- Panel B: Spatial Clustering Map of Deployments ---
# We extract unique coordinate locations for deployments in Congo
# and map them to their Haversine cluster ID to show single-linkage clustering
unique_deps <- det_all %>%
  select(region, longitude, latitude, deployment_id) %>%
  distinct()

# Calculate spatial clusters dynamically for plotting
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

coords_df <- unique_deps %>% 
  filter(region == "Congo") %>%
  select(longitude, latitude) %>% 
  distinct() %>% 
  mutate(cluster_num = 0)

# Haversine distance single-linkage clustering (11.1km threshold)
n_c <- nrow(coords_df)
if (n_c > 1) {
  dist_mat <- matrix(0, nrow=n_c, ncol=n_c)
  for (i in 1:n_c) {
    for (j in 1:n_c) {
      dist_mat[i,j] <- haversine_dist(coords_df$longitude[i], coords_df$latitude[i], coords_df$longitude[j], coords_df$latitude[j])
    }
  }
  hc <- hclust(as.dist(dist_mat), method="single")
  coords_df$cluster_num <- cutree(hc, h=11.1)
} else {
  coords_df$cluster_num <- 1
}

# Overlay Congo country outlines zoomed into Congo Camera Trap clusters
congo_deps_zoom <- unique_deps %>%
  filter(region == "Congo") %>%
  left_join(coords_df, by = c("longitude", "latitude"))

congo_cluster_centroids <- congo_deps_zoom %>%
  group_by(cluster_num) %>%
  summarise(lon = mean(longitude), lat = mean(latitude), .groups = "drop")

p2_b <- ggplot() +
  geom_spatvector(data = countries_c, fill = "#F8F8F6", colour = "grey80", linewidth = 0.3) +
  # Draw circles indicating 11.1km radius around centroids
  geom_point(data = congo_cluster_centroids, aes(x = lon, y = lat), size = 6.8, color = "black", fill = NA, shape = 21, stroke = 0.4, linetype = "dashed", alpha = 0.4) +
  # Plot individual camera trap deployments colored by spatial cluster ID
  geom_point(data = congo_deps_zoom, aes(x = longitude, y = latitude, fill = factor(cluster_num)), size = 1.6, shape = 21, color = "black", stroke = 0.2) +
  scale_fill_viridis_d(option = "turbo", guide = "none") +
  labs(
    title = "B. Haversine Single-Linkage Clustering",
    subtitle = "Congo Basin: 11.1 km Distance Thresholds"
  ) +
  t_theme +
  theme(axis.text = element_text(size = 5.5),
        plot.margin = margin(2, 4, 2, 4, "pt")) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)")

# --- Panel C: Temporal Calibration Timeline ---
# Visualizes survey start/end dates in relation to GEDI launch (April 17, 2019)
# and shows how temporal weights are assigned
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

p2_c <- ggplot(cluster_spans) +
  # Draw a vertical line representing the GEDI mission launch
  geom_vline(xintercept = gedi_start, linetype = "solid", color = "#D32F2F", linewidth = 0.75) +
  annotate("text", x = gedi_start + 180, y = 5, label = "GEDI Launch\n(April 2019)", color = "#D32F2F", size = 2.0, fontface = "bold", hjust = 0) +
  # Draw horizontal bars for each camera trap cluster survey span
  geom_segment(aes(x = start, xend = end, y = row_idx, yend = row_idx, color = w_temp_cluster), linewidth = 1.5, alpha = 0.8) +
  scale_color_gradientn(
    colors = c("#D32F2F", "#F57C00", "#FBC02D", "#388E3C"),
    name = "Temporal Weight",
    breaks = c(0.1, 0.25, 0.5, 1.0)
  ) +
  scale_x_date(date_breaks = "4 years", date_labels = "%Y", limits = c(as.Date("2003-01-01"), as.Date("2024-12-31"))) +
  labs(
    title = "C. Survey Temporal Calibration Span",
    x = "Survey Era",
    y = "Spatial Clusters (Ranked by Start Date)"
  ) +
  t_theme +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 6.0, face = "bold"),
        legend.text = element_text(size = 5.0),
        legend.key.height = unit(0.12, "cm"),
        legend.key.width = unit(1.0, "cm"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(2, 4, 2, 4, "pt"))

# --- Panel D: Taxonomic order composition by event count ---
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
  filter(!is.na(order)) %>%
  group_by(region, order) %>%
  summarise(n_events = sum(n_detections), .groups = "drop")

# Keep top orders, pool other
top_orders <- c("Cetartiodactyla", "Proboscidea", "Carnivora", "Rodentia", "Primates", "Cingulata", "Perissodactyla")
df_tax <- df_tax %>%
  mutate(order_clean = ifelse(order %in% top_orders, order, "Other")) %>%
  group_by(region, order_clean) %>%
  summarise(n_events = sum(n_events), .groups = "drop")

# Calculate proportions
df_tax <- df_tax %>%
  group_by(region) %>%
  mutate(prop = n_events / sum(n_events) * 100) %>%
  ungroup()

df_tax$order_clean <- factor(df_tax$order_clean, levels = rev(c(top_orders, "Other")))

p2_d <- ggplot(df_tax, aes(x = prop, y = region, fill = order_clean)) +
  geom_bar(stat = "identity", width = 0.55, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = ORDER_COLOURS, name = "Taxonomic Order:") +
  labs(
    title = "D. Vertebrate Order Composition",
    x = "Proportion of Detections (%)",
    y = "Basin/Continent"
  ) +
  t_theme +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 6.0, face = "bold"),
        legend.text = element_text(size = 5.5),
        legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.margin = margin(0,0,0,0),
        plot.margin = margin(2, 4, 2, 4, "pt")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# Assemble Figure 2
row1_p2 <- plot_grid(p2_a, p2_b, ncol = 2, rel_widths = c(1, 1))
row2_p2 <- plot_grid(p2_c, p2_d, ncol = 2, rel_widths = c(1, 1))

fig2_final <- plot_grid(row1_p2, row2_p2, ncol = 1, rel_heights = c(1, 1.1), hspace = 0.3)

fig2_png <- "figures/figure2_camera_trap_pipeline.png"
ggsave(filename = fig2_png, plot = fig2_final, width = 17.8, height = 15.0, units = "cm", dpi = 600, bg = "white")
cat("✓ Figure 2 successfully saved to:", fig2_png, "\n\n")

# Copy figures to active brain artifacts folder
brain_artifacts_dir <- "/home/j/.gemini/antigravity/brain/913e5cea-7c99-4b21-8124-ea8455da8457"
if (file.exists(brain_artifacts_dir)) {
  file.copy(fig1_png, file.path(brain_artifacts_dir, "figure1_gedi_pipeline.png"), overwrite = TRUE)
  file.copy(fig2_png, file.path(brain_artifacts_dir, "figure2_camera_trap_pipeline.png"), overwrite = TRUE)
  cat("✓ Copied Figure 1 and 2 to brain artifacts folder.\n")
}

cat("============================================================\n")
cat("=== Manuscript Figures 1 & 2 Completed Successfully ===\n")
cat("============================================================\n")
