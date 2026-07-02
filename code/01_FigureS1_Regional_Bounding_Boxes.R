# =============================================================================
# code/01_FigureS1_Regional_Bounding_Boxes.R
#
# Beautiful tidyterra map script producing a column of 3 maps (Amazon, Congo, 
# Southeast Asia) running exactly from 15°S to 15°N latitude and spanning 
# exactly 50° longitude (centered on their respective regional centroids), 
# showing country outlines, highlighting their bounding boxes, and overlaying 
# the geographical range maps for all wild proboscid species, categorized 
# by IUCN status (Extant, Possibly Extant, and Possibly Extinct).
#
# Outputs:
#   - figures/figureS1.png (and .pdf)
#   - outputs/elephant_ranges.gpkg
# =============================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(cowplot)

cat("=== Generating Figure S1: Regional Bounding Boxes Map with Elephant Range Overlays ===\n\n")

# Ensure directories exist
dir.create("outputs", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# 1. Load Country Boundaries (uses the local world-administrative-boundaries)
global_countries <- vect("data/world-administrative-boundaries")

# 2. Load Proboscidea range maps efficiently using SQL query
prob_path <- "data/MAMMALS_TERRESTRIAL_ONLY(1)/MAMMALS_TERRESTRIAL_ONLY.shp"
proboscids <- NULL
if (file.exists(prob_path)) {
  cat("Loading Proboscidea range maps...\n")
  proboscids <- terra::vect(prob_path, query = "SELECT * FROM MAMMALS_TERRESTRIAL_ONLY WHERE order_ = 'PROBOSCIDEA'")
  
  # Categorise range status cleanly in base R (100% robust)
  status_vec <- rep("Unknown", nrow(proboscids))
  status_vec[proboscids$presence == 1] <- "Extant"
  status_vec[proboscids$presence == 3] <- "Possibly Extant"
  status_vec[proboscids$presence == 4] <- "Possibly Extinct"
  proboscids$status <- factor(status_vec, levels = c("Extant", "Possibly Extant", "Possibly Extinct"))
} else {
  stop("Error: Mammal shapefile not found at data/MAMMALS_TERRESTRIAL_ONLY(1)/\n")
}

# 3. Define the 3 Bounding Boxes: Exactly 50° Wide × 30° High (-15 to 15 Lat)
# Amazon: Centered at Lon -60, spanning -85 to -35
ext_amazon <- ext(-85, -35, -15, 15)
# Congo: Centered at Lon 20, spanning -5 to 45
ext_congo  <- ext(-5, 45, -15, 15)
# SE Asia: Centered at Lon 115, spanning 90 to 140
ext_sea    <- ext(90, 140, -15, 15)

# Convert extents to SpatVector polygons to plot them as bounding box overlays
bbox_amazon_v <- as.polygons(ext_amazon, crs = "EPSG:4326")
bbox_congo_v  <- as.polygons(ext_congo, crs = "EPSG:4326")
bbox_sea_v    <- as.polygons(ext_sea, crs = "EPSG:4326")

# Crop country boundaries to each region's bounding box
countries_amazon <- crop(global_countries, ext_amazon)
countries_congo  <- crop(global_countries, ext_congo)
countries_sea    <- crop(global_countries, ext_sea)

# Crop Proboscidea range maps to respective basins
prob_congo <- NULL
prob_sea <- NULL
if (!is.null(proboscids)) {
  # Congo: Keep Loxodonta africana and Loxodonta cyclotis
  prob_congo_raw <- crop(proboscids, ext_congo)
  prob_congo <- prob_congo_raw[prob_congo_raw$sci_name %in% c("Loxodonta africana", "Loxodonta cyclotis"), ]
  
  # SE Asia: Keep Elephas maximus
  prob_sea_raw <- crop(proboscids, ext_sea)
  prob_sea <- prob_sea_raw[prob_sea_raw$sci_name == "Elephas maximus", ]
}

# Combine and export elephant range SpatVector to GeoPackage
combined_proboscids <- NULL
if (!is.null(prob_congo) && !is.null(prob_sea)) {
  combined_proboscids <- rbind(prob_congo, prob_sea)
} else if (!is.null(prob_congo)) {
  combined_proboscids <- prob_congo
} else if (!is.null(prob_sea)) {
  combined_proboscids <- prob_sea
}

if (!is.null(combined_proboscids)) {
  gpkg_path <- "outputs/elephant_ranges.gpkg"
  writeVector(combined_proboscids, gpkg_path, overwrite = TRUE)
  cat("✓ Successfully saved combined elephant range vector to outputs/elephant_ranges.gpkg\n")
}

# Define a premium, clean PNAS-style map theme
map_theme <- theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(size = 10, face = "bold", margin = margin(b = 4)),
    panel.background = element_rect(fill = "#EBF5FB", color = NA), # Soft light blue ocean fill
    panel.grid.major = element_line(color = "white", linewidth = 0.2), # Soft white grid lines
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 6, color = "grey50"),
    plot.margin = margin(2, 2, 2, 2, "pt"),
    panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5)
  )

# Plot Panel A: Amazon
p_amazon <- ggplot() +
  geom_spatvector(data = countries_amazon, fill = "#F4F6F7", color = "grey70", linewidth = 0.25) +
  # Draw a bold bounding box outline inside the plot
  geom_spatvector(data = bbox_amazon_v, fill = NA, color = "#E65100", linewidth = 1.0) +
  coord_sf(xlim = c(-85, -35), ylim = c(-15, 15), expand = FALSE) +
  labs(title = "A") +
  map_theme +
  # Amazon has no wild proboscids; place a small text overlay acknowledging this
  annotate("text", x = -60, y = -12, label = "Neotropics: Megafauna Depleted (No Proboscids)", 
           fontface = "italic", size = 2.4, color = "grey40")

# Plot Panel B: Congo
p_congo <- ggplot() +
  geom_spatvector(data = countries_congo, fill = "#F4F6F7", color = "grey70", linewidth = 0.25)

# Overlay Proboscids in Congo if present, colored by status
if (!is.null(prob_congo) && nrow(prob_congo) > 0) {
  p_congo <- p_congo +
    geom_spatvector(data = prob_congo, aes(fill = status, color = status), linewidth = 0.3, alpha = 0.22) +
    scale_fill_manual(
      values = c("Extant" = "#26A69A", "Possibly Extant" = "#FFA726", "Possibly Extinct" = "#EF5350"),
      name = "Elephant Status:",
      drop = FALSE
    ) +
    scale_color_manual(
      values = c("Extant" = "#00695C", "Possibly Extant" = "#EF6C00", "Possibly Extinct" = "#C62828"),
      name = "Elephant Status:",
      drop = FALSE
    ) +
    theme(
      legend.position = c(0.18, 0.24),
      legend.title = element_text(size = 5.5, face = "bold"),
      legend.text = element_text(size = 5.0),
      legend.background = element_rect(fill = alpha("white", 0.8), color = "grey80", linewidth = 0.2),
      legend.key.size = unit(0.2, "cm"),
      legend.margin = margin(2, 2, 2, 2, "pt")
    )
}

p_congo <- p_congo +
  # Draw a bold bounding box outline inside the plot
  geom_spatvector(data = bbox_congo_v, fill = NA, color = "#1B5E20", linewidth = 1.0) +
  coord_sf(xlim = c(-5, 45), ylim = c(-15, 15), expand = FALSE) +
  labs(title = "B") +
  map_theme

# Plot Panel C: Southeast Asia
p_sea <- ggplot() +
  geom_spatvector(data = countries_sea, fill = "#F4F6F7", color = "grey70", linewidth = 0.25)

# Overlay Proboscids in SE Asia if present, colored by status
if (!is.null(prob_sea) && nrow(prob_sea) > 0) {
  p_sea <- p_sea +
    geom_spatvector(data = prob_sea, aes(fill = status, color = status), linewidth = 0.3, alpha = 0.22) +
    scale_fill_manual(
      values = c("Extant" = "#26A69A", "Possibly Extant" = "#FFA726", "Possibly Extinct" = "#EF5350"),
      name = "Elephant Status:",
      drop = FALSE
    ) +
    scale_color_manual(
      values = c("Extant" = "#00695C", "Possibly Extant" = "#EF6C00", "Possibly Extinct" = "#C62828"),
      name = "Elephant Status:",
      drop = FALSE
    ) +
    theme(
      legend.position = c(0.18, 0.24),
      legend.title = element_text(size = 5.5, face = "bold"),
      legend.text = element_text(size = 5.0),
      legend.background = element_rect(fill = alpha("white", 0.8), color = "grey80", linewidth = 0.2),
      legend.key.size = unit(0.2, "cm"),
      legend.margin = margin(2, 2, 2, 2, "pt")
    )
}

p_sea <- p_sea +
  # Draw a bold bounding box outline inside the plot
  geom_spatvector(data = bbox_sea_v, fill = NA, color = "#0D47A1", linewidth = 1.0) +
  coord_sf(xlim = c(90, 140), ylim = c(-15, 15), expand = FALSE) +
  labs(title = "C") +
  map_theme

# Combine into a single vertical column of 3 figures
col_figure <- plot_grid(
  p_amazon,
  p_congo,
  p_sea,
  ncol = 1,
  align = "v"
)

# Save high-resolution outputs
fig_path_png_local <- "figures/figureS1.png"
fig_path_pdf_local <- "figures/figureS1.pdf"

# Save high-resolution publication-quality PNG and PDF (300 DPI)
ggsave(fig_path_png_local, plot = col_figure, width = 11.5, height = 16.5, units = "cm", dpi = 300, bg = "white")
ggsave(fig_path_pdf_local, plot = col_figure, width = 11.5, height = 16.5, units = "cm", dpi = 300, bg = "white")

cat("✓ Successfully saved Figure S1 regional bounding box column figure to:\n  -", fig_path_png_local, "\n  -", fig_path_pdf_local, "\n")
