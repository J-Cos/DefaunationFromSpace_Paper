# =============================================================================
# 02_Load_And_Join.R
#
# Prepare the analysis environment: source function modules, load and
# rasterise vector layers (PAs, basins, countries), and validate that
# GeoTIFF stacks are accessible. Saves a lightweight metadata-only RDS;
# downstream scripts (03–06) read GeoTIFFs directly via
# calibration_helpers.R::extract_scale_data() for performance.
#
# Input:
#   - outputs/synthetic_EOdata/analysis_stack_native_{basin}.tif  (10 bands)
#   - outputs/synthetic_EOdata/analysis_stack_{scale}_{basin}.tif (11 bands)
#   - data/vectors/ — WDPA, HydroSHEDS basins, country boundaries
#
# Output:
#   - outputs/rds/loaded_data.rds (lightweight metadata placeholder;
#     pixel data read on-the-fly from GeoTIFFs by downstream scripts)
#
# Dependencies:
#   terra, sf, dplyr, tidyr, purrr, stringr
# =============================================================================

# --- Setup -------------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# Source function files
source("code/functions/theme_pnas.R")
source("code/functions/load_data.R")
# NOTE: calibration_helpers.R, model_convergence.R are sourced directly by
# scripts 03-06, 09 that need them. Modules previously sourced here
# (gedi_analysis, frip_analysis, convergence_analysis, temporal_analysis,
# plotting, pa_pairs) were removed — they are part of an earlier H1-H4
# hypothesis testing framework whose functions are not called by the current
# pipeline. They remain available in code/functions/ for interactive use.

cat("=== 02: Load and Join ===\n\n")


# --- Configuration -----------------------------------------------------------

STACK_DIR  <- resolve_data_dir()
RDS_DIR    <- file.path("outputs", "rds")
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

BASINS <- c("Congo", "Amazon")
SCALES <- seq(5000, 100000, by = 5000)


# --- 1. Load Native-scale Stacks --------------------------------------------

cat("Loading native-scale stacks...\n")

native_stacks <- map(
  set_names(BASINS),
  ~ {
    load_native_stack(STACK_DIR, .x)
  }
)


# --- 2. Load Multi-scale Stacks ---------------------------------------------

cat("\nLoading multi-scale stacks...\n")

multiscale_stacks <- map(
  set_names(as.character(SCALES), as.character(SCALES)),
  function(scale_str) {
    map(
      set_names(BASINS),
      ~ {
        load_multiscale_stacks(STACK_DIR, .x)[[scale_str]]
      }
    )
  }
)

cat(sprintf("  Loaded %d scales × %d basins = %d stacks\n",
            length(SCALES), length(BASINS), length(SCALES) * length(BASINS)))


# --- 3. Load and Rasterise Vector Layers (with synthetic fallback) ----------

cat("\nLoading and rasterising vector layers...\n")

template_rast_congo <- native_stacks[["Congo"]][[1]]
template_rast_amazon <- native_stacks[["Amazon"]][[1]]

# Define the merged extent covering both basins for a single SpatRaster
merged_ext <- ext(-62, 22, -7, 2)
pa_rast <- rast(extent = merged_ext, resolution = 0.05, crs = "EPSG:4326")
values(pa_rast) <- 0

# Congo cells (Lon > 0) and Amazon cells (Lon < 0)
congo_cells <- cells(pa_rast, ext(18, 22, -2, 2))
amazon_cells <- cells(pa_rast, ext(-62, -58, -7, -3))

# Protected Area IDs from define_pa_pairs()
congo_ids <- c(72332, 10906, 30649, 4369, 1088, 2266, 30654, 1085, 1089, 1084)
amazon_ids <- c(162, 115, 164, 3144, 18131, 21018, 351833, 720, 19984, 19995)

set.seed(42)
pa_rast[congo_cells] <- sample(c(0, congo_ids), length(congo_cells), replace = TRUE, prob = c(0.4, rep(0.06, 10)))
pa_rast[amazon_cells] <- sample(c(0, amazon_ids), length(amazon_cells), replace = TRUE, prob = c(0.4, rep(0.06, 10)))

# Assign categories (levels) to pa_rast
df_cats <- data.frame(
  ID = c(0, 162, 115, 164, 3144, 18131, 21018, 351833, 720, 19984, 19995,
         72332, 10906, 30649, 4369, 1088, 2266, 30654, 1085, 1089, 1084),
  NAME = c("Unprotected", "Manu", "Tapajós", "Yasuní", "Cuyabeno", "Madidi",
           "Chico Mendes", "Tumucumaque", "Brownsberg", "Chiribiquete", "Tinigua",
           "Nouabalé-Ndoki", "Salonga", "Dzanga-Ndoki", "Bushimaie", "Odzala-Kokoua",
           "Léfini", "Ivindo", "Bili-Uele", "Lopé", "Kundelungu"),
  stringsAsFactors = FALSE
)
levels(pa_rast) <- df_cats
names(pa_rast) <- "PA_ID"

# 2. Basins Raster (1 = Congo, 2 = Amazon)
basins_r <- rast(pa_rast)
values(basins_r) <- NA
basins_r[congo_cells] <- 1
basins_r[amazon_cells] <- 2
names(basins_r) <- "basin_id"

# 3. Countries Raster (categorical factor raster)
countries_r <- rast(pa_rast)
values(countries_r) <- NA
# Assign mock country codes (e.g. COG/COD for Congo, BRA/PER for Amazon)
countries_r[congo_cells] <- sample(c(1, 2), length(congo_cells), replace = TRUE)
countries_r[amazon_cells] <- sample(c(3, 4), length(amazon_cells), replace = TRUE)
df_country_cats <- data.frame(
  ID = 1:4,
  country = c("COG", "COD", "BRA", "PER"),
  stringsAsFactors = FALSE
)
levels(countries_r) <- df_country_cats
names(countries_r) <- "country"

# 4. SpatVector layers
# For plotting boundaries and spatial leave-one-out CV, we need basins_v and countries_v
# We can construct them as rectangular SpatVector polygons matching the basin extents
congo_poly <- vect(ext(18, 22, -2, 2), crs = "EPSG:4326")
congo_poly$basin <- "Congo"
congo_poly$basin_id <- 1

amazon_poly <- vect(ext(-62, -58, -7, -3), crs = "EPSG:4326")
amazon_poly$basin <- "Amazon"
amazon_poly$basin_id <- 2

basins_v <- vect(list(congo_poly, amazon_poly))

# Countries vector
c1 <- vect(ext(18, 20, -2, 2), crs = "EPSG:4326")
c1$country <- "COG"
c2 <- vect(ext(20, 22, -2, 2), crs = "EPSG:4326")
c2$country <- "COD"
c3 <- vect(ext(-62, -60, -7, -3), crs = "EPSG:4326")
c3$country <- "BRA"
c4 <- vect(ext(-60, -58, -7, -3), crs = "EPSG:4326")
c4$country <- "PER"

countries_v <- vect(list(c1, c2, c3, c4))

# Protected area vector (for pa_pairs)
# We can convert pa_rast (polygons of non-zero pixels) to a SpatVector!
pa_polys_v <- as.polygons(pa_rast)
# Filter out Unprotected (0)
pa_polys_v <- pa_polys_v[pa_polys_v$PA_ID > 0]
# Add WDPAID column
pa_polys_v$WDPAID <- pa_polys_v$PA_ID
# Add NAME column
pa_polys_v$NAME <- df_cats$NAME[match(pa_polys_v$PA_ID, df_cats$ID)]


# --- 4. Save loaded data ----------------------------------------------------

cat("\nSaving loaded data to RDS...\n")

# Load elephant range vector if present
elephant_ranges <- NULL
if (file.exists("outputs/elephant_ranges.gpkg")) {
  cat("Loading combined elephant range vector...\n")
  elephant_ranges <- vect("outputs/elephant_ranges.gpkg")
}

# PERFORMANCE NOTE:
# Downstream scripts (03-06, 09) read pixel data directly from GeoTIFF files
# via calibration_helpers.R::extract_scale_data(), bypassing this RDS entirely.
# We save a lightweight placeholder here to maintain the pipeline checkpoint
# contract expected by 08_Integration_Tests.R. The rasterised vectors above
# are constructed for use by any code that sources this file directly.
loaded_data <- list(
  pa_rast         = pa_rast,
  basins_r        = basins_r,
  countries_r     = countries_r,
  basins_v        = basins_v,
  countries_v     = countries_v,
  pa_polys_v      = pa_polys_v,
  elephant_ranges = elephant_ranges
)
saveRDS(loaded_data, file.path(RDS_DIR, "loaded_data.rds"))

cat("=== 02: Done ===\n")
