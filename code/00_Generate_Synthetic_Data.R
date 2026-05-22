# =============================================================================
# 00_Generate_Synthetic_Data.R
#
# Generates synthetic multi-band GeoTIFFs that mimic the output of
# 03_FRIP_Signals_And_Drive_Exports_GEE.ipynb (Notebook 3).
#
# Creates data at MODIS native scale (~463m) for two small 200x200km
# basin subsets, then aggregates to each of 20 multi-scale resolutions
# (5km to 100km in 5km steps).
#
# Output naming convention matches NB3:
#   data/analysis_stack_native_{basin}.tif   (10 bands)
#   data/analysis_stack_{scale}_{basin}.tif  (11 bands)
#
# This allows the full R analysis pipeline to be developed and tested
# before the real GEE exports complete.
# =============================================================================

library(terra)

set.seed(42)
cat("=== Synthetic Data Generator for Defaunation Pipeline ===\n\n")

# -----------------------------------------------------------------------------
# 1. CONFIGURATION
# -----------------------------------------------------------------------------

# MODIS equatorial pixel size in degrees (~463.3m at equator)
MODIS_SCALE_M  <- 463.3127165279165
MODIS_SCALE_DEG <- MODIS_SCALE_M / 111320  # ~0.00416°

# Multi-scale target resolutions (meters)
SCALES <- seq(5000, 100000, by = 5000)

# Full-scale basin bounding boxes
# Congo basin extent
CONGO_BBOX  <- ext(8, 35, -12, 8)
# Amazon basin extent
AMAZON_BBOX <- ext(-73, -44, -18, 8)

BASINS <- list(
  Congo  = CONGO_BBOX,
  Amazon = AMAZON_BBOX
)

# Band definitions (must match NB3 output exactly)
MULTISCALE_BANDS <- c(
  "frip", "frip_mk_tau",
  "uoi", "rh98", "gedi_n",
  "elevation", "slope", "hnd", "precip", "clay", "forest_fraction"
)

NATIVE_BANDS <- c(
  "uoi", "rh98", "gedi_n",
  "elevation", "slope", "hnd", "precip", "clay", "forest_fraction",
  "Npp_median"
)

# Output directory (real NB3 exports go to outputs/EOdata/)
OUT_DIR <- file.path("outputs", "synthetic_EOdata")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("  MODIS pixel: %.4f° (~%.1fm)\n", MODIS_SCALE_DEG, MODIS_SCALE_M))
cat(sprintf("  Scales: %d–%dkm (%d scales)\n",
            SCALES[1] / 1000, SCALES[length(SCALES)] / 1000, length(SCALES)))
cat(sprintf("  Basins: %s\n", paste(names(BASINS), collapse = ", ")))
cat(sprintf("  Output: %s/\n\n", OUT_DIR))


# -----------------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# -----------------------------------------------------------------------------

make_template <- function(bbox, res_deg) {
  #' Create an empty single-band raster template at a given resolution.
  rast(extent = bbox, resolution = res_deg, crs = "EPSG:4326")
}

add_spatial_structure <- function(r, window_size = 15) {
  #' Add spatial autocorrelation to white noise by applying a focal mean.
  #' This creates ecologically plausible spatial patterns.
  w <- matrix(1, nrow = window_size, ncol = window_size)
  focal(r, w = w, fun = "mean", na.rm = TRUE, pad = TRUE)
}

generate_band <- function(template, band_name, basin_name) {
  #' Generate a single band of spatially autocorrelated synthetic data.
  #'

  #' Values are calibrated to ecologically plausible ranges.
  #' Inter-basin offsets encode the study hypotheses:
  #'   H1: Congo > Amazon for UOI (intact structure)
  #'   H2: Amazon > Congo for FRIP (broken nutrient pump)

  r <- template
  n <- ncell(r)
  is_congo <- basin_name == "Congo"

  # Generate spatially correlated base noise
  values(r) <- rnorm(n)
  r <- add_spatial_structure(r)

  # Normalise to [0, 1] range
  v <- values(r, mat = FALSE)
  v <- (v - min(v, na.rm = TRUE)) / (max(v, na.rm = TRUE) - min(v, na.rm = TRUE))

  # Map to ecologically plausible ranges with inter-basin shifts
  values(r) <- switch(band_name,
    # --- FRIP signals (H2: Amazon > Congo) ---
    frip = {
      mu <- ifelse(is_congo, -0.15, 0.10)
      qnorm(pmin(0.999, pmax(0.001, v)), mean = mu, sd = 0.20)
    },
    frip_mk_tau = {
      mu <- ifelse(is_congo, -0.02, 0.08)
      qnorm(pmin(0.999, pmax(0.001, v)), mean = mu, sd = 0.12)
    },

    # --- GEDI structural signals (H1: Congo > Amazon for UOI) ---
    uoi = {
      mu <- ifelse(is_congo, 0.65, 0.45)
      pmin(1, pmax(0, qnorm(pmin(0.999, pmax(0.001, v)), mean = mu, sd = 0.12)))
    },
    rh98 = {
      mu <- ifelse(is_congo, 35, 28)
      pmax(5, qnorm(pmin(0.999, pmax(0.001, v)), mean = mu, sd = 7))
    },
    gedi_n = {
      pmax(1, round(qnorm(pmin(0.999, pmax(0.001, v)), mean = 80, sd = 40)))
    },

    # --- Environmental covariates (broadly similar between basins) ---
    elevation = pmax(0, qnorm(pmin(0.999, pmax(0.001, v)), mean = 250, sd = 120)),
    slope     = pmax(0, qnorm(pmin(0.999, pmax(0.001, v)), mean = 3, sd = 2)),
    hnd       = pmax(0, qnorm(pmin(0.999, pmax(0.001, v)), mean = 15, sd = 10)),
    precip    = pmax(500, qnorm(pmin(0.999, pmax(0.001, v)), mean = 2000, sd = 400)),
    clay      = pmin(80, pmax(5, qnorm(pmin(0.999, pmax(0.001, v)), mean = 35, sd = 12))),
    forest_fraction = pmin(1, pmax(0, qnorm(pmin(0.999, pmax(0.001, v)), mean = 0.88, sd = 0.08))),

    # --- NPP (native-scale only) ---
    Npp_median = pmax(0, qnorm(pmin(0.999, pmax(0.001, v)), mean = 12000, sd = 3000)),

    # Fallback
    v
  )

  return(r)
}

aggregate_stack <- function(stack, target_scale_m, band_names) {
  #' Aggregate a MODIS-scale stack to a coarser target scale.
  #'
  #' Uses mean for all bands except gedi_n which uses sum,
  #' matching the GEE reduceResolution reducers in NB3.

  agg_factor <- max(2, round(target_scale_m / MODIS_SCALE_M))

  result_bands <- list()
  for (i in seq_along(band_names)) {
    bname <- band_names[i]
    b <- stack[[i]]

    if (bname == "gedi_n") {
      result_bands[[i]] <- aggregate(b, fact = agg_factor, fun = "sum", na.rm = TRUE)
    } else {
      result_bands[[i]] <- aggregate(b, fact = agg_factor, fun = "mean", na.rm = TRUE)
    }
  }

  out <- rast(result_bands)
  names(out) <- band_names
  return(out)
}


# -----------------------------------------------------------------------------
# 3. GENERATE AND EXPORT
# -----------------------------------------------------------------------------

for (basin_name in names(BASINS)) {
  bbox <- BASINS[[basin_name]]
  cat(sprintf("--- %s basin ---\n", basin_name))

  # Create MODIS-scale template
  template <- make_template(bbox, MODIS_SCALE_DEG)
  dims <- dim(template)
  cat(sprintf("  Template: %d x %d pixels at %.4f° (~%.0fm)\n",
              dims[1], dims[2], MODIS_SCALE_DEG, MODIS_SCALE_M))

  # Generate all 12 unique bands at MODIS scale
  all_band_names <- unique(c(MULTISCALE_BANDS, NATIVE_BANDS))
  all_bands <- list()
  for (bname in all_band_names) {
    all_bands[[bname]] <- generate_band(template, bname, basin_name)
  }
  full_stack <- rast(all_bands)
  names(full_stack) <- all_band_names

  # --- Export native-scale stack (10 bands) ---
  native_stack <- full_stack[[NATIVE_BANDS]]
  native_path <- file.path(OUT_DIR, sprintf("analysis_stack_native_%s.tif", basin_name))
  writeRaster(native_stack, native_path, overwrite = TRUE)
  cat(sprintf("  ✓ Native stack: %d bands → %s\n", nlyr(native_stack), basename(native_path)))

  # --- Export multi-scale stacks (11 bands at each of 20 scales) ---
  multiscale_source <- full_stack[[MULTISCALE_BANDS]]

  for (scale in SCALES) {
    agg_stack <- aggregate_stack(multiscale_source, scale, MULTISCALE_BANDS)
    agg_dims <- dim(agg_stack)
    out_path <- file.path(OUT_DIR, sprintf("analysis_stack_%d_%s.tif", scale, basin_name))
    writeRaster(agg_stack, out_path, overwrite = TRUE)
    cat(sprintf("  ✓ %3dkm: %d bands, %dx%d px → %s\n",
                scale / 1000, nlyr(agg_stack), agg_dims[1], agg_dims[2], basename(out_path)))
  }

  cat("\n")
}

# --- Summary ---
tif_files <- list.files(OUT_DIR, pattern = "\\.tif$")
cat(sprintf("=== Done! %d GeoTIFFs written to %s/ ===\n", length(tif_files), OUT_DIR))
cat("Files:\n")
for (f in sort(tif_files)) {
  sz <- file.size(file.path(OUT_DIR, f))
  cat(sprintf("  %s (%.1f KB)\n", f, sz / 1024))
}
