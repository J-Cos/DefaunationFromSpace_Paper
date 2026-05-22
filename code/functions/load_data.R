# =============================================================================
# load_data.R
#
# Data loading functions for the Defaunation-from-Space analysis pipeline.
#
# Loads multi-scale and native-resolution GeoTIFF stacks exported by
# GEE Notebook 3, basin/country/PA vector layers, and external defaunation
# indices. All paths and constants are defined here for single-source
# configuration.
# =============================================================================

library(terra)
library(dplyr)
library(stringr)

# =============================================================================
# PATH & SCALE CONSTANTS
# =============================================================================

#' Root directory for analysis-ready rasters (GEE exports or synthetic data)
#' @export
DATA_DIR <- file.path("outputs", "synthetic_EOdata")

#' Root directory for legacy FRIP-era data (DefaunationFromSpace_Paper)
#' @export
LEGACY_DATA_DIR <- file.path("legacy", "DefaunationFromSpace_Paper", "Data")

#' Multi-scale target resolutions in metres (5 km to 100 km, 5 km steps)
#' @export
SCALES <- seq(5000, 100000, by = 5000)

#' Basin names used throughout the pipeline
#' @export
BASINS <- c("Congo", "Amazon")

#' Minimum cover fraction for rasterise-and-mask (pixels must be ≥ 99%
#' within a single polygon to be retained; see legacy Functions.r)
#' @export
COVER_THRESHOLD <- 0.99

#' Minimum PA area (km²) for focal Amazon PAs
#' @export
PA_MIN_AREA_AMAZON_KM2 <- 40000

#' Minimum PA area (km²) for focal Congo PAs
#' @export
PA_MIN_AREA_CONGO_KM2 <- 20000


# =============================================================================
# MULTI-SCALE STACKS
# =============================================================================

#' Load multi-scale analysis stacks for a single basin
#'
#' Reads all \code{analysis_stack_{scale}_{basin}.tif} files from
#' \code{data_dir} and returns a named list of \code{SpatRaster} objects
#' keyed by scale in metres (e.g. \code{"5000"}, \code{"10000"}, …).
#'
#' @param data_dir Character. Path to the directory containing GeoTIFFs.
#' @param basin Character. Basin name (\code{"Congo"} or \code{"Amazon"}).
#'
#' @return A named list of 20 \code{terra::SpatRaster} objects, one per scale.
#'
#' @details File pattern: \code{analysis_stack_\\d+_{basin}.tif}.
#'   Band names are read from the file and assumed to match
#'   \code{c("frip", "frip_mk_tau", "uoi", "rh98", "gedi_n",
#'   "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")}.
#'
#' @export
load_multiscale_stacks <- function(data_dir = DATA_DIR, basin) {
  stacks <- list()
  for (scale in SCALES) {
    filename <- file.path(data_dir, sprintf("analysis_stack_%d_%s.tif", scale, basin))
    if (!file.exists(filename)) {
      stop("File does not exist: ", filename)
    }
    r <- rast(filename)
    names(r) <- c("frip", "frip_mk_tau", "uoi", "rh98", "gedi_n",
                  "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
    stacks[[as.character(scale)]] <- r
  }
  return(stacks)
}


#' Load the native-resolution analysis stack for a single basin
#'
#' Reads \code{analysis_stack_native_{basin}.tif} from \code{data_dir}.
#'
#' @param data_dir Character. Path to the directory containing GeoTIFFs.
#' @param basin Character. Basin name (\code{"Congo"} or \code{"Amazon"}).
#'
#' @return A \code{terra::SpatRaster} with bands:
#'   \code{c("uoi", "rh98", "gedi_n", "elevation", "slope", "hnd",
#'   "precip", "clay", "forest_fraction", "Npp_median")}.
#'
#' @export
load_native_stack <- function(data_dir = DATA_DIR, basin) {
  filename <- file.path(data_dir, sprintf("analysis_stack_native_%s.tif", basin))
  if (!file.exists(filename)) {
    stop("File does not exist: ", filename)
  }
  r <- rast(filename)
  names(r) <- c("uoi", "rh98", "gedi_n", "elevation", "slope", "hnd",
                "precip", "clay", "forest_fraction", "Npp_median")
  return(r)
}


#' Load all raster data for both basins
#'
#' Orchestrates \code{load_multiscale_stacks()} and
#' \code{load_native_stack()} for Congo and Amazon, returning a nested list.
#'
#' @param data_dir Character. Path to the directory containing GeoTIFFs.
#'
#' @return A named list with elements \code{Congo} and \code{Amazon}, each
#'   containing \code{$multiscale} (list of 20 SpatRasters) and
#'   \code{$native} (single SpatRaster).
#'
#' @export
load_all_data <- function(data_dir = DATA_DIR) {
  res <- list()
  for (basin in BASINS) {
    res[[basin]] <- list(
      multiscale = load_multiscale_stacks(data_dir, basin),
      native = load_native_stack(data_dir, basin)
    )
  }
  return(res)
}


# =============================================================================
# VECTOR LAYERS
# =============================================================================

#' Load HydroSHEDS Level-2 basin boundaries
#'
#' Reads HydroSHEDS Level-2 shapefiles for Africa and South America,
#' selects the Congo basin (\code{HYBAS_ID == 1020018110}) and the
#' Amazon basin (largest \code{SUB_AREA}), and returns them as a
#' merged \code{SpatVector}.
#'
#' @return A \code{terra::SpatVector} with two features (Congo, Amazon).
#'
#' @details Mirrors legacy \code{01_ReprojectAndCropData.r} lines 11-16.
#'
#' @export
load_basins <- function() {
  af_path <- file.path("data", "hybas_af_lev02_v1c")
  sa_path <- file.path("data", "hybas_sa_lev02_v1c")
  
  if (!file.exists(af_path) && !file.exists(file.path(af_path, "hybas_af_lev02_v1c.shp"))) {
    stop("Congo basin HydroSHEDS file missing.")
  }
  if (!file.exists(sa_path) && !file.exists(file.path(sa_path, "hybas_sa_lev02_v1c.shp"))) {
    stop("Amazon basin HydroSHEDS file missing.")
  }
  
  af <- vect(af_path)
  congo <- af[af$HYBAS_ID == 1020018110]
  
  sa <- vect(sa_path)
  amazon <- sa[sa$SUB_AREA == max(sa$SUB_AREA)]
  
  basins <- vect(list(congo, amazon))
  basins$basin <- c("Congo", "Amazon")
  return(basins)
}


#' Load country boundaries cropped to study basins
#'
#' Reads world administrative boundaries, reprojects to the analysis CRS,
#' and crops to each basin extent.
#'
#' @param basins A \code{terra::SpatVector} with two basin features
#'   (as returned by \code{load_basins()}).
#'
#' @return A \code{terra::SpatVector} of country polygons covering
#'   both basins.
#'
#' @details Mirrors legacy \code{01_ReprojectAndCropData.r} lines 29-35.
#'
#' @export
load_countries <- function(basins) {
  world_admin_path <- file.path("data", "world-administrative-boundaries")
  if (!file.exists(world_admin_path)) {
    stop("World administrative boundaries file missing.")
  }
  
  countries_raw <- vect(world_admin_path)
  countries_proj <- project(countries_raw, crs(basins))
  
  countries_congo <- crop(countries_proj, basins[basins$basin == "Congo"])
  countries_amazon <- crop(countries_proj, basins[basins$basin == "Amazon"])
  
  countries <- vect(list(countries_congo, countries_amazon))
  return(countries)
}


#' Load WDPA protected areas cropped and filtered by size
#'
#' Reads WDPA shapefiles, merges parts, reprojects, crops to basin
#' extents, and filters to PAs exceeding area thresholds
#' (\code{PA_MIN_AREA_AMAZON_KM2} for Amazon,
#' \code{PA_MIN_AREA_CONGO_KM2} for Congo).
#'
#' @param basins A \code{terra::SpatVector} with two basin features
#'   (as returned by \code{load_basins()}).
#'
#' @return A \code{terra::SpatVector} of merged, filtered protected areas.
#'
#' @details Mirrors legacy \code{06_PA_analysis_Multiscale.r} lines 20-28.
#'
#' @export
load_pas <- function(basins) {
  pa0_path <- file.path("data", "WDPA_Nov2024_Public_shp_0")
  pa1_path <- file.path("data", "WDPA_Nov2024_Public_shp_1")
  pa2_path <- file.path("data", "WDPA_Nov2024_Public_shp_2")
  
  if (!file.exists(pa0_path) || !file.exists(pa1_path) || !file.exists(pa2_path)) {
    stop("WDPA protected area shapefiles missing.")
  }
  
  pa0 <- vect(pa0_path)
  pa1 <- vect(pa1_path)
  pa2 <- vect(pa2_path)
  pas_raw <- c(pa0, pa1, pa2)
  
  pas_proj <- project(pas_raw, crs(basins))
  
  pas_congo <- crop(pas_proj, basins[basins$basin == "Congo"])
  pas_amazon <- crop(pas_proj, basins[basins$basin == "Amazon"])
  
  pas_congo_filtered <- pas_congo[expanse(pas_congo, unit = "km") > PA_MIN_AREA_CONGO_KM2]
  pas_amazon_filtered <- pas_amazon[expanse(pas_amazon, unit = "km") > PA_MIN_AREA_AMAZON_KM2]
  
  pas_all <- vect(list(pas_congo_filtered, pas_amazon_filtered))
  return(pas_all)
}


# =============================================================================
# RASTERISE & MASK
# =============================================================================

#' Rasterize a vector layer and mask by coverage fraction
#'
#' Converts a \code{SpatVector} to a raster using \code{terra::rasterize},
#' then masks pixels where the maximum single-polygon cover fraction is
#' below \code{COVER_THRESHOLD} (default 0.99). This ensures only pixels
#' unambiguously within a single polygon are retained.
#'
#' @param v A \code{terra::SpatVector} to rasterize.
#' @param rast_template A \code{terra::SpatRaster} defining the target grid.
#' @param col Character. Column name in \code{v} used as the raster field.
#'
#' @return A \code{terra::SpatRaster} with values from \code{col}, masked
#'   to high-coverage pixels.
#'
#' @details Ported from legacy \code{Functions.r::rasteriseAndMask()}.
#'   Formula: \code{any(cover_by_polygon > 0.99) == TRUE} to retain pixel.
#'
#' @export
rasterise_and_mask <- function(v, rast_template, col) {
  r_unmasked <- terra::rasterize(v, rast_template, field = col)
  r_cover <- terra::rasterize(v, rast_template, cover = TRUE, by = col)
  r_coverMask <- any(r_cover > COVER_THRESHOLD, na.rm = TRUE)
  r <- mask(r_unmasked, r_coverMask, maskvalues = FALSE)
  return(r)
}


# =============================================================================
# DEFAUNATION INDICES
# =============================================================================

#' Load and reproject external defaunation indices
#'
#' Reads three raster-format defaunation indices from the legacy data
#' directory, reprojects each to the analysis grid CRS and extent, and
#' returns them as bands in a single \code{SpatRaster}.
#'
#' Indices loaded:
#' \describe{
#'   \item{DefInd}{Defaunation Index (Mollweide projection → analysis CRS)}
#'   \item{Bogoni}{Bogoni et al. defaunation index}
#'   \item{largeDefInd}{Large-bodied defaunation index (Mollweide → analysis CRS)}
#' }
#'
#' @param rast_template A \code{terra::SpatRaster} defining the target grid
#'   (CRS, extent, and resolution). Typically one scale from the
#'   multi-scale stack.
#'
#' @return A \code{terra::SpatRaster} with three bands:
#'   \code{c("DI", "DI_Bogoni", "DIlarge")}.
#'
#' @details Mirrors legacy \code{01_ReprojectAndCropData.r} lines 64-88.
#'   Note: DefInd and largeDefInd ship in Mollweide projection; their CRS
#'   must be set to \code{"+proj=moll"} before reprojecting.
#'
#' @export
load_defaunation_indices <- function(rast_template) {
  di_path <- file.path("data", "DefInd.tif")
  di_large_path <- file.path("data", "largeDefInd.tif")
  bogoni_path <- file.path("data", "Bogoni_DI.tif")
  
  if (!file.exists(di_path) || !file.exists(di_large_path) || !file.exists(bogoni_path)) {
    stop("One or more external defaunation rasters missing from data/.")
  }
  
  DI <- rast(di_path)
  crs(DI) <- "+proj=moll"
  DI_proj <- project(DI, rast_template)
  ext(DI_proj) <- ext(rast_template)
  names(DI_proj) <- "DI"
  
  DIlarge <- rast(di_large_path)
  crs(DIlarge) <- "+proj=moll"
  DIlarge_proj <- project(DIlarge, rast_template)
  ext(DIlarge_proj) <- ext(rast_template)
  names(DIlarge_proj) <- "DIlarge"
  
  Bogoni <- rast(bogoni_path)
  Bogoni_proj <- project(Bogoni, rast_template)
  ext(Bogoni_proj) <- ext(rast_template)
  names(Bogoni_proj) <- "DI_Bogoni"
  
  combined <- c(DI_proj, Bogoni_proj, DIlarge_proj)
  return(combined)
}
