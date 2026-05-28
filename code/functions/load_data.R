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

#' Automatically resolves the path to the analysis stack directory
#' @details Checks environment variables, standard local sync folders, and GVFS mounts.
#' @export
resolve_data_dir <- function(default_dir = file.path("outputs", "synthetic_EOdata")) {
  # 1. Environment variable override
  env_dir <- Sys.getenv("GEE_DRIVE_DIR")
  if (env_dir != "") {
    if (dir.exists(env_dir)) {
      message(sprintf("Using data directory from GEE_DRIVE_DIR: %s", env_dir))
      return(env_dir)
    } else {
      warning(sprintf("GEE_DRIVE_DIR is set to %s, but directory does not exist.", env_dir))
    }
  }

  # 2. Local sync directory candidates
  home_dir <- Sys.getenv("HOME")
  candidates <- c(
    file.path(home_dir, "GoogleDrive", "DefaunationSynthesis", "AnalysisStack"),
    file.path(home_dir, "Google Drive", "DefaunationSynthesis", "AnalysisStack"),
    file.path(home_dir, "gdrive", "DefaunationSynthesis", "AnalysisStack"),
    file.path(home_dir, "GoogleDrive-MyDrive", "DefaunationSynthesis", "AnalysisStack")
  )
  
  # 3. GVFS google-drive mount detection (Linux Gnome Online Accounts)
  run_dir <- "/run/user"
  if (dir.exists(run_dir)) {
    uids <- list.files(run_dir)
    for (uid in uids) {
      gvfs_dir <- file.path(run_dir, uid, "gvfs")
      if (dir.exists(gvfs_dir)) {
        mounts <- list.files(gvfs_dir, pattern = "^google-drive")
        for (m in mounts) {
          path_with_my_drive <- file.path(gvfs_dir, m, "My Drive", "DefaunationSynthesis", "AnalysisStack")
          path_without_my_drive <- file.path(gvfs_dir, m, "DefaunationSynthesis", "AnalysisStack")
          candidates <- c(candidates, path_with_my_drive, path_without_my_drive)
        }
      }
    }
  }
  
  # 4. Standard Mac and Windows Google Drive paths
  candidates <- c(candidates,
    file.path("G:", "My Drive", "DefaunationSynthesis", "AnalysisStack"),
    file.path("/Volumes", "GoogleDrive", "My Drive", "DefaunationSynthesis", "AnalysisStack")
  )

  # Standardize and filter candidates
  candidates <- unique(path.expand(candidates))
  
  for (cand in candidates) {
    if (dir.exists(cand)) {
      test_files <- list.files(cand, pattern = "^analysis_stack_.*\\.tif$")
      if (length(test_files) > 0) {
        message(sprintf("✓ Automatically detected Google Drive folder with analysis stacks: %s", cand))
        return(cand)
      }
    }
  }

  # 5. Fallback to default
  message(sprintf("Using default data directory: %s", default_dir))
  return(default_dir)
}

#' Root directory for analysis-ready rasters (GEE exports or synthetic data)
#' @export
DATA_DIR <- resolve_data_dir()

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
    basename <- sprintf("analysis_stack_%d_%s.tif", scale, basin)
    
    # Path checking order:
    # 1. Real GEE exports directory (outputs/EOdata)
    # 2. Directly in passed data_dir
    # 3. Synthetic fallback directory (outputs/synthetic_EOdata)
    candidates <- c(
      file.path("outputs", "EOdata", basename),
      file.path(data_dir, basename),
      file.path("outputs", "synthetic_EOdata", basename)
    )
    
    filename <- NULL
    for (cand in candidates) {
      if (file.exists(cand)) {
        filename <- cand
        break
      }
    }
    
    # Dynamic aggregation fallback from 5,000m real stack if target scale real file is missing
    r_5000_basename <- sprintf("analysis_stack_5000_%s.tif", basin)
    r_5000_candidates <- c(
      file.path("outputs", "EOdata", r_5000_basename),
      file.path(data_dir, r_5000_basename)
    )
    r_5000_filename <- NULL
    for (cand in r_5000_candidates) {
      if (file.exists(cand)) {
        r_5000_filename <- cand
        break
      }
    }
    
    is_real_file_missing <- !file.exists(file.path("outputs", "EOdata", basename)) && 
                            !file.exists(file.path(data_dir, basename))
    
    if (is_real_file_missing && !is.null(r_5000_filename) && scale > 5000) {
      fact <- scale / 5000
      message(sprintf("✓ Dynamically aggregating real 5,000m stack by factor of %d to %d m: %s", fact, scale, basename))
      r_5000 <- rast(r_5000_filename)
      r <- terra::aggregate(r_5000, fact = fact, fun = "mean", na.rm = TRUE)
      
      num_layers <- nlyr(r)
      if (num_layers == 12) {
        names(r) <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                      "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
      } else if (num_layers == 11) {
        names(r) <- c("frip", "frip_mk_tau", "uoi", "rh98", "gedi_n",
                      "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
      }
      
      stacks[[as.character(scale)]] <- r
      next
    }
    
    if (is.null(filename)) {
      stop("File does not exist in any candidate location: ", basename)
    }
    
    # Notify if loading real GEE export
    if (grepl("outputs/EOdata", filename, fixed = TRUE)) {
      message(sprintf("✓ Loading real GEE GeoTIFF stack: %s", filename))
    }
    
    r <- rast(filename)
    num_layers <- nlyr(r)
    if (num_layers == 12) {
      names(r) <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                    "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
    } else if (num_layers == 11) {
      names(r) <- c("frip", "frip_mk_tau", "uoi", "rh98", "gedi_n",
                    "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
    } else {
      stop(sprintf("Unexpected number of bands (%d) in multiscale stack: %s", num_layers, filename))
    }
    
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
  basename <- sprintf("analysis_stack_native_%s.tif", basin)
  
  # Path checking order:
  # 1. Real GEE exports directory (outputs/EOdata)
  # 2. Directly in passed data_dir
  # 3. Synthetic fallback directory (outputs/synthetic_EOdata)
  candidates <- c(
    file.path("outputs", "EOdata", basename),
    file.path(data_dir, basename),
    file.path("outputs", "synthetic_EOdata", basename)
  )
  
  filename <- NULL
  for (cand in candidates) {
    if (file.exists(cand)) {
      filename <- cand
      break
    }
  }
  
  if (is.null(filename)) {
    stop("File does not exist in any candidate location: ", basename)
  }
  
  # Notify if loading real GEE export
  if (grepl("outputs/EOdata", filename, fixed = TRUE)) {
    message(sprintf("✓ Loading real native GeoTIFF stack: %s", filename))
  }
  
  r <- rast(filename)
  num_layers <- nlyr(r)
  if (num_layers == 11) {
    names(r) <- c("uoi", "uoi_sd", "rh98", "gedi_n", "elevation", "slope", "hnd",
                  "precip", "clay", "forest_fraction", "Npp_median")
  } else if (num_layers == 10) {
    names(r) <- c("uoi", "rh98", "gedi_n", "elevation", "slope", "hnd",
                  "precip", "clay", "forest_fraction", "Npp_median")
  } else {
    stop(sprintf("Unexpected number of bands (%d) in native stack: %s", num_layers, filename))
  }
  
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
