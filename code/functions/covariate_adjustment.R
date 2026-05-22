# =============================================================================
# covariate_adjustment.R
#
# OLS covariate adjustment for the Understory Openness Index (UOI).
#
# Rationale: UOI covaries with environmental gradients (precipitation,
# clay content, elevation, slope, hand, forest fraction) that are unrelated to
# the megafauna signal. We fit an OLS model and subtract the covariate-predicted
# component, yielding a residual UOI that is orthogonal to these drivers.
#
# Model formula:
#   uoi ~ elevation + slope + hnd + precip + clay + forest_fraction
#
# Adjustment applied per-pixel:
#   adjusted_uoi = raw_uoi - Σ(β_j × z_j)
#   where z_j = (x_j - mean_j) / sd_j
#
# This file contains the fitting and application functions.
# =============================================================================

library(terra)
library(dplyr)

# =============================================================================
# MODEL FITTING
# =============================================================================

#' Fit the OLS covariate adjustment model
#'
#' Extracts \code{uoi}, \code{elevation}, \code{slope}, \code{hnd}, \code{precip}, \code{clay}, and
#' \code{forest_fraction} bands from a native-resolution \code{SpatRaster},
#' standardises (z-scores) each predictor, and fits an OLS model:
#' \code{uoi ~ elevation_z + slope_z + hnd_z + precip_z + clay_z + forest_fraction_z}.
#'
#' @param native_stack A \code{terra::SpatRaster} with at least the bands
#'   \code{uoi}, \code{elevation}, \code{slope}, \code{hnd}, \code{precip}, \code{clay},
#'   \code{forest_fraction}.
#'
#' @return A named list:
#'   \describe{
#'     \item{coefficients}{Named numeric vector of OLS slopes (excluding
#'       intercept), one per z-scored predictor.}
#'     \item{means}{Named numeric vector of predictor means used for
#'       standardisation.}
#'     \item{sds}{Named numeric vector of predictor standard deviations.}
#'     \item{intercept}{Numeric scalar. OLS intercept.}
#'     \item{lm_fit}{The \code{lm} object for diagnostics.}
#'   }
#'
#' @export
fit_adjustment_model <- function(native_stack) {
  # Extract SpatRaster values to data frame, ignoring NAs
  df <- as.data.frame(native_stack, na.rm = TRUE)
  
  # Covariate list (identical for GEDI and FRIP)
  covariates <- c("elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  
  for (cov in covariates) {
    if (!cov %in% names(df)) {
      stop(sprintf("Required covariate band '%s' missing from stack.", cov))
    }
  }
  
  if (!"uoi" %in% names(df)) {
    stop("Required response band 'uoi' missing from stack.")
  }
  
  # Compute means and standard deviations for standardisation
  means <- sapply(df[covariates], mean, na.rm = TRUE)
  sds <- sapply(df[covariates], sd, na.rm = TRUE)
  
  # Add standardized versions of covariates
  df_scaled <- df
  for (cov in covariates) {
    df_scaled[[paste0(cov, "_z")]] <- (df[[cov]] - means[cov]) / sds[cov]
  }
  
  # Fit linear model
  formula_str <- paste("uoi ~", paste(paste0(covariates, "_z"), collapse = " + "))
  fit <- lm(as.formula(formula_str), data = df_scaled)
  
  coefs <- coef(fit)
  
  list(
    coefficients = coefs[-1],  # Exclude intercept
    intercept = coefs[1],
    means = means,
    sds = sds,
    lm_fit = fit
  )
}


# =============================================================================
# ADJUSTMENT APPLICATION
# =============================================================================

#' Apply OLS covariate adjustment to UOI
#'
#' Subtracts the covariate-predicted UOI component from the raw UOI band,
#' returning a \code{SpatRaster} with an \code{adjusted_uoi} layer.
#'
#' @param native_stack A \code{terra::SpatRaster} with at least the bands
#'   \code{uoi}, \code{elevation}, \code{slope}, \code{hnd}, \code{precip}, \code{clay},
#'   \code{forest_fraction}.
#' @param model Optional. A model list as returned by
#'   \code{fit_adjustment_model()}. If \code{NULL} (default), the model
#'   is fitted on \code{native_stack} before applying.
#'
#' @return A \code{terra::SpatRaster} with one band (\code{adjusted_uoi}).
#'
#' @details
#' Per-pixel computation:
#' \preformatted{
#'   z_j        = (x_j - mean_j) / sd_j
#'   predicted  = intercept + Σ(β_j × z_j)
#'   adjusted   = raw_uoi - (predicted - intercept)
#'              = raw_uoi - Σ(β_j × z_j)
#' }
#' The intercept cancels so the adjusted values are centred on the raw
#' mean (not zero), preserving interpretability.
#'
#' @export
adjust_uoi_covariates <- function(native_stack, model = NULL) {
  if (is.null(model)) {
    model <- fit_adjustment_model(native_stack)
  }
  
  covariates <- c("elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  
  coefs <- model$coefficients
  means <- model$means
  sds <- model$sds
  
  uoi_raw <- native_stack[["uoi"]]
  
  # Start with an adjustment of 0
  adjustment <- rast(uoi_raw, vals = 0)
  
  for (cov in covariates) {
    coef_name <- paste0(cov, "_z")
    beta_j <- coefs[coef_name]
    if (is.na(beta_j)) {
      beta_j <- 0
    }
    
    # Calculate standardised covariate raster
    x_j <- native_stack[[cov]]
    z_j <- (x_j - means[cov]) / sds[cov]
    
    # Accumulate predicted component
    adjustment <- adjustment + (beta_j * z_j)
  }
  
  adjusted_uoi <- uoi_raw - adjustment
  names(adjusted_uoi) <- "adjusted_uoi"
  return(adjusted_uoi)
}
