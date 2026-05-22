# =============================================================================
# frip_analysis.R
#
# H2 (Functional) statistical functions.
#
# Tests whether flooding more strongly predicts productivity where megafauna
# are depleted (nutrient pump broken).
# Prediction: Amazon > Congo for FRIP (Flooding Role in Productivity).
#
# Dependencies:
#   terra, dplyr, tibble, broom, MuMIn, multcompView
# =============================================================================

# --- Imports -----------------------------------------------------------------

library(terra)
library(dplyr)
library(tibble)
library(broom)
library(MuMIn)
library(multcompView)


# --- Functions ---------------------------------------------------------------

#' ANOVA: FRIP ~ basin + country + protection
#'
#' Extracts FRIP values from the multi-scale stack and fits a three-way ANOVA
#' with basin, country, and protection status as factors. Performs TukeyHSD
#' post-hoc comparisons and derives compact letter display.
#'
#' @param stack SpatRaster. Multi-scale stack containing at least the `frip` band.
#' @param basins_r SpatRaster. Rasterised basin IDs (1 = Congo, 2 = Amazon).
#' @param countries_r SpatRaster. Rasterised country codes.
#' @param pa_r SpatRaster. Rasterised protected-area status (1 = protected,
#'   0 = unprotected).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{aov}{The `aov` model object.}
#'     \item{tukey}{The TukeyHSD result object.}
#'     \item{letters}{Named character vector of compact letter groups.}
#'   }
#'
#' @examples
#' # res <- test_frip_by_basin_country_pa(stack, basins_r, countries_r, pa_r)
#' # summary(res$aov)
test_frip_by_basin_country_pa <- function(stack, basins_r, countries_r, pa_r) {
  frip_vals <- values(stack[["frip"]], mat = FALSE)
  
  basins_aligned <- resample(basins_r, stack, method = "near")
  basin_vals <- values(basins_aligned, mat = FALSE)
  
  countries_aligned <- resample(countries_r, stack, method = "near")
  country_vals <- values(countries_aligned, mat = FALSE)
  
  pa_aligned <- resample(pa_r, stack, method = "near")
  pa_vals <- values(pa_aligned, mat = FALSE)
  
  df <- data.frame(
    frip = frip_vals,
    basin = ifelse(basin_vals == 1, "Congo", ifelse(basin_vals == 2, "Amazon", NA_character_)),
    country = as.character(country_vals),
    protection = ifelse(pa_vals > 0, "Protected", ifelse(pa_vals == 0, "Unprotected", NA_character_)),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(frip), !is.na(basin), !is.na(country), !is.na(protection))
  
  if (nrow(df) < 5) {
    stop("Insufficient data for 3-way ANOVA.")
  }
  
  df$basin <- as.factor(df$basin)
  df$country <- as.factor(df$country)
  df$protection <- as.factor(df$protection)
  
  fit <- aov(frip ~ basin + country + protection, data = df)
  tuk <- TukeyHSD(fit)
  
  p_vals <- setNames(tuk$protection[, "p adj"], rownames(tuk$protection))
  mc_let <- multcompView::multcompLetters(p_vals)
  
  list(
    aov = fit,
    tukey = tuk,
    letters = mc_let$Letters
  )
}


#' OLS regression: FRIP ~ Defaunation Index
#'
#' Fits an OLS model predicting FRIP from a continuous defaunation index.
#' Returns tidy model summary via \code{broom::tidy} and \code{broom::glance}.
#'
#' @param stack SpatRaster. Multi-scale stack containing the `frip` band.
#' @param di_rast SpatRaster. Continuous defaunation index raster, aligned to
#'   the stack.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{tidy}{A tibble from \code{broom::tidy()} with coefficient estimates.}
#'     \item{glance}{A tibble from \code{broom::glance()} with model-level stats
#'       (R², adjusted R², AIC, etc.).}
#'     \item{model}{The fitted \code{lm} object.}
#'   }
#'
#' @examples
#' # res <- test_frip_vs_di(stack, di_rast)
#' # res$glance$r.squared
test_frip_vs_di <- function(stack, di_rast) {
  frip_vals <- values(stack[["frip"]], mat = FALSE)
  di_aligned <- resample(di_rast, stack, method = "bilinear")
  di_vals <- values(di_aligned, mat = FALSE)
  
  df <- data.frame(frip = frip_vals, di = di_vals) %>%
    filter(!is.na(frip), !is.na(di))
  
  if (nrow(df) < 3) {
    stop("Insufficient data for FRIP vs DI regression.")
  }
  
  fit <- lm(frip ~ di, data = df)
  
  list(
    tidy = broom::tidy(fit),
    glance = broom::glance(fit),
    model = fit
  )
}


#' AIC model selection: compare defaunation indices
#'
#' Fits competing OLS models (FRIP ~ DI_Benitez-Lopez vs FRIP ~ DI_Bogoni)
#' and uses \code{MuMIn::model.sel} for AIC-based model comparison.
#'
#' @param stack SpatRaster. Multi-scale stack containing the `frip` band.
#' @param di_bl SpatRaster. Benitez-Lopez defaunation index raster.
#' @param di_bogoni SpatRaster. Bogoni defaunation index raster.
#'
#' @return A data.frame of model selection results with AIC, delta AIC,
#'   and Akaike weights from \code{MuMIn::model.sel()}.
#'
#' @examples
#' # weights <- compare_di_indices(stack, di_bl, di_bogoni)
#' # weights
compare_di_indices <- function(stack, di_bl, di_bogoni) {
  frip_vals <- values(stack[["frip"]], mat = FALSE)
  bl_aligned <- resample(di_bl, stack, method = "bilinear")
  bl_vals <- values(bl_aligned, mat = FALSE)
  
  bogoni_aligned <- resample(di_bogoni, stack, method = "bilinear")
  bogoni_vals <- values(bogoni_aligned, mat = FALSE)
  
  df <- data.frame(frip = frip_vals, di_bl = bl_vals, di_bogoni = bogoni_vals) %>%
    filter(!is.na(frip), !is.na(di_bl), !is.na(di_bogoni))
  
  if (nrow(df) < 5) {
    stop("Insufficient data for AIC model comparison.")
  }
  
  fit_bl <- lm(frip ~ di_bl, data = df)
  fit_bogoni <- lm(frip ~ di_bogoni, data = df)
  
  sel <- MuMIn::model.sel(fit_bl, fit_bogoni)
  return(as.data.frame(sel))
}


#' Spatial cross-validation: 3×3 tile leave-one-out
#'
#' Divides each basin into a 3×3 spatial tile grid. For each held-out tile,
#' fits a covariate-adjustment model on the remaining 8 tiles, predicts
#' residuals for the held-out tile, and computes R² before and after denoising.
#'
#' @param stack SpatRaster. Multi-scale stack with `frip` and covariate bands.
#' @param di_rast SpatRaster. Defaunation index raster.
#' @param basins_v SpatVector. Basin boundary polygons.
#'
#' @return A \code{\link[tibble]{tibble}} with columns:
#'   \describe{
#'     \item{tile_id}{Integer. Tile index (1–9).}
#'     \item{r2_raw}{Numeric. R² of FRIP ~ DI before denoising.}
#'     \item{r2_denoised}{Numeric. R² of FRIP_residual ~ DI after covariate adjustment.}
#'     \item{ratio}{Numeric. r2_denoised / r2_raw.}
#'   }
#'
#' @examples
#' # cv_results <- run_denoising_cv(stack, di_rast, basins_v)
#' # mean(cv_results$ratio)
run_denoising_cv <- function(stack, di_rast, basins_v) {
  e <- ext(stack)
  x_breaks <- seq(e$xmin, e$xmax, length.out = 4)
  y_breaks <- seq(e$ymin, e$ymax, length.out = 4)
  
  pts <- crds(stack, df = TRUE, na.rm = FALSE)
  col_idx <- cut(pts$x, breaks = x_breaks, labels = FALSE, include.lowest = TRUE)
  row_idx <- cut(pts$y, breaks = y_breaks, labels = FALSE, include.lowest = TRUE)
  tile_id <- (row_idx - 1) * 3 + col_idx
  
  df_full <- as.data.frame(stack, na.rm = FALSE)
  di_aligned <- resample(di_rast, stack, method = "bilinear")
  df_full$di <- values(di_aligned, mat = FALSE)
  df_full$tile_id <- tile_id
  
  covariates <- c("elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  
  df <- df_full %>%
    filter(!is.na(frip), !is.na(di), !is.na(tile_id),
           across(all_of(covariates), ~!is.na(.x)))
  
  results <- tibble(
    tile_id = 1:9,
    r2_raw = NA_real_,
    r2_denoised = NA_real_,
    ratio = NA_real_
  )
  
  if (nrow(df) < 20) {
    return(results)
  }
  
  for (t in 1:9) {
    train_df <- df %>% filter(tile_id != t)
    test_df <- df %>% filter(tile_id == t)
    
    if (nrow(train_df) < 10 || nrow(test_df) < 10) next
    
    means <- sapply(train_df[covariates], mean, na.rm = TRUE)
    sds <- sapply(train_df[covariates], sd, na.rm = TRUE)
    sds[sds == 0] <- 1e-6
    
    train_scaled <- train_df
    test_scaled <- test_df
    for (cov in covariates) {
      train_scaled[[paste0(cov, "_z")]] <- (train_df[[cov]] - means[cov]) / sds[cov]
      test_scaled[[paste0(cov, "_z")]] <- (test_df[[cov]] - means[cov]) / sds[cov]
    }
    
    formula_cov <- as.formula(paste("frip ~", paste(paste0(covariates, "_z"), collapse = " + ")))
    cov_fit <- lm(formula_cov, data = train_scaled)
    
    pred_frip <- predict(cov_fit, newdata = test_scaled)
    test_df$frip_residual <- test_df$frip - pred_frip
    
    lm_raw <- lm(frip ~ di, data = test_df)
    r2_raw <- summary(lm_raw)$r.squared
    
    lm_denoised <- lm(frip_residual ~ di, data = test_df)
    r2_denoised <- summary(lm_denoised)$r.squared
    
    results$r2_raw[results$tile_id == t] <- r2_raw
    results$r2_denoised[results$tile_id == t] <- r2_denoised
    results$ratio[results$tile_id == t] <- r2_denoised / max(r2_raw, 1e-6)
  }
  
  return(results)
}


#' Basin-split denoising validation
#'
#' Train covariate-adjustment model on basin A, validate on basin B, and
#' vice versa. Assesses whether covariate effects are transferable between
#' Congo and Amazon.
#'
#' @param stack SpatRaster. Multi-scale stack with `frip` and covariate bands.
#' @param di_rast SpatRaster. Defaunation index raster.
#' @param basins_v SpatVector. Basin boundary polygons with a `basin` column.
#'
#' @return A \code{\link[tibble]{tibble}} with columns:
#'   \describe{
#'     \item{train_basin}{Character. Training basin name.}
#'     \item{test_basin}{Character. Validation basin name.}
#'     \item{r2_raw}{Numeric. R² before denoising on test basin.}
#'     \item{r2_denoised}{Numeric. R² after denoising on test basin.}
#'     \item{improvement}{Numeric. Proportional R² improvement.}
#'   }
#'
#' @examples
#' # split_res <- run_basin_split_denoising(stack, di_rast, basins_v)
#' # split_res
run_basin_split_denoising <- function(stack, di_rast, basins_v) {
  basin_rast <- terra::rasterize(basins_v, stack, field = "basin")
  basin_vals <- values(basin_rast, mat = FALSE)
  
  df_full <- as.data.frame(stack, na.rm = FALSE)
  df_full$basin <- basin_vals
  
  di_aligned <- resample(di_rast, stack, method = "bilinear")
  df_full$di <- values(di_aligned, mat = FALSE)
  
  covariates <- c("elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  df <- df_full %>%
    filter(!is.na(frip), !is.na(di), !is.na(basin),
           across(all_of(covariates), ~!is.na(.x)))
  
  df_congo <- df %>% filter(basin == "Congo")
  df_amazon <- df %>% filter(basin == "Amazon")
  
  results <- tibble(
    train_basin = c("Congo", "Amazon"),
    test_basin = c("Amazon", "Congo"),
    r2_raw = NA_real_,
    r2_denoised = NA_real_,
    improvement = NA_real_
  )
  
  for (i in 1:2) {
    train_b <- results$train_basin[i]
    test_b <- results$test_basin[i]
    
    train_set <- if (train_b == "Congo") df_congo else df_amazon
    test_set <- if (test_b == "Congo") df_congo else df_amazon
    
    if (nrow(train_set) < 10 || nrow(test_set) < 10) next
    
    means <- sapply(train_set[covariates], mean, na.rm = TRUE)
    sds <- sapply(train_set[covariates], sd, na.rm = TRUE)
    sds[sds == 0] <- 1e-6
    
    train_scaled <- train_set
    test_scaled <- test_set
    for (cov in covariates) {
      train_scaled[[paste0(cov, "_z")]] <- (train_set[[cov]] - means[cov]) / sds[cov]
      test_scaled[[paste0(cov, "_z")]] <- (test_set[[cov]] - means[cov]) / sds[cov]
    }
    
    formula_cov <- as.formula(paste("frip ~", paste(paste0(covariates, "_z"), collapse = " + ")))
    cov_fit <- lm(formula_cov, data = train_scaled)
    
    pred_frip <- predict(cov_fit, newdata = test_scaled)
    test_set$frip_residual <- test_set$frip - pred_frip
    
    lm_raw <- lm(frip ~ di, data = test_set)
    r2_raw <- summary(lm_raw)$r.squared
    
    lm_denoised <- lm(frip_residual ~ di, data = test_set)
    r2_denoised <- summary(lm_denoised)$r.squared
    
    results$r2_raw[i] <- r2_raw
    results$r2_denoised[i] <- r2_denoised
    results$improvement[i] <- (r2_denoised - r2_raw) / max(r2_raw, 1e-6)
  }
  
  return(results)
}


#' Full multi-scale H2 analysis loop
#'
#' Iterates over all 20 aggregation scales, running the core H2 analyses
#' at each scale: FRIP ANOVA, FRIP ~ DI regression, DI index comparison,
#' and spatial cross-validation.
#'
#' @param multiscale_stacks Named list of SpatRaster stacks keyed by scale.
#' @param ... Additional arguments passed to individual test functions
#'   (e.g., `basins_r`, `countries_r`, `pa_r`, `di_rast`, `basins_v`).
#'
#' @return A \code{\link[tibble]{tibble}} with one row per scale and columns
#'   summarising each analysis (R², p-values, AIC weights, etc.).
#'
#' @examples
#' # h2_ms <- run_h2_multiscale(multiscale_stacks, basins_r = br, di_rast = di)
run_h2_multiscale <- function(multiscale_stacks, ...) {
  args <- list(...)
  scales <- names(multiscale_stacks)
  
  results <- list()
  for (scale_str in scales) {
    scale_val <- as.integer(scale_str)
    # Merge Congo and Amazon stacks at this scale
    stack_at_scale <- merge(multiscale_stacks[[scale_str]][["Congo"]],
                            multiscale_stacks[[scale_str]][["Amazon"]])
    
    r2_raw <- NA_real_
    p_val_di <- NA_real_
    if ("di_rast" %in% names(args)) {
      res_di <- test_frip_vs_di(stack_at_scale, args$di_rast)
      r2_raw <- res_di$glance$r.squared
      p_val_di <- res_di$tidy$p.value[res_di$tidy$term == "di"]
      if (length(p_val_di) == 0) p_val_di <- NA_real_
    }
    
    results[[scale_str]] <- tibble(
      scale = scale_val,
      r2_raw = r2_raw,
      p_value_di = p_val_di,
      signif = !is.na(p_val_di) && p_val_di < 0.05
    )
  }
  bind_rows(results)
}
