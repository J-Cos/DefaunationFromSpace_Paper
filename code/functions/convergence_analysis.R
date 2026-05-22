# =============================================================================
# convergence_analysis.R
#
# H3 (Convergence) statistical functions.
#
# Tests whether the structural (UOI) and functional (FRIP) defaunation
# signals spatially converge: protected areas and pixels with high UOI
# (intact structure) should have low FRIP (intact function).
#
# Dependencies:
#   terra, dplyr, tibble
# =============================================================================

# --- Imports -----------------------------------------------------------------

library(terra)
library(dplyr)
library(tibble)


# --- Functions ---------------------------------------------------------------

#' PA-scale convergence: Spearman r between mean UOI and mean FRIP
#'
#' Computes the mean UOI and mean FRIP within each protected area, then
#' calculates the Spearman rank correlation between them. Expects a negative
#' correlation if the signals converge (high UOI = intact structure → low FRIP
#' = intact nutrient pump).
#'
#' @param native_stacks Named list of SpatRaster stacks ("Congo", "Amazon"),
#'   each containing `uoi` and `frip` bands.
#' @param pa_rast SpatRaster. Raster with integer PA IDs (0 = unprotected).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{rho}{Numeric. Spearman correlation coefficient.}
#'     \item{p_value}{Numeric. p-value for the correlation test.}
#'     \item{pa_df}{A tibble with columns `pa_id`, `mean_uoi`, `mean_frip`,
#'       `basin`.}
#'   }
#'
#' @examples
#' # res <- test_pa_convergence(native_stacks, pa_rast)
#' # res$rho
test_pa_convergence <- function(native_stacks, pa_rast) {
  frip_band <- if ("frip" %in% names(native_stacks[["Congo"]])) "frip" else if ("frip_mean" %in% names(native_stacks[["Congo"]])) "frip_mean" else "Npp_median"
  
  congo_uoi <- values(native_stacks[["Congo"]][["uoi"]], mat = FALSE)
  congo_frip <- values(native_stacks[["Congo"]][[frip_band]], mat = FALSE)
  pa_congo_r <- resample(pa_rast, native_stacks[["Congo"]], method = "near")
  pa_congo <- values(pa_congo_r, mat = FALSE)
  
  amazon_uoi <- values(native_stacks[["Amazon"]][["uoi"]], mat = FALSE)
  amazon_frip <- values(native_stacks[["Amazon"]][[frip_band]], mat = FALSE)
  pa_amazon_r <- resample(pa_rast, native_stacks[["Amazon"]], method = "near")
  pa_amazon <- values(pa_amazon_r, mat = FALSE)
  
  df_c <- data.frame(uoi = congo_uoi, frip = congo_frip, pa_id = pa_congo, basin = "Congo", stringsAsFactors = FALSE)
  df_a <- data.frame(uoi = amazon_uoi, frip = amazon_frip, pa_id = pa_amazon, basin = "Amazon", stringsAsFactors = FALSE)
  
  df <- rbind(df_c, df_a) %>%
    filter(!is.na(uoi), !is.na(frip), !is.na(pa_id), pa_id > 0)
  
  if (nrow(df) == 0) {
    stop("No protected area pixels found.")
  }
  
  pa_df <- df %>%
    group_by(pa_id, basin) %>%
    summarise(mean_uoi = mean(uoi), mean_frip = mean(frip), .groups = "drop")
  
  if (nrow(pa_df) < 3) {
    stop("Insufficient PAs for correlation test.")
  }
  
  ct <- cor.test(pa_df$mean_uoi, pa_df$mean_frip, method = "spearman", exact = FALSE)
  
  list(
    rho = as.numeric(ct$estimate),
    p_value = as.numeric(ct$p.value),
    pa_df = as_tibble(pa_df)
  )
}


#' Pixel-scale convergence: Spearman UOI × FRIP per basin
#'
#' Computes the pixel-wise Spearman correlation between UOI and FRIP
#' separately for each basin. Uses a random subsample if pixel count
#' exceeds a threshold to manage computation time.
#'
#' @param native_stacks Named list of SpatRaster stacks ("Congo", "Amazon"),
#'   each containing `uoi` and `frip` bands.
#'
#' @return A list with one element per basin, each containing:
#'   \describe{
#'     \item{rho}{Numeric. Spearman r.}
#'     \item{p_value}{Numeric. p-value.}
#'     \item{n_pixels}{Integer. Number of pixels used.}
#'   }
#'
#' @examples
#' # res <- test_pixel_convergence(native_stacks)
#' # res$Congo$rho
test_pixel_convergence <- function(native_stacks) {
  res <- list()
  frip_band <- if ("frip" %in% names(native_stacks[[1]])) "frip" else if ("frip_mean" %in% names(native_stacks[[1]])) "frip_mean" else "Npp_median"
  
  for (basin in names(native_stacks)) {
    uoi_vals <- values(native_stacks[[basin]][["uoi"]], mat = FALSE)
    frip_vals <- values(native_stacks[[basin]][[frip_band]], mat = FALSE)
    
    df <- data.frame(uoi = uoi_vals, frip = frip_vals) %>%
      filter(!is.na(uoi), !is.na(frip))
    
    n_pixels <- nrow(df)
    if (n_pixels == 0) {
      res[[basin]] <- list(rho = NA_real_, p_value = NA_real_, n_pixels = 0L)
      next
    }
    
    if (n_pixels > 10000) {
      set.seed(42)
      df <- df[sample(seq_len(n_pixels), 10000), ]
    }
    
    ct <- cor.test(df$uoi, df$frip, method = "spearman", exact = FALSE)
    
    res[[basin]] <- list(
      rho = as.numeric(ct$estimate),
      p_value = as.numeric(ct$p.value),
      n_pixels = n_pixels
    )
  }
  return(res)
}


#' Bivariate classification: 2×2 high/low UOI × high/low FRIP
#'
#' Classifies each pixel into one of four categories based on whether
#' UOI and FRIP are above or below their respective basin medians:
#'
#' | Class | UOI   | FRIP  | Interpretation            |
#' |-------|-------|-------|---------------------------|
#' |   1   | High  | Low   | Intact (convergent)       |
#' |   2   | High  | High  | Structure intact, func.   |
#' |   3   | Low   | Low   | Structure lost, func. ok  |
#' |   4   | Low   | High  | Degraded (convergent)     |
#'
#' @param native_stacks Named list of SpatRaster stacks ("Congo", "Amazon"),
#'   each containing `uoi` and `frip` bands.
#'
#' @return A SpatRaster with integer values 1–4 corresponding to the
#'   bivariate classes defined above. Contains both basins mosaicked.
#'
#' @examples
#' # bivar <- classify_bivariate(native_stacks)
#' # plot(bivar)
classify_bivariate <- function(native_stacks) {
  frip_band <- if ("frip" %in% names(native_stacks[[1]])) "frip" else if ("frip_mean" %in% names(native_stacks[[1]])) "frip_mean" else "Npp_median"
  
  basin_rasters <- list()
  for (basin in names(native_stacks)) {
    stack <- native_stacks[[basin]]
    uoi_r <- stack[["uoi"]]
    frip_r <- stack[[frip_band]]
    
    uoi_vals <- values(uoi_r, mat = FALSE)
    frip_vals <- values(frip_r, mat = FALSE)
    
    uoi_med <- median(uoi_vals, na.rm = TRUE)
    frip_med <- median(frip_vals, na.rm = TRUE)
    
    bivar_vals <- rep(NA_integer_, length(uoi_vals))
    
    is_high_uoi <- !is.na(uoi_vals) & uoi_vals >= uoi_med
    is_low_uoi <- !is.na(uoi_vals) & uoi_vals < uoi_med
    is_high_frip <- !is.na(frip_vals) & frip_vals >= frip_med
    is_low_frip <- !is.na(frip_vals) & frip_vals < frip_med
    
    bivar_vals[is_high_uoi & is_low_frip] <- 1L
    bivar_vals[is_high_uoi & is_high_frip] <- 2L
    bivar_vals[is_low_uoi & is_low_frip] <- 3L
    bivar_vals[is_low_uoi & is_high_frip] <- 4L
    
    bivar_r <- rast(uoi_r)
    values(bivar_r) <- bivar_vals
    names(bivar_r) <- "bivariate_class"
    basin_rasters[[basin]] <- bivar_r
  }
  
  mosaicked <- do.call(terra::mosaic, unname(basin_rasters))
  return(mosaicked)
}


#' Multi-scale UOI–FRIP convergence
#'
#' Computes the Spearman correlation between UOI and FRIP at each of the
#' 20 aggregation scales to assess scale-dependence of convergence.
#'
#' @param multiscale_stacks Named list of lists. Outer names are scales,
#'   inner names are basins. Each element is a SpatRaster with `uoi` and
#'   `frip` bands.
#' @param pa_rast SpatRaster. Protected-area raster (for PA-level convergence
#'   at each scale).
#'
#' @return A \code{\link[tibble]{tibble}} with columns:
#'   \describe{
#'     \item{scale}{Integer. Scale in metres.}
#'     \item{rho_pixel}{Numeric. Pixel-wise Spearman r (pooled).}
#'     \item{rho_pa}{Numeric. PA-level Spearman r.}
#'     \item{p_pixel}{Numeric. p-value for pixel-wise correlation.}
#'     \item{p_pa}{Numeric. p-value for PA-level correlation.}
#'   }
#'
#' @examples
#' # ms_conv <- test_convergence_multiscale(multiscale_stacks, pa_rast)
#' # ms_conv |> filter(p_pixel < 0.05)
test_convergence_multiscale <- function(multiscale_stacks, pa_rast) {
  scales <- names(multiscale_stacks)
  results <- list()
  
  for (scale_str in scales) {
    scale_val <- as.integer(scale_str)
    stack_at_scale <- multiscale_stacks[[scale_str]]
    
    uoi_all <- c()
    frip_all <- c()
    
    frip_band <- if ("frip" %in% names(stack_at_scale[["Congo"]])) "frip" else if ("frip_mean" %in% names(stack_at_scale[["Congo"]])) "frip_mean" else "Npp_median"
    
    for (basin in names(stack_at_scale)) {
      uoi_all <- c(uoi_all, values(stack_at_scale[[basin]][["uoi"]], mat = FALSE))
      frip_all <- c(frip_all, values(stack_at_scale[[basin]][[frip_band]], mat = FALSE))
    }
    
    df_pix <- data.frame(uoi = uoi_all, frip = frip_all) %>%
      filter(!is.na(uoi), !is.na(frip))
    
    rho_pixel <- NA_real_
    p_pixel <- NA_real_
    if (nrow(df_pix) >= 3) {
      if (nrow(df_pix) > 10000) {
        set.seed(42)
        df_pix <- df_pix[sample(seq_len(nrow(df_pix)), 10000), ]
      }
      ct_pix <- cor.test(df_pix$uoi, df_pix$frip, method = "spearman", exact = FALSE)
      rho_pixel <- as.numeric(ct_pix$estimate)
      p_pixel <- as.numeric(ct_pix$p.value)
    }
    
    rho_pa <- NA_real_
    p_pa <- NA_real_
    try({
      res_pa <- test_pa_convergence(stack_at_scale, pa_rast)
      rho_pa <- res_pa$rho
      p_pa <- res_pa$p_value
    }, silent = TRUE)
    
    results[[scale_str]] <- tibble(
      scale = scale_val,
      rho_pixel = rho_pixel,
      rho_pa = rho_pa,
      p_pixel = p_pixel,
      p_pa = p_pa
    )
  }
  bind_rows(results)
}


#' Dual-signal protection test
#'
#' Tests whether both UOI and FRIP differ significantly between protected
#' and unprotected pixels. Performs Welch's t-test for each signal within
#' each basin, expecting:
#'   - UOI: protected > unprotected (structure preserved)
#'   - FRIP: protected < unprotected (function preserved)
#'
#' @param native_stacks Named list of SpatRaster stacks ("Congo", "Amazon"),
#'   each containing `uoi` and `frip` bands.
#' @param pa_rast SpatRaster. Binary protected-area raster (1 = protected,
#'   0 = unprotected).
#'
#' @return A \code{\link[tibble]{tibble}} with columns:
#'   \describe{
#'     \item{basin}{Character. Basin name.}
#'     \item{signal}{Character. "uoi" or "frip".}
#'     \item{mean_protected}{Numeric. Mean value in protected pixels.}
#'     \item{mean_unprotected}{Numeric. Mean value in unprotected pixels.}
#'     \item{t_stat}{Numeric. Welch's t-statistic.}
#'     \item{p_value}{Numeric. p-value (two-tailed).}
#'     \item{direction}{Character. "higher_protected" or "lower_protected".}
#'   }
#'
#' @examples
#' # dual <- test_protection_dual_signal(native_stacks, pa_rast)
#' # dual |> filter(signal == "uoi")
test_protection_dual_signal <- function(native_stacks, pa_rast) {
  frip_band <- if ("frip" %in% names(native_stacks[[1]])) "frip" else if ("frip_mean" %in% names(native_stacks[[1]])) "frip_mean" else "Npp_median"
  
  results <- list()
  for (basin in names(native_stacks)) {
    stack <- native_stacks[[basin]]
    uoi_vals <- values(stack[["uoi"]], mat = FALSE)
    frip_vals <- values(stack[[frip_band]], mat = FALSE)
    
    pa_aligned <- resample(pa_rast, stack, method = "near")
    pa_vals <- values(pa_aligned, mat = FALSE)
    
    df <- data.frame(
      uoi = uoi_vals,
      frip = frip_vals,
      protection = ifelse(pa_vals > 0, "Protected", ifelse(pa_vals == 0, "Unprotected", NA_character_)),
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(protection))
    
    # Welch's t-test for UOI
    df_uoi <- df %>% filter(!is.na(uoi))
    t_uoi <- NA_real_
    p_uoi <- NA_real_
    mean_uoi_prot <- NA_real_
    mean_uoi_unprot <- NA_real_
    dir_uoi <- "none"
    
    if (nrow(df_uoi %>% filter(protection == "Protected")) >= 2 &&
        nrow(df_uoi %>% filter(protection == "Unprotected")) >= 2) {
      tt <- t.test(uoi ~ protection, data = df_uoi, var.equal = FALSE)
      t_uoi <- as.numeric(tt$statistic)
      p_uoi <- as.numeric(tt$p.value)
      mean_uoi_unprot <- as.numeric(tt$estimate[1])
      mean_uoi_prot <- as.numeric(tt$estimate[2])
      dir_uoi <- ifelse(mean_uoi_prot > mean_uoi_unprot, "higher_protected", "lower_protected")
    }
    
    # Welch's t-test for FRIP
    df_frip <- df %>% filter(!is.na(frip))
    t_frip <- NA_real_
    p_frip <- NA_real_
    mean_frip_prot <- NA_real_
    mean_frip_unprot <- NA_real_
    dir_frip <- "none"
    
    if (nrow(df_frip %>% filter(protection == "Protected")) >= 2 &&
        nrow(df_frip %>% filter(protection == "Unprotected")) >= 2) {
      tt <- t.test(frip ~ protection, data = df_frip, var.equal = FALSE)
      t_frip <- as.numeric(tt$statistic)
      p_frip <- as.numeric(tt$p.value)
      mean_frip_unprot <- as.numeric(tt$estimate[1])
      mean_frip_prot <- as.numeric(tt$estimate[2])
      dir_frip <- ifelse(mean_frip_prot > mean_frip_unprot, "higher_protected", "lower_protected")
    }
    
    results[[paste0(basin, "_uoi")]] <- tibble(
      basin = basin,
      signal = "uoi",
      mean_protected = mean_uoi_prot,
      mean_unprotected = mean_uoi_unprot,
      t_stat = t_uoi,
      p_value = p_uoi,
      direction = dir_uoi
    )
    
    results[[paste0(basin, "_frip")]] <- tibble(
      basin = basin,
      signal = "frip",
      mean_protected = mean_frip_prot,
      mean_unprotected = mean_frip_unprot,
      t_stat = t_frip,
      p_value = p_frip,
      direction = dir_frip
    )
  }
  bind_rows(results)
}
