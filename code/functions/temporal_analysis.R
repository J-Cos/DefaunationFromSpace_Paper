# =============================================================================
# temporal_analysis.R
#
# H4 (Temporal) statistical functions.
#
# Tests whether the functional signal (FRIP) has strengthened over time
# where defaunation is increasing. Uses the GEE pre-computed `frip_mk_tau`
# band (Mann-Kendall trend τ of annual FRIP across 2001–2023).
#
# Dependencies:
#   terra, dplyr, tibble, broom
# =============================================================================

# --- Imports -----------------------------------------------------------------

library(terra)
library(dplyr)
library(tibble)
library(broom)


# --- Functions ---------------------------------------------------------------

#' Summarise MK-tau by basin
#'
#' Extracts the `frip_mk_tau` band from each basin's multi-scale stack
#' and computes summary statistics (mean, median, sd, quantiles) per basin.
#'
#' @param stacks Named list of SpatRaster stacks ("Congo", "Amazon"),
#'   each containing the `frip_mk_tau` band.
#'
#' @return A \code{\link[tibble]{tibble}} with columns:
#'   \describe{
#'     \item{basin}{Character. Basin name.}
#'     \item{mean_tau}{Numeric. Mean Mann-Kendall tau.}
#'     \item{median_tau}{Numeric. Median Mann-Kendall tau.}
#'     \item{sd_tau}{Numeric. Standard deviation of tau.}
#'     \item{q25}{Numeric. 25th percentile.}
#'     \item{q75}{Numeric. 75th percentile.}
#'     \item{n_pixels}{Integer. Number of non-NA pixels.}
#'   }
#'
#' @examples
#' # tau_summary <- summarise_mk_tau_by_basin(stacks)
#' # tau_summary
summarise_mk_tau_by_basin <- function(stacks) {
  if (length(stacks) == 0) {
    stop("stacks list is empty")
  }
  
  results <- list()
  for (basin_name in names(stacks)) {
    stack <- stacks[[basin_name]]
    if (!"frip_mk_tau" %in% names(stack)) {
      stop(sprintf("frip_mk_tau band missing from stack for %s", basin_name))
    }
    tau_vals <- values(stack[["frip_mk_tau"]], mat = FALSE)
    tau_vals <- tau_vals[!is.na(tau_vals)]
    
    if (length(tau_vals) == 0) {
      results[[basin_name]] <- tibble(
        basin = basin_name,
        mean_tau = NA_real_,
        median_tau = NA_real_,
        sd_tau = NA_real_,
        q25 = NA_real_,
        q75 = NA_real_,
        n_pixels = 0L
      )
    } else {
      results[[basin_name]] <- tibble(
        basin = basin_name,
        mean_tau = mean(tau_vals),
        median_tau = median(tau_vals),
        sd_tau = sd(tau_vals),
        q25 = as.numeric(quantile(tau_vals, probs = 0.25)),
        q75 = as.numeric(quantile(tau_vals, probs = 0.75)),
        n_pixels = length(tau_vals)
      )
    }
  }
  bind_rows(results)
}


#' MK-tau intercept-free model by PA
#'
#' Fits an intercept-free linear model `tau ~ 0 + PA_NAME` to estimate
#' the mean MK-tau trend for each protected area. Classifies each PA as
#' "Increasing", "Decreasing", or "None" based on the sign and significance
#' of its coefficient.
#'
#' @param stacks Named list of SpatRaster stacks containing `frip_mk_tau`.
#' @param pa_rast SpatRaster. Raster with integer PA IDs or PA name factor.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{model}{The fitted \code{lm} object.}
#'     \item{coefficients}{A tibble from \code{broom::tidy()} with an added
#'       `trend_class` column ("Increasing", "Decreasing", or "None").}
#'     \item{summary}{A tibble from \code{broom::glance()}.}
#'   }
#'
#' @examples
#' # res <- test_mk_tau_by_pa(stacks, pa_rast)
#' # res$coefficients |> filter(trend_class == "Increasing")
test_mk_tau_by_pa <- function(stacks, pa_rast) {
  if (length(stacks) == 0) {
    stop("stacks list is empty")
  }
  
  dfs <- list()
  for (basin_name in names(stacks)) {
    stack <- stacks[[basin_name]]
    if (!"frip_mk_tau" %in% names(stack)) {
      stop(sprintf("frip_mk_tau band missing from stack for %s", basin_name))
    }
    tau_vals <- values(stack[["frip_mk_tau"]], mat = FALSE)
    pa_r <- resample(pa_rast, stack, method = "near")
    pa_vals <- values(pa_r, mat = FALSE)
    
    df_b <- data.frame(
      frip_mk_tau = tau_vals,
      pa_id = pa_vals,
      region = basin_name,
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(frip_mk_tau), !is.na(pa_id), pa_id > 0)
    
    dfs[[basin_name]] <- df_b
  }
  
  df <- bind_rows(dfs)
  
  if (nrow(df) == 0) {
    stop("No protected area pixels found with valid MK-tau.")
  }
  
  # Map categories to PA names if available
  pa_cats <- cats(pa_rast)
  if (!is.null(pa_cats) && length(pa_cats) > 0 && !is.null(pa_cats[[1]])) {
    cat_df <- pa_cats[[1]]
    name_col <- grep("NAME", names(cat_df), value = TRUE, ignore.case = TRUE)
    if (length(name_col) > 0) {
      df <- df %>%
        left_join(cat_df, by = c("pa_id" = names(cat_df)[1])) %>%
        rename(PA_NAME = !!sym(name_col[1]))
    } else {
      df$PA_NAME <- as.factor(df$pa_id)
    }
  } else {
    df$PA_NAME <- as.factor(df$pa_id)
  }
  
  df$PA_NAME <- as.factor(df$PA_NAME)
  
  # Fit linear model (intercept-free)
  fit <- lm(frip_mk_tau ~ 0 + PA_NAME, data = df)
  
  coefs_tidy <- broom::tidy(fit) %>%
    mutate(
      trend_class = case_when(
        estimate > 0 & p.value < 0.05 ~ "Increasing",
        estimate < 0 & p.value < 0.05 ~ "Decreasing",
        TRUE ~ "None"
      )
    )
  
  summary_glance <- broom::glance(fit)
  
  list(
    model = fit,
    coefficients = coefs_tidy,
    summary = summary_glance
  )
}


#' Protection effect on MK-tau
#'
#' Tests whether the Mann-Kendall tau is more positive (i.e. stronger
#' increasing FRIP trend) in unprotected areas compared to protected areas.
#' Uses a one-tailed Welch's t-test (H_a: unprotected tau > protected tau).
#'
#' @param stacks Named list of SpatRaster stacks containing `frip_mk_tau`.
#' @param pa_rast SpatRaster. Binary protected-area raster (1 = protected,
#'   0 = unprotected).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{t_stat}{Numeric. Welch's t-statistic.}
#'     \item{p_value}{Numeric. One-tailed p-value.}
#'     \item{mean_protected}{Numeric. Mean tau in protected pixels.}
#'     \item{mean_unprotected}{Numeric. Mean tau in unprotected pixels.}
#'     \item{cohens_d}{Numeric. Cohen's d effect size.}
#'   }
#'
#' @examples
#' # res <- test_mk_tau_protection(stacks, pa_rast)
#' # res$p_value
test_mk_tau_protection <- function(stacks, pa_rast) {
  if (length(stacks) == 0) {
    stop("stacks list is empty")
  }
  
  dfs <- list()
  for (basin_name in names(stacks)) {
    stack <- stacks[[basin_name]]
    if (!"frip_mk_tau" %in% names(stack)) {
      stop(sprintf("frip_mk_tau band missing from stack for %s", basin_name))
    }
    tau_vals <- values(stack[["frip_mk_tau"]], mat = FALSE)
    pa_r <- resample(pa_rast, stack, method = "near")
    pa_vals <- values(pa_r, mat = FALSE)
    
    df_b <- data.frame(
      tau = tau_vals,
      pa_id = pa_vals,
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(tau), !is.na(pa_id))
    
    dfs[[basin_name]] <- df_b
  }
  
  df <- bind_rows(dfs)
  
  protected_tau <- df$tau[df$pa_id > 0]
  unprotected_tau <- df$tau[df$pa_id == 0]
  
  if (length(protected_tau) < 2 || length(unprotected_tau) < 2) {
    stop("Insufficient non-NA data in protected or unprotected pixels.")
  }
  
  # Ha: unprotected > protected
  tt <- t.test(unprotected_tau, protected_tau, alternative = "greater", var.equal = FALSE)
  
  d_res <- effsize::cohen.d(unprotected_tau, protected_tau)
  
  list(
    t_stat = as.numeric(tt$statistic),
    p_value = as.numeric(tt$p.value),
    mean_protected = mean(protected_tau),
    mean_unprotected = mean(unprotected_tau),
    cohens_d = as.numeric(d_res$estimate)
  )
}


#' Multi-scale MK-tau analysis loop
#'
#' Runs the temporal analyses at each of 20 aggregation scales (5 km to
# 100 km). Summarises basin-level tau, per-PA trends, and protection
# effects at each scale.
#'
#' @param multiscale_stacks Named list of lists. Outer names are scales
#'   (character), inner names are basins. Each element is a SpatRaster
#'   with the `frip_mk_tau` band.
#' @param pa_rast SpatRaster. Protected-area raster.
#'
#' @return A \code{\link[tibble]{tibble}} with columns:
#'   \describe{
#'     \item{scale}{Integer. Aggregation scale in metres.}
#'     \item{congo_mean_tau}{Numeric. Congo mean MK-tau at this scale.}
#'     \item{amazon_mean_tau}{Numeric. Amazon mean MK-tau at this scale.}
#'     \item{protection_p}{Numeric. p-value for protection effect.}
#'     \item{protection_d}{Numeric. Cohen's d for protection effect.}
#'     \item{n_increasing_pa}{Integer. Number of PAs with increasing trend.}
#'     \item{n_decreasing_pa}{Integer. Number of PAs with decreasing trend.}
#'   }
#'
#' @examples
#' # h4_ms <- run_h4_multiscale(multiscale_stacks, pa_rast)
#' # h4_ms |> filter(protection_p < 0.05)
run_h4_multiscale <- function(multiscale_stacks, pa_rast) {
  scales <- names(multiscale_stacks)
  results <- list()
  
  for (scale_str in scales) {
    scale_val <- as.integer(scale_str)
    stack_at_scale <- multiscale_stacks[[scale_str]]
    
    # 1. Summarise by basin
    tau_sum <- summarise_mk_tau_by_basin(stack_at_scale)
    congo_mean <- tau_sum$mean_tau[tau_sum$basin == "Congo"]
    amazon_mean <- tau_sum$mean_tau[tau_sum$basin == "Amazon"]
    if (length(congo_mean) == 0) congo_mean <- NA_real_
    if (length(amazon_mean) == 0) amazon_mean <- NA_real_
    
    # 2. Protection effect
    prot_res <- tryCatch(
      test_mk_tau_protection(stack_at_scale, pa_rast),
      error = function(e) list(p_value = NA_real_, cohens_d = NA_real_)
    )
    
    # 3. PA-level trends
    pa_res <- tryCatch(
      test_mk_tau_by_pa(stack_at_scale, pa_rast),
      error = function(e) NULL
    )
    
    n_inc <- 0L
    n_dec <- 0L
    if (!is.null(pa_res)) {
      n_inc <- sum(pa_res$coefficients$trend_class == "Increasing", na.rm = TRUE)
      n_dec <- sum(pa_res$coefficients$trend_class == "Decreasing", na.rm = TRUE)
    }
    
    results[[scale_str]] <- tibble(
      scale = scale_val,
      congo_mean_tau = congo_mean,
      amazon_mean_tau = amazon_mean,
      protection_p = prot_res$p_value,
      protection_d = prot_res$cohens_d,
      n_increasing_pa = n_inc,
      n_decreasing_pa = n_dec
    )
  }
  bind_rows(results)
}
