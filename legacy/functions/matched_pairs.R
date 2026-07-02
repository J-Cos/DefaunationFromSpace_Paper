# =============================================================================
# matched_pairs.R
#
# Matched-pairs analysis for comparing Congo vs Amazon pixel-level UOI.
#
# Ported from: legacy/GEDI_openness/code/matching_algorithm.py
#
# Algorithm:
#   For each Congo pixel (target), find the Amazon pixel (pool) with the
#   smallest normalised distance in (elevation, slope) space. Matching is
#   without replacement — each Amazon pixel can be used at most once.
#
# Distance metric:
#   d = |Δelevation| / elev_tol + |Δslope| / slope_tol
#
# Default tolerances: elevation ±50 m, slope ±2°.
#
# Candidates whose absolute difference exceeds the tolerance in either
# dimension are excluded before distance ranking.
# =============================================================================

library(dplyr)
library(tibble)

# =============================================================================
# PIXEL MATCHING
# =============================================================================

#' Match Congo pixels to Amazon pixels by elevation and slope
#'
#' For each row in \code{df_congo}, finds the closest available row in
#' \code{df_amazon} using a normalised distance metric. Matching is
#' performed without replacement (each Amazon pixel is used at most once).
#' Pixels exceeding the tolerance in either dimension are excluded.
#'
#' @param df_congo A data.frame with at least columns \code{elevation},
#'   \code{slope}, and \code{uoi}. One row per Congo pixel.
#' @param df_amazon A data.frame with at least columns \code{elevation},
#'   \code{slope}, and \code{uoi}. One row per Amazon pixel.
#' @param elev_tol Numeric. Maximum allowed elevation difference in metres
#'   (default 50).
#' @param slope_tol Numeric. Maximum allowed slope difference in degrees
#'   (default 2).
#'
#' @return A \code{tibble} with columns:
#'   \describe{
#'     \item{pair_id}{Integer. Sequential pair identifier.}
#'     \item{congo_elevation, congo_slope, congo_uoi}{Congo pixel values.}
#'     \item{amazon_elevation, amazon_slope, amazon_uoi}{Matched Amazon
#'       pixel values.}
#'     \item{elevation_diff}{Absolute elevation difference (m).}
#'     \item{slope_diff}{Absolute slope difference (°).}
#'     \item{uoi_diff}{Congo UOI minus Amazon UOI.}
#'     \item{distance}{Normalised matching distance.}
#'   }
#'
#' @export
match_pixels <- function(df_congo, df_amazon, elev_tol = 50, slope_tol = 2) {
  # Filter out rows with missing variables we need for matching
  df_c <- df_congo %>%
    filter(!is.na(elevation), !is.na(slope), !is.na(uoi))
  df_a <- df_amazon %>%
    filter(!is.na(elevation), !is.na(slope), !is.na(uoi))
  
  n_congo <- nrow(df_c)
  if (n_congo == 0 || nrow(df_a) == 0) {
    return(tibble(
      pair_id = integer(0),
      congo_elevation = numeric(0),
      congo_slope = numeric(0),
      congo_uoi = numeric(0),
      amazon_elevation = numeric(0),
      amazon_slope = numeric(0),
      amazon_uoi = numeric(0),
      elevation_diff = numeric(0),
      slope_diff = numeric(0),
      uoi_diff = numeric(0),
      distance = numeric(0)
    ))
  }
  
  # Initialize vectors to hold results for speed
  res_pair_id <- integer(n_congo)
  res_c_elev <- numeric(n_congo)
  res_c_slope <- numeric(n_congo)
  res_c_uoi <- numeric(n_congo)
  res_a_elev <- numeric(n_congo)
  res_a_slope <- numeric(n_congo)
  res_a_uoi <- numeric(n_congo)
  res_elev_diff <- numeric(n_congo)
  res_slope_diff <- numeric(n_congo)
  res_uoi_diff <- numeric(n_congo)
  res_distance <- numeric(n_congo)
  
  # Keep track of indices in df_a that are still available
  available <- rep(TRUE, nrow(df_a))
  
  # Pre-extract vectors from df_a for fast operations
  a_elevs <- df_a$elevation
  a_slopes <- df_a$slope
  a_uois <- df_a$uoi
  
  pair_count <- 0
  
  for (i in seq_len(n_congo)) {
    c_elev <- df_c$elevation[i]
    c_slope <- df_c$slope[i]
    c_uoi <- df_c$uoi[i]
    
    # Fast filtering of available candidates within tolerances
    cand_idx <- which(
      available & 
      abs(a_elevs - c_elev) <= elev_tol & 
      abs(a_slopes - c_slope) <= slope_tol
    )
    
    if (length(cand_idx) > 0) {
      # Calculate distances
      elev_diffs <- abs(a_elevs[cand_idx] - c_elev)
      slope_diffs <- abs(a_slopes[cand_idx] - c_slope)
      dists <- (elev_diffs / elev_tol) + (slope_diffs / slope_tol)
      
      # Find best match
      best_rel_idx <- which.min(dists)
      best_abs_idx <- cand_idx[best_rel_idx]
      
      pair_count <- pair_count + 1
      
      # Record match
      res_pair_id[pair_count] <- pair_count
      res_c_elev[pair_count] <- c_elev
      res_c_slope[pair_count] <- c_slope
      res_c_uoi[pair_count] <- c_uoi
      res_a_elev[pair_count] <- a_elevs[best_abs_idx]
      res_a_slope[pair_count] <- a_slopes[best_abs_idx]
      res_a_uoi[pair_count] <- a_uois[best_abs_idx]
      res_elev_diff[pair_count] <- elev_diffs[best_rel_idx]
      res_slope_diff[pair_count] <- slope_diffs[best_rel_idx]
      res_uoi_diff[pair_count] <- c_uoi - a_uois[best_abs_idx]
      res_distance[pair_count] <- dists[best_rel_idx]
      
      # Remove matched Amazon pixel from pool
      available[best_abs_idx] <- FALSE
    }
  }
  
  # Trim arrays to actual pair count
  if (pair_count == 0) {
    return(tibble(
      pair_id = integer(0),
      congo_elevation = numeric(0),
      congo_slope = numeric(0),
      congo_uoi = numeric(0),
      amazon_elevation = numeric(0),
      amazon_slope = numeric(0),
      amazon_uoi = numeric(0),
      elevation_diff = numeric(0),
      slope_diff = numeric(0),
      uoi_diff = numeric(0),
      distance = numeric(0)
    ))
  }
  
  tibble(
    pair_id = res_pair_id[1:pair_count],
    congo_elevation = res_c_elev[1:pair_count],
    congo_slope = res_c_slope[1:pair_count],
    congo_uoi = res_c_uoi[1:pair_count],
    amazon_elevation = res_a_elev[1:pair_count],
    amazon_slope = res_a_slope[1:pair_count],
    amazon_uoi = res_a_uoi[1:pair_count],
    elevation_diff = res_elev_diff[1:pair_count],
    slope_diff = res_slope_diff[1:pair_count],
    uoi_diff = res_uoi_diff[1:pair_count],
    distance = res_distance[1:pair_count]
  )
}

# =============================================================================
# MATCH QUALITY SUMMARY
# =============================================================================

#' Summarise matching quality statistics
#'
#' Reports the number of matched pairs and descriptive statistics for
#' elevation, slope, and UOI differences.
#'
#' @param matched_df A \code{tibble} as returned by \code{match_pixels()}.
#'
#' @return A named list:
#'   \describe{
#'     \item{n_pairs}{Integer. Number of matched pairs.}
#'     \item{elev_diff_mean, elev_diff_sd}{Mean and SD of absolute
#'       elevation differences.}
#'     \item{slope_diff_mean, slope_diff_sd}{Mean and SD of absolute
#'       slope differences.}
#'     \item{uoi_diff_mean, uoi_diff_sd}{Mean and SD of UOI differences
#'       (Congo − Amazon).}
#'   }
#'
#' @export
summarise_matches <- function(matched_df) {
  if (nrow(matched_df) == 0) {
    return(list(
      n_pairs = 0L,
      elev_diff_mean = NA_real_,
      elev_diff_sd = NA_real_,
      slope_diff_mean = NA_real_,
      slope_diff_sd = NA_real_,
      uoi_diff_mean = NA_real_,
      uoi_diff_sd = NA_real_
    ))
  }
  
  list(
    n_pairs = nrow(matched_df),
    elev_diff_mean = mean(matched_df$elevation_diff, na.rm = TRUE),
    elev_diff_sd = sd(matched_df$elevation_diff, na.rm = TRUE),
    slope_diff_mean = mean(matched_df$slope_diff, na.rm = TRUE),
    slope_diff_sd = sd(matched_df$slope_diff, na.rm = TRUE),
    uoi_diff_mean = mean(matched_df$uoi_diff, na.rm = TRUE),
    uoi_diff_sd = sd(matched_df$uoi_diff, na.rm = TRUE)
  )
}

# =============================================================================
# STATISTICAL TEST
# =============================================================================

#' Paired t-test on matched Congo–Amazon UOI differences
#'
#' Tests whether the mean UOI difference (Congo − Amazon) across matched
#' pairs is significantly different from zero.
#'
#' @param matched_df A \code{tibble} as returned by \code{match_pixels()}.
#'
#' @return A named list:
#'   \describe{
#'     \item{t}{Test statistic (t-value).}
#'     \item{p}{Two-sided p-value.}
#'     \item{cohens_d}{Cohen's d effect size: mean_diff / sd_diff.}
#'     \item{mean_diff}{Mean of UOI differences.}
#'     \item{ci}{Numeric vector of length 2: 95\% confidence interval
#'       for the mean difference.}
#'     \item{n_pairs}{Number of pairs used.}
#'   }
#'
#' @export
test_matched_uoi <- function(matched_df) {
  n_pairs <- nrow(matched_df)
  if (n_pairs < 2) {
    return(list(
      t = NA_real_,
      p = NA_real_,
      cohens_d = NA_real_,
      mean_diff = NA_real_,
      ci = c(NA_real_, NA_real_),
      n_pairs = n_pairs
    ))
  }
  
  # Paired t-test is equivalent to one-sample t-test on differences
  tt <- t.test(matched_df$uoi_diff, mu = 0)
  
  # Cohen's d for paired differences is mean(diff) / sd(diff)
  d_val <- mean(matched_df$uoi_diff, na.rm = TRUE) / sd(matched_df$uoi_diff, na.rm = TRUE)
  
  list(
    t = as.numeric(tt$statistic),
    p = as.numeric(tt$p.value),
    cohens_d = d_val,
    mean_diff = as.numeric(tt$estimate),
    ci = as.numeric(tt$conf.int),
    n_pairs = n_pairs
  )
}
