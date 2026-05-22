# =============================================================================
# gedi_analysis.R
#
# H1 (Structural) statistical functions.
#
# Tests whether forest understories are more open where megafauna are intact.
# Prediction: Congo > Amazon for UOI (Understory Openness Index).
#
# Dependencies:
#   terra, dplyr, tibble, broom, effsize, multcompView
# =============================================================================

# --- Imports -----------------------------------------------------------------

library(terra)
library(dplyr)
library(tibble)
library(broom)
library(effsize)
library(multcompView)


# --- Functions ---------------------------------------------------------------

#' One-tailed Welch's t-test: Congo UOI > Amazon UOI
#'
#' Extracts the `uoi` band from each basin's native-scale stack and performs
#' a one-tailed Welch's t-test (H1: Congo mean > Amazon mean). Also computes
#' Cohen's d effect size.
#'
#' @param native_stacks Named list of SpatRaster stacks (names: "Congo", "Amazon"),
#'   each containing at least the `uoi` band.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{t_stat}{Numeric. The t-statistic.}
#'     \item{p_value}{Numeric. One-tailed p-value (H_a: Congo > Amazon).}
#'     \item{cohens_d}{Numeric. Cohen's d effect size.}
#'     \item{congo_mean}{Numeric. Mean UOI for Congo.}
#'     \item{amazon_mean}{Numeric. Mean UOI for Amazon.}
#'   }
#'
#' @examples
#' # result <- test_regional_uoi(list(Congo = congo_stack, Amazon = amazon_stack))
#' # result$p_value
test_regional_uoi <- function(native_stacks) {
  if (length(native_stacks) == 0 || !all(c("Congo", "Amazon") %in% names(native_stacks))) {
    stop("Not yet implemented: native_stacks must contain Congo and Amazon")
  }
  
  congo_uoi <- values(native_stacks[["Congo"]][["uoi"]], mat = FALSE)
  amazon_uoi <- values(native_stacks[["Amazon"]][["uoi"]], mat = FALSE)
  
  congo_uoi <- congo_uoi[!is.na(congo_uoi)]
  amazon_uoi <- amazon_uoi[!is.na(amazon_uoi)]
  
  if (length(congo_uoi) < 2 || length(amazon_uoi) < 2) {
    stop("Insufficient non-NA data in one or both basins.")
  }
  
  tt <- t.test(congo_uoi, amazon_uoi, alternative = "greater", var.equal = FALSE)
  
  d_res <- effsize::cohen.d(congo_uoi, amazon_uoi)
  
  list(
    t_stat = as.numeric(tt$statistic),
    p_value = as.numeric(tt$p.value),
    cohens_d = as.numeric(d_res$estimate),
    congo_mean = mean(congo_uoi),
    amazon_mean = mean(amazon_uoi)
  )
}


#' Two-way ANOVA: UOI ~ Region × Protection
#'
#' Rasterises the protected-area layer, extracts UOI values, assigns each pixel
#' a Region (Congo/Amazon) and Protection (Protected/Unprotected) label, then
#' fits a two-way ANOVA with interaction.
#'
#' @param native_stacks Named list of SpatRaster stacks ("Congo", "Amazon").
#' @param pa_rast SpatRaster. Binary raster of protected-area status (1 = protected,
#'   0 = unprotected), aligned to the native stacks.
#'
#' @return An `aov` model object for `uoi ~ region * protection`.
#'
#' @examples
#' # model <- test_uoi_region_protection(native_stacks, pa_rast)
#' # summary(model)
test_uoi_region_protection <- function(native_stacks, pa_rast) {
  if (length(native_stacks) == 0 || !all(c("Congo", "Amazon") %in% names(native_stacks))) {
    stop("native_stacks must contain Congo and Amazon")
  }
  
  congo_uoi <- values(native_stacks[["Congo"]][["uoi"]], mat = FALSE)
  pa_congo_r <- resample(pa_rast, native_stacks[["Congo"]], method = "near")
  pa_congo <- values(pa_congo_r, mat = FALSE)
  
  amazon_uoi <- values(native_stacks[["Amazon"]][["uoi"]], mat = FALSE)
  pa_amazon_r <- resample(pa_rast, native_stacks[["Amazon"]], method = "near")
  pa_amazon <- values(pa_amazon_r, mat = FALSE)
  
  df_c <- data.frame(
    uoi = congo_uoi,
    protection = ifelse(pa_congo > 0, "Protected", "Unprotected"),
    region = "Congo",
    stringsAsFactors = FALSE
  )
  
  df_a <- data.frame(
    uoi = amazon_uoi,
    protection = ifelse(pa_amazon > 0, "Protected", "Unprotected"),
    region = "Amazon",
    stringsAsFactors = FALSE
  )
  
  df <- rbind(df_c, df_a) %>%
    filter(!is.na(uoi), !is.na(protection))
  
  df$region <- as.factor(df$region)
  df$protection <- as.factor(df$protection)
  
  model <- aov(uoi ~ region * protection, data = df)
  return(model)
}


#' Per-PA ANOVA + TukeyHSD + compact letter display
#'
#' Groups pixels by individual protected area name, fits a one-way ANOVA
#' (`uoi ~ PA_NAME`), runs TukeyHSD post-hoc, and converts to compact letter
#' display using \code{multcompView::multcompLetters}.
#'
#' @param native_stacks Named list of SpatRaster stacks ("Congo", "Amazon").
#' @param pa_rast SpatRaster. Raster with integer PA IDs (0 = unprotected).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{aov}{The `aov` model object.}
#'     \item{tukey}{The TukeyHSD result object.}
#'     \item{letters}{A named character vector of compact letter display groups.}
#'   }
#'
#' @examples
#' # res <- test_uoi_by_pa(native_stacks, pa_rast)
#' # res$letters
test_uoi_by_pa <- function(native_stacks, pa_rast) {
  if (length(native_stacks) == 0 || !all(c("Congo", "Amazon") %in% names(native_stacks))) {
    stop("native_stacks must contain Congo and Amazon")
  }
  
  congo_uoi <- values(native_stacks[["Congo"]][["uoi"]], mat = FALSE)
  pa_congo_r <- resample(pa_rast, native_stacks[["Congo"]], method = "near")
  pa_congo <- values(pa_congo_r, mat = FALSE)
  
  amazon_uoi <- values(native_stacks[["Amazon"]][["uoi"]], mat = FALSE)
  pa_amazon_r <- resample(pa_rast, native_stacks[["Amazon"]], method = "near")
  pa_amazon <- values(pa_amazon_r, mat = FALSE)
  
  df_c <- data.frame(uoi = congo_uoi, pa_id = pa_congo, region = "Congo", stringsAsFactors = FALSE)
  df_a <- data.frame(uoi = amazon_uoi, pa_id = pa_amazon, region = "Amazon", stringsAsFactors = FALSE)
  
  df <- rbind(df_c, df_a) %>%
    filter(!is.na(uoi), !is.na(pa_id), pa_id > 0)
  
  if (nrow(df) == 0) {
    stop("No protected area pixels found.")
  }
  
  # Map category table if available
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
  
  df$PA_NAME <- gsub("-", " ", df$PA_NAME)
  df$PA_NAME <- as.factor(df$PA_NAME)
  
  fit <- aov(uoi ~ PA_NAME, data = df)
  tuk <- TukeyHSD(fit)
  
  p_vals <- setNames(tuk[[1]][, "p adj"], rownames(tuk[[1]]))
  mc_let <- multcompView::multcompLetters(p_vals)
  
  list(
    aov = fit,
    tukey = tuk,
    letters = mc_let$Letters
  )
}


#' Multi-scale regional UOI t-test loop
#'
#' Runs \code{test_regional_uoi} at each of 20 aggregation scales
#' (5 km to 100 km in 5 km steps). Returns a tidy tibble of results with
#' significance flags.
#'
#' @param multiscale_stacks Named list of lists. Outer names are scales
#'   (e.g. "5000", "10000"), inner names are basins ("Congo", "Amazon").
#'   Each element is a SpatRaster with the `uoi` band.
#' @param pa_rast SpatRaster. Protected-area raster (used for filtering if needed).
#'
#' @return A \code{\link[tibble]{tibble}} with columns:
#'   \describe{
#'     \item{scale}{Integer. Aggregation scale in metres.}
#'     \item{t_stat}{Numeric. t-statistic.}
#'     \item{p_value}{Numeric. One-tailed p-value.}
#'     \item{cohens_d}{Numeric. Effect size.}
#'     \item{signif}{Logical. TRUE if p < 0.05.}
#'   }
#'
#' @examples
#' # ms_results <- run_h1_multiscale(multiscale_stacks, pa_rast)
#' # ms_results |> filter(signif)
run_h1_multiscale <- function(multiscale_stacks, pa_rast) {
  scales <- names(multiscale_stacks)
  results <- list()
  for (scale_str in scales) {
    scale_val <- as.integer(scale_str)
    stack_at_scale <- multiscale_stacks[[scale_str]]
    res <- test_regional_uoi(stack_at_scale)
    results[[scale_str]] <- tibble(
      scale = scale_val,
      t_stat = res$t_stat,
      p_value = res$p_value,
      cohens_d = res$cohens_d,
      signif = res$p_value < 0.05
    )
  }
  bind_rows(results)
}
