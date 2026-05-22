# =============================================================================
# pa_pairs.R
#
# Named Protected Area (PA) pair definitions and ranked-park lists for
# the megafauna effect analysis.
#
# Data sources:
#   - PA pairs: legacy/GEDI_openness/code/analyze_pa_pairs.py (PAIRS dict)
#   - Ranked parks: legacy/GEDI_openness/code/plot_regional_boxplots.py
#     (AMAZON_PARKS, CONGO_PARKS)
#
# Each pair contrasts a "full" PA (with intact megafauna) against an
# "empty" PA (defaunated), using WDPA polygon IDs for lookup.
# =============================================================================

library(terra)
library(dplyr)
library(tibble)

# =============================================================================
# PA PAIR DEFINITIONS
# =============================================================================

#' Define the 10 named PA pairs for megafauna effect analysis
#'
#' Returns a tibble with one row per pair. Each pair contrasts a "full" PA
#' (with relatively intact megafauna populations) against an "empty" PA
#' (defaunated or degraded), identified by their WDPA site IDs.
#'
#' @return A \code{tibble} with columns:
#'   \describe{
#'     \item{id}{Integer. Pair identifier (1–10).}
#'     \item{full_name}{Character. Name of the "full" (intact) PA.}
#'     \item{empty_name}{Character. Name of the "empty" (defaunated) PA.}
#'     \item{full_wdpa}{Integer. WDPA site ID of the full PA.}
#'     \item{empty_wdpa}{Integer. WDPA site ID of the empty PA.}
#'     \item{region}{Character. Basin: \code{"Amazon"} or \code{"Congo"}.}
#'   }
#'
#' @details
#' WDPA IDs are taken directly from legacy \code{analyze_pa_pairs.py}:
#' \preformatted{
#'   Amazon:
#'     1. Manu (162)           vs Tapajós (115)
#'     4. Yasuní (164)         vs Cuyabeno (3144)
#'     9. Madidi (18131)       vs Chico Mendes (21018)
#'     5. Tumucumaque (351833) vs Brownsberg (720)
#'     7. Chiribiquete (19984) vs Tinigua (19995)
#'   Congo:
#'     2. Nouabalé-Ndoki (72332) vs Salonga (10906)
#    10. Dzanga-Ndoki (30649)   vs Bushimaie (4369)
#     6. Odzala-Kokoua (1088)   vs Léfini (2266)
#     8. Ivindo (30654)         vs Bili-Uele (1085)
#     3. Lopé (1089)            vs Kundelungu (1084)
#' }
#'
#' @export
define_pa_pairs <- function() {
  tibble::tribble(
    ~id, ~full_name,          ~empty_name,      ~full_wdpa, ~empty_wdpa, ~region,
    1L,  "Manu",              "Tapajós",         162L,       115L,        "Amazon",
    4L,  "Yasuní",            "Cuyabeno",        164L,       3144L,       "Amazon",
    9L,  "Madidi",            "Chico Mendes",    18131L,     21018L,      "Amazon",
    5L,  "Tumucumaque",       "Brownsberg",      351833L,    720L,        "Amazon",
    7L,  "Chiribiquete",      "Tinigua",         19984L,     19995L,      "Amazon",
    2L,  "Nouabalé-Ndoki",    "Salonga",         72332L,     10906L,      "Congo",
    10L, "Dzanga-Ndoki",      "Bushimaie",       30649L,     4369L,       "Congo",
    6L,  "Odzala-Kokoua",     "Léfini",          1088L,      2266L,       "Congo",
    8L,  "Ivindo",            "Bili-Uele",       30654L,     1085L,       "Congo",
    3L,  "Lopé",              "Kundelungu",      1089L,      1084L,       "Congo"
  )
}


# =============================================================================
# RANKED PARK DEFINITIONS
# =============================================================================

#' Define the 20 ranked parks for regional boxplot analysis
#'
#' Returns a tibble with 10 Amazon and 10 Congo parks, ranked from most
#' intact megafauna populations (rank 1) to most defaunated (rank 10).
#' Each park has a search name used for WDPA polygon lookup.
#'
#' @return A \code{tibble} with columns:
#'   \describe{
#'     \item{rank}{Integer. Fauna intactness rank (1 = best).}
#'     \item{name}{Character. Display name of the park.}
#'     \item{search_name}{Character. WDPA search string for polygon lookup.}
#'     \item{status}{Character. Qualitative fauna status label.}
#'     \item{region}{Character. Basin: \code{"Amazon"} or \code{"Congo"}.}
#'   }
#'
#' @details
#' Park lists taken directly from legacy \code{plot_regional_boxplots.py}:
#' \preformatted{
#'   Amazon: Manu, Chiribiquete, Yasuní, Madidi, Tumucumaque,
#'           Jaú, Cuyabeno, Tapajós, Chico Mendes, Tinigua
#'   Congo:  Nouabalé-Ndoki, Lopé, Dzanga-Ndoki, Ivindo,
#'           Odzala-Kokoua, Kahuzi-Biega, Salonga, Bili-Uere,
#'           Léfini, Kundelungu
#' }
#'
#' @export
define_ranked_parks <- function() {
  amazon <- tibble::tribble(
    ~rank, ~name,            ~search_name,       ~status,        ~region,
    1L,    "Manu",           "Manu",             "Pristine",     "Amazon",
    2L,    "Chiribiquete",   "Chiribiquete",     "Pristine",     "Amazon",
    3L,    "Yasuní",         "Yasuní",           "High",         "Amazon",
    4L,    "Madidi",         "Madidi",           "High",         "Amazon",
    5L,    "Tumucumaque",    "Tumucumaque",      "Mod-High",     "Amazon",
    6L,    "Jaú",            "Jau",              "Moderate",     "Amazon",
    7L,    "Cuyabeno",       "Cuyabeno",         "Low-Mod",      "Amazon",
    8L,    "Tapajós",        "Tapajós",          "Low (Silent)", "Amazon",
    9L,    "Chico Mendes",   "Chico Mendes",     "Very Low",     "Amazon",
    10L,   "Tinigua",        "Tinigua",          "Critical",     "Amazon"
  )

  congo <- tibble::tribble(
    ~rank, ~name,              ~search_name,       ~status,        ~region,
    1L,    "Nouabalé-Ndoki",   "Nouabalé-Ndoki",   "Pristine",     "Congo",
    2L,    "Lopé",             "Lopé",             "High",         "Congo",
    3L,    "Dzanga-Ndoki",     "Dzanga",           "High",         "Congo",
    4L,    "Ivindo",           "Ivindo",           "High",         "Congo",
    5L,    "Odzala-Kokoua",    "Odzala",           "High",         "Congo",
    6L,    "Kahuzi-Biega",     "Kahuzi",           "Low-Mod",      "Congo",
    7L,    "Salonga",          "Salonga",          "Low (Empty)",  "Congo",
    8L,    "Bili-Uere",        "Uere",             "Very Low",     "Congo",
    9L,    "Léfini",           "Léfini",           "Critical",     "Congo",
    10L,   "Kundelungu",       "Kundelungu",       "Zero",         "Congo"
  )

  dplyr::bind_rows(amazon, congo)
}


# =============================================================================
# EXTRACTION
# =============================================================================

#' Extract raster values within a PA polygon
#'
#' Wrapper around \code{terra::extract()} for extracting band values
#' from a raster stack within a single PA polygon.
#'
#' @param pa_poly A \code{terra::SpatVector} with one feature (the PA
#'   polygon).
#' @param rast A \code{terra::SpatRaster} to extract from.
#' @param band Character or integer. Band name or index to extract.
#'   If \code{NULL}, all bands are extracted.
#'
#' @return A numeric vector of extracted values (NAs removed), or a
#'   \code{data.frame} if multiple bands are extracted.
#'
#' @export
extract_pa_values <- function(pa_poly, rast, band = NULL) {
  if (!is.null(band)) {
    r <- rast[[band]]
  } else {
    r <- rast
  }
  
  # Align projection of polygon to raster if needed
  if (crs(pa_poly) != crs(r)) {
    pa_poly <- project(pa_poly, crs(r))
  }
  
  # Extract values
  vals <- terra::extract(r, pa_poly, ID = FALSE)
  
  if (is.null(band) || length(names(r)) > 1) {
    # Return data frame, removing rows with all NA
    vals <- vals[rowSums(is.na(vals)) < ncol(vals), , drop = FALSE]
    return(as_tibble(vals))
  } else {
    # Return numeric vector
    vec <- vals[[1]]
    return(vec[!is.na(vec)])
  }
}


# =============================================================================
# PA PAIR ANALYSIS
# =============================================================================

#' Analyse PA pairs: full vs empty comparison
#'
#' For each of the 10 PA pairs (from \code{define_pa_pairs()}), extracts
#' raster values within each polygon, computes raw and (optionally)
#' covariate-adjusted differences, and runs a Welch's t-test.
#'
#' @param stacks A named list of \code{terra::SpatRaster} objects, keyed
#'   by basin (\code{"Congo"}, \code{"Amazon"}).
#' @param pa_polys A \code{terra::SpatVector} of WDPA polygons with a
#'   \code{WDPAID} (or equivalent) column.
#' @param band Character. Band name to analyse (default \code{"uoi"}).
#'
#' @return A \code{tibble} with one row per pair and columns:
#'   \describe{
#'     \item{id, full_name, empty_name, region}{From pair definition.}
#'     \item{full_mean, empty_mean}{Mean band value in each PA.}
#'     \item{raw_diff}{full_mean − empty_mean.}
#'     \item{n_full, n_empty}{Number of extracted pixels.}
#'     \item{t_stat, p_value}{Welch's t-test results.}
#'     \item{ci_low, ci_high}{95\% confidence interval for the difference.}
#'   }
#'
#' @details
#' Mirrors logic from legacy \code{analyze_pa_pairs.py::main()} including
#' the SE approximation: \code{SE = sqrt(s1²/n1 + s2²/n2)}.
#'
#' @export
analyse_pa_pairs <- function(stacks, pa_polys, band = "uoi") {
  pairs_df <- define_pa_pairs()
  results <- list()
  
  # Find WDPA ID column in pa_polys
  id_col <- NULL
  for (col in c("WDPAID", "wdpaid", "WDPA_ID", "pa_id", "PA_ID", "id")) {
    if (col %in% names(pa_polys)) {
      id_col <- col
      break
    }
  }
  if (is.null(id_col)) {
    id_col <- names(pa_polys)[1]
  }
  
  for (i in 1:nrow(pairs_df)) {
    pair <- pairs_df[i, ]
    region_stack <- stacks[[pair$region]]
    
    # Find polygons
    full_poly <- pa_polys[pa_polys[[id_col]] == pair$full_wdpa, ]
    empty_poly <- pa_polys[pa_polys[[id_col]] == pair$empty_wdpa, ]
    
    full_vals <- numeric(0)
    empty_vals <- numeric(0)
    
    if (nrow(full_poly) > 0 && !is.null(region_stack)) {
      full_vals <- tryCatch(
        extract_pa_values(full_poly, region_stack, band),
        error = function(e) numeric(0)
      )
    }
    if (nrow(empty_poly) > 0 && !is.null(region_stack)) {
      empty_vals <- tryCatch(
        extract_pa_values(empty_poly, region_stack, band),
        error = function(e) numeric(0)
      )
    }
    
    full_vals <- full_vals[!is.na(full_vals)]
    empty_vals <- empty_vals[!is.na(empty_vals)]
    
    if (length(full_vals) >= 2 && length(empty_vals) >= 2) {
      tt <- t.test(full_vals, empty_vals, var.equal = FALSE)
      
      results[[i]] <- tibble(
        id = pair$id,
        full_name = pair$full_name,
        empty_name = pair$empty_name,
        region = pair$region,
        full_mean = mean(full_vals),
        empty_mean = mean(empty_vals),
        raw_diff = mean(full_vals) - mean(empty_vals),
        n_full = length(full_vals),
        n_empty = length(empty_vals),
        t_stat = as.numeric(tt$statistic),
        p_value = as.numeric(tt$p.value),
        ci_low = tt$conf.int[1],
        ci_high = tt$conf.int[2]
      )
    } else {
      results[[i]] <- tibble(
        id = pair$id,
        full_name = pair$full_name,
        empty_name = pair$empty_name,
        region = pair$region,
        full_mean = if(length(full_vals) > 0) mean(full_vals) else NA_real_,
        empty_mean = if(length(empty_vals) > 0) mean(empty_vals) else NA_real_,
        raw_diff = NA_real_,
        n_full = length(full_vals),
        n_empty = length(empty_vals),
        t_stat = NA_real_,
        p_value = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      )
    }
  }
  bind_rows(results)
}


# =============================================================================
# RANKED PARK ANALYSIS
# =============================================================================

#' Extract values for ranked parks and return long-form data
#'
#' For each of the 20 ranked parks (from \code{define_ranked_parks()}),
#' extracts raster values within the park polygon and returns a tidy
#' long-form data frame suitable for boxplot visualisation.
#'
#' @param stacks A named list of \code{terra::SpatRaster} objects, keyed
#'   by basin (\code{"Congo"}, \code{"Amazon"}).
#' @param pa_polys A \code{terra::SpatVector} of WDPA polygons.
#' @param band Character. Band name to extract (default \code{"uoi"}).
#'
#' @return A \code{tibble} in long form with columns:
#'   \describe{
#'     \item{rank, name, status, region}{From park definition.}
#'     \item{value}{Extracted raster value for this pixel.}
#'   }
#'
#' @details
#' Mirrors logic from legacy \code{plot_regional_boxplots.py::main()}.
#' Parks with no extractable pixels will appear with zero rows.
#'
#' @export
analyse_ranked_parks <- function(stacks, pa_polys, band = "uoi") {
  parks_df <- define_ranked_parks()
  results <- list()
  
  # Find name column in pa_polys
  name_col <- NULL
  for (col in c("PA_ID", "WDPAID", "NAME", "name", "ORIG_NAME", "orig_name", "DESIG", "desig")) {
    if (col %in% names(pa_polys)) {
      vals <- as.data.frame(pa_polys)[[col]]
      if (any(!is.na(vals))) {
        name_col <- col
        break
      }
    }
  }
  if (is.null(name_col)) {
    name_col <- names(pa_polys)[1]
  }
  
  for (i in 1:nrow(parks_df)) {
    park <- parks_df[i, ]
    region_stack <- stacks[[park$region]]
    
    # Find matching polygons
    match_idx <- which(grepl(park$search_name, pa_polys[[name_col]], ignore.case = TRUE))
    
    if (length(match_idx) > 0 && !is.null(region_stack)) {
      poly_match <- pa_polys[match_idx, ]
      vals <- tryCatch(
        extract_pa_values(poly_match, region_stack, band),
        error = function(e) numeric(0)
      )
      vals <- vals[!is.na(vals)]
      
      if (length(vals) > 0) {
        results[[i]] <- tibble(
          rank = park$rank,
          name = park$name,
          status = park$status,
          region = park$region,
          value = vals
        )
      }
    }
  }
  bind_rows(results)
}
