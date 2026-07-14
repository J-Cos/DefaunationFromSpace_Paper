# =============================================================================
# code/functions/framework_helpers.R
#
# Shared helper definitions for Framework 1 and Framework 2 analysis and figure
# generation scripts. Consolidates formula lists, variable maps, and the
# plotmath label formatter that were previously duplicated across 03, 04, and 04b.
# =============================================================================

# --- Framework 1 Candidate Formulas (Beta Regression: UOI ~ Biomass + covariates) ---
framework1_formulas <- list(
  # --- Base Models ---
  "M1: Biomass Only"                       = uoi ~ B_H_index,

  # --- Biomass + Environmental Covariate Set ---
  "M4: Biomass + Elev"                     = uoi ~ B_H_index + elevation,
  "M5: Biomass + Slope"                    = uoi ~ B_H_index + slope,
  "M6: Biomass + HAND"                     = uoi ~ B_H_index + hnd,
  "M7: Biomass + Precip"                   = uoi ~ B_H_index + precip,
  "M8: Biomass + Clay"                     = uoi ~ B_H_index + clay,
  "M9: Biomass + Forest"                   = uoi ~ B_H_index + forest_fraction,

  # --- Elephant Possible Set ---
  "M24_possible: ElephantPossible Only"          = uoi ~ elephant_present_possible,
  "M25_possible: Biomass + ElephantPossible"     = uoi ~ B_H_index + elephant_present_possible,
  "M26_possible: Biomass + ElephantPossible + Elev" = uoi ~ B_H_index + elephant_present_possible + elevation,
  "M27_possible: Biomass + ElephantPossible + Slope" = uoi ~ B_H_index + elephant_present_possible + slope,
  "M28_possible: Biomass + ElephantPossible + HAND" = uoi ~ B_H_index + elephant_present_possible + hnd,
  "M29_possible: Biomass + ElephantPossible + Precip" = uoi ~ B_H_index + elephant_present_possible + precip,
  "M30_possible: Biomass + ElephantPossible + Clay" = uoi ~ B_H_index + elephant_present_possible + clay,
  "M31_possible: Biomass + ElephantPossible + Forest" = uoi ~ B_H_index + elephant_present_possible + forest_fraction,

  # --- Elephant Possible Interaction Set ---
  "M32_possible: Biomass * ElephantPossible"     = uoi ~ B_H_index * elephant_present_possible,
  "M33_possible: Biomass * ElephantPossible + Elev" = uoi ~ B_H_index * elephant_present_possible + elevation,
  "M34_possible: Biomass * ElephantPossible + Slope" = uoi ~ B_H_index * elephant_present_possible + slope,
  "M35_possible: Biomass * ElephantPossible + HAND" = uoi ~ B_H_index * elephant_present_possible + hnd,
  "M36_possible: Biomass * ElephantPossible + Precip" = uoi ~ B_H_index * elephant_present_possible + precip,
  "M37_possible: Biomass * ElephantPossible + Clay" = uoi ~ B_H_index * elephant_present_possible + clay,
  "M38_possible: Biomass * ElephantPossible + Forest" = uoi ~ B_H_index * elephant_present_possible + forest_fraction,

  # --- Elephant Strict Set ---
  "M24_strict: ElephantStrict Only"              = uoi ~ elephant_present_strict,
  "M25_strict: Biomass + ElephantStrict"         = uoi ~ B_H_index + elephant_present_strict,
  "M26_strict: Biomass + ElephantStrict + Elev"     = uoi ~ B_H_index + elephant_present_strict + elevation,
  "M27_strict: Biomass + ElephantStrict + Slope"     = uoi ~ B_H_index + elephant_present_strict + slope,
  "M28_strict: Biomass + ElephantStrict + HAND"      = uoi ~ B_H_index + elephant_present_strict + hnd,
  "M29_strict: Biomass + ElephantStrict + Precip"    = uoi ~ B_H_index + elephant_present_strict + precip,
  "M30_strict: Biomass + ElephantStrict + Clay"      = uoi ~ B_H_index + elephant_present_strict + clay,
  "M31_strict: Biomass + ElephantStrict + Forest"    = uoi ~ B_H_index + elephant_present_strict + forest_fraction,

  # --- Elephant Strict Interaction Set ---
  "M32_strict: Biomass * ElephantStrict"             = uoi ~ B_H_index * elephant_present_strict,
  "M33_strict: Biomass * ElephantStrict + Elev"     = uoi ~ B_H_index * elephant_present_strict + elevation,
  "M34_strict: Biomass * ElephantStrict + Slope"     = uoi ~ B_H_index * elephant_present_strict + slope,
  "M35_strict: Biomass * ElephantStrict + HAND"      = uoi ~ B_H_index * elephant_present_strict + hnd,
  "M36_strict: Biomass * ElephantStrict + Precip"    = uoi ~ B_H_index * elephant_present_strict + precip,
  "M37_strict: Biomass * ElephantStrict + Clay"      = uoi ~ B_H_index * elephant_present_strict + clay,
  "M38_strict: Biomass * ElephantStrict + Forest"    = uoi ~ B_H_index * elephant_present_strict + forest_fraction
)


# --- Framework 2 Candidate Formulas (Tweedie GLM: Biomass ~ UOI + covariates) ---
framework2_formulas <- list(
  # --- Base Models ---
  "M2.1: UOI Only"                           = B_H_index ~ uoi,
  "M2.2p: Elephant Possible Only"            = B_H_index ~ elephant_present_possible,
  "M2.2s: Elephant Strict Only"              = B_H_index ~ elephant_present_strict,
  "M2.3: UOI + Elev"                         = B_H_index ~ uoi + elevation,
  "M2.4: UOI + Slope"                        = B_H_index ~ uoi + slope,
  "M2.5: UOI + HAND"                         = B_H_index ~ uoi + hnd,
  "M2.6: UOI + Precip"                       = B_H_index ~ uoi + precip,
  "M2.7: UOI + Clay"                         = B_H_index ~ uoi + clay,
  "M2.8: UOI + Forest"                       = B_H_index ~ uoi + forest_fraction,
  "M2.9p: UOI + Elephant Possible"           = B_H_index ~ uoi + elephant_present_possible,
  "M2.10p: UOI + Elephant Possible + Elev"   = B_H_index ~ uoi + elephant_present_possible + elevation,
  "M2.11p: UOI + Elephant Possible + Slope"  = B_H_index ~ uoi + elephant_present_possible + slope,
  "M2.12p: UOI + Elephant Possible + HAND"   = B_H_index ~ uoi + elephant_present_possible + hnd,
  "M2.13p: UOI + Elephant Possible + Precip" = B_H_index ~ uoi + elephant_present_possible + precip,
  "M2.14p: UOI + Elephant Possible + Clay"   = B_H_index ~ uoi + elephant_present_possible + clay,
  "M2.15p: UOI + Elephant Possible + Forest" = B_H_index ~ uoi + elephant_present_possible + forest_fraction,
  "M2.9s: UOI + Elephant Strict"             = B_H_index ~ uoi + elephant_present_strict,
  "M2.10s: UOI + Elephant Strict + Elev"     = B_H_index ~ uoi + elephant_present_strict + elevation,
  "M2.11s: UOI + Elephant Strict + Slope"     = B_H_index ~ uoi + elephant_present_strict + slope,
  "M2.12s: UOI + Elephant Strict + HAND"      = B_H_index ~ uoi + elephant_present_strict + hnd,
  "M2.13s: UOI + Elephant Strict + Precip"    = B_H_index ~ uoi + elephant_present_strict + precip,
  "M2.14s: UOI + Elephant Strict + Clay"      = B_H_index ~ uoi + elephant_present_strict + clay,
  "M2.15s: UOI + Elephant Strict + Forest"    = B_H_index ~ uoi + elephant_present_strict + forest_fraction,
  "M2.16p: UOI * Elephant Possible"           = B_H_index ~ uoi * elephant_present_possible,
  "M2.17p: UOI * Elephant Possible + Elev"   = B_H_index ~ uoi * elephant_present_possible + elevation,
  "M2.18p: UOI * Elephant Possible + Slope"  = B_H_index ~ uoi * elephant_present_possible + slope,
  "M2.19p: UOI * Elephant Possible + HAND"   = B_H_index ~ uoi * elephant_present_possible + hnd,
  "M2.20p: UOI * Elephant Possible + Precip" = B_H_index ~ uoi * elephant_present_possible + precip,
  "M2.21p: UOI * Elephant Possible + Clay"   = B_H_index ~ uoi * elephant_present_possible + clay,
  "M2.22p: UOI * Elephant Possible + Forest" = B_H_index ~ uoi * elephant_present_possible + forest_fraction,
  "M2.16s: UOI * Elephant Strict"             = B_H_index ~ uoi * elephant_present_strict,
  "M2.17s: UOI * Elephant Strict + Elev"     = B_H_index ~ uoi * elephant_present_strict + elevation,
  "M2.18s: UOI * Elephant Strict + Slope"     = B_H_index ~ uoi * elephant_present_strict + slope,
  "M2.19s: UOI * Elephant Strict + HAND"      = B_H_index ~ uoi * elephant_present_strict + hnd,
  "M2.20s: UOI * Elephant Strict + Precip"    = B_H_index ~ uoi * elephant_present_strict + precip,
  "M2.21s: UOI * Elephant Strict + Clay"      = B_H_index ~ uoi * elephant_present_strict + clay,
  "M2.22s: UOI * Elephant Strict + Forest"    = B_H_index ~ uoi * elephant_present_strict + forest_fraction
)


# --- Variable Maps for plotmath label formatting ---

#' Variable map for Framework 1 (Beta Regression: UOI ~ Biomass + covariates)
VAR_MAP_FRAMEWORK1 <- list(
  "Biomass"          = "B_H_index",
  "Basin"            = c("basinCongo", "basinSE_Asia"),
  "ElephantPossible" = "elephant_present_possiblePresent",
  "ElephantStrict"   = "elephant_present_strictPresent",
  "Elephant"         = "elephant_presentPresent",
  "UOI"              = "uoi",
  "Elevation"        = "elevation",
  "Elev"             = "elevation",
  "Slope"            = "slope",
  "HAND"             = "hnd",
  "Precipitation"    = "precip",
  "Precip"           = "precip",
  "Clay"             = "clay",
  "Forest"           = "forest_fraction",
  "UOI:Basin"        = "uoi:basinCongo",
  "Biomass:Basin"    = c("B_H_index:basinCongo", "B_H_index:basinSE_Asia")
)

#' Variable map for Framework 2 (Tweedie GLM: Biomass ~ UOI + covariates)
VAR_MAP_FRAMEWORK2 <- list(
  "ElephantPossible" = "elephant_present_possiblePresent",
  "ElephantStrict"   = "elephant_present_strictPresent",
  "Elephant"         = c("elephant_presentPresent", "elephant_present_possiblePresent", "elephant_present_strictPresent"),
  "Basin"            = c("basinCongo", "basinSE_Asia"),
  "UOI"              = "uoi",
  "Elevation"        = "elevation",
  "Elev"             = "elevation",
  "Slope"            = "slope",
  "HAND"             = "hnd",
  "Precipitation"    = "precip",
  "Precip"           = "precip",
  "Clay"             = "clay",
  "Forest"           = "forest_fraction",
  "UOI:Elephant"     = c("uoi:elephant_presentPresent", "uoi:elephant_present_possiblePresent", "uoi:elephant_present_strictPresent"),
  "UOI:Basin"        = c("uoi:basinCongo", "uoi:basinSE_Asia")
)


#' Format model labels with plotmath bolding of significant variables
#'
#' Creates plotmath expressions for y-axis labels in model selection bar plots.
#' Variables that are statistically significant (p < 0.05) are rendered in bold.
#'
#' @param model_name Character. The model name (e.g., "M2.1: UOI Only")
#' @param model_obj A fitted GAM model object (used to extract p-values)
#' @param var_map Named list mapping display tokens to coefficient names in the model
#' @param prefix_pattern Regex pattern to strip the model number prefix
#'   (default strips both "M1:", "M24_strict:", "M2.1:", etc.)
#'
#' @return A character string containing a plotmath expression
format_model_label <- function(model_name, model_obj, var_map,
                               prefix_pattern = "^M[0-9\\.]+[a-z_]*: ") {
  clean_name <- gsub(prefix_pattern, "", model_name)
  clean_name <- gsub(" (Shared)", "", clean_name, fixed = TRUE)

  tokens <- strsplit(clean_name, "\\s+")[[1]]
  tokens <- tokens[tokens != ""]

  p_table <- summary(model_obj)$p.table

  plotmath_tokens <- sapply(tokens, function(tok) {
    if (tok %in% c("+", "*", ":")) return(sprintf("plain(\" %s \")", tok))
    if (tok == "Only") return(sprintf("plain(\" %s\")", tok))

    matched_terms <- var_map[[tok]]
    if (!is.null(matched_terms)) {
      is_sig <- FALSE
      for (term in matched_terms) {
        if (term %in% rownames(p_table)) {
          p_val <- p_table[term, ncol(p_table)]
          if (!is.na(p_val) && p_val < 0.05) {
            is_sig <- TRUE
            break
          }
        }
      }
      return(ifelse(is_sig, sprintf("bold(\"%s\")", tok), sprintf("plain(\"%s\")", tok)))
    } else {
      return(sprintf("plain(\"%s\")", tok))
    }
  })

  paste(plotmath_tokens, collapse = " * ")
}
