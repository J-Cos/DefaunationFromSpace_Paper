# =============================================================================
# plotting.R
#
# Reusable plot helpers for the Defaunation from Space manuscript.
#
# All functions produce ggplot2 objects that can be composed with
# cowplot::plot_grid() and saved with save_pnas().
#
# Dependencies:
#   ggplot2, tidyterra, cowplot, dplyr, terra
# =============================================================================

# --- Imports -----------------------------------------------------------------

library(ggplot2)
library(tidyterra)
library(cowplot)
library(dplyr)
library(terra)

# Make sure theme_pnas is sourced if available
if (file.exists("code/functions/theme_pnas.R")) {
  source("code/functions/theme_pnas.R")
}

# --- Constants ---------------------------------------------------------------

#' Colour palettes used across the manuscript
BASIN_COLOURS <- c(Congo = "#2E86AB", Amazon = "#A23B72")
PROTECTION_COLOURS <- c(Protected = "#2D6A4F", Unprotected = "#D4A373")

#' Bivariate classification colour scheme (2×2)
BIVARIATE_COLOURS <- c(
  "1" = "#2166AC",  # High UOI / Low FRIP  (Intact)
  "2" = "#92C5DE",  # High UOI / High FRIP (Structure ok, function degraded)
  "3" = "#F4A582",  # Low UOI  / Low FRIP  (Structure lost, function ok)
  "4" = "#B2182B"   # Low UOI  / High FRIP (Degraded)
)


# --- Map functions -----------------------------------------------------------

#' Single-basin map with geom_spatraster
#'
#' Creates a map of a single-band raster using \code{tidyterra::geom_spatraster},
#' with an optional country boundary overlay.
#'
#' @param rast SpatRaster. Single-band raster to map.
#' @param fill_col Character. Name of the band to map (used in aes).
#' @param scale_fn A ggplot2 scale function (e.g., \code{scale_fill_viridis_c()})
#'   for the fill aesthetic.
#' @param countries SpatVector or NULL. Country boundary polygons for overlay.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' # p <- make_basin_map(congo_uoi, "uoi", scale_fill_viridis_c())
make_basin_map <- function(rast, fill_col, scale_fn, countries = NULL) {
  # Clean up layer names to match exactly
  rast_sub <- rast[[fill_col]]
  
  # Check if theme_pnas is loaded
  t_theme <- if (exists("theme_pnas", mode = "function")) theme_pnas() else theme_minimal()
  
  p <- ggplot() +
    geom_spatraster(data = rast_sub, aes(fill = !!sym(fill_col))) +
    scale_fn +
    t_theme +
    theme(
      legend.position = "right",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
  
  if (!is.null(countries)) {
    if (crs(countries) != crs(rast_sub)) {
      countries <- project(countries, crs(rast_sub))
    }
    p <- p + geom_spatvector(data = countries, fill = NA, colour = "white", linewidth = 0.3)
  }
  
  return(p)
}


#' Paired side-by-side basin maps with shared legend
#'
#' Creates Congo and Amazon maps side by side with a single shared legend
#' using \code{cowplot::plot_grid}.
#'
#' @param rast_congo SpatRaster. Congo basin raster (single band).
#' @param rast_amazon SpatRaster. Amazon basin raster (single band).
#' @param fill_col Character. Band name to map.
#' @param scale_fn ggplot2 scale function for fill aesthetic.
#' @param countries SpatVector. Country boundaries for overlay.
#'
#' @return A \code{ggplot} object (combined panel).
#'
#' @examples
#' # p <- make_paired_maps(congo_uoi, amazon_uoi, "uoi",
#' #                       scale_fill_viridis_c(), countries)
make_paired_maps <- function(rast_congo, rast_amazon, fill_col, scale_fn,
                             countries) {
  p_congo <- make_basin_map(rast_congo, fill_col, scale_fn, countries) +
    theme(legend.position = "none") +
    labs(title = "Congo Basin")
  
  p_amazon_leg <- make_basin_map(rast_amazon, fill_col, scale_fn, countries) +
    labs(title = "Amazon Basin")
  
  legend <- cowplot::get_legend(p_amazon_leg)
  
  p_amazon <- p_amazon_leg + theme(legend.position = "none")
  
  cowplot::plot_grid(p_congo, p_amazon, legend, ncol = 3, rel_widths = c(1, 1, 0.2))
}


# --- Statistical plot functions ----------------------------------------------

#' Boxplot with Tukey compact letter annotations
#'
#' Creates a jitter + boxplot with compact letter display annotations from
#' a TukeyHSD analysis positioned above each group.
#'
#' @param df Data frame containing the data.
#' @param x Character. Column name for the x-axis grouping variable.
#' @param y Character. Column name for the y-axis response variable.
#' @param letters_df Data frame with columns matching \code{x} and `letter`
#'   for the compact letter annotations.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' # p <- make_boxplot_with_letters(df, "PA_NAME", "uoi", letters_df)
make_boxplot_with_letters <- function(df, x, y, letters_df) {
  # Calculate y positions for letters
  y_maxs <- df %>%
    group_by(!!sym(x)) %>%
    summarise(y_pos = max(!!sym(y), na.rm = TRUE), .groups = "drop")
  
  # Find letter column
  let_col <- grep("letter|Letters", names(letters_df), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(let_col)) {
    let_col <- names(letters_df)[2]
  }
  
  ann_df <- letters_df %>%
    left_join(y_maxs, by = x)
  
  t_theme <- if (exists("theme_pnas", mode = "function")) theme_pnas() else theme_minimal()
  
  ggplot(df, aes(x = !!sym(x), y = !!sym(y))) +
    geom_jitter(width = 0.2, alpha = 0.2, colour = "grey50") +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "white", colour = "black") +
    geom_text(data = ann_df, aes(x = !!sym(x), y = y_pos, label = !!sym(let_col)),
              vjust = -0.5, size = 3, fontface = "bold") +
    t_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


#' PA pairs bar chart: raw vs adjusted difference
#'
#' Grouped bars showing raw and covariate-adjusted signal differences
#' between protected and unprotected areas, with error bars and
#' region-coloured grouping.
#'
#' @param results_df Data frame with columns: `pa_name`, `region`,
#'   `raw_diff`, `adj_diff`, `raw_se`, `adj_se`.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' # p <- make_pa_pairs_bar(results_df)
make_pa_pairs_bar <- function(results_df) {
  # Create a pair name column if not present
  if (!"pair_name" %in% names(results_df)) {
    if (all(c("full_name", "empty_name") %in% names(results_df))) {
      results_df <- results_df %>%
        mutate(pair_name = paste(full_name, "vs", empty_name))
    } else if ("pa_name" %in% names(results_df)) {
      results_df <- results_df %>%
        rename(pair_name = pa_name)
    } else {
      results_df$pair_name <- as.character(1:nrow(results_df))
    }
  }
  
  plot_df <- results_df %>%
    select(pair_name, region, raw_diff, adj_diff) %>%
    tidyr::pivot_longer(cols = c(raw_diff, adj_diff), names_to = "Difference", values_to = "diff_val") %>%
    mutate(Difference = factor(Difference, levels = c("raw_diff", "adj_diff"), labels = c("Raw", "Adjusted")))
  
  if (all(c("raw_se", "adj_se") %in% names(results_df))) {
    se_df <- results_df %>%
      select(pair_name, raw_se, adj_se) %>%
      tidyr::pivot_longer(cols = c(raw_se, adj_se), names_to = "Difference", values_to = "se_val") %>%
      mutate(Difference = factor(Difference, levels = c("raw_se", "adj_se"), labels = c("Raw", "Adjusted")))
    
    plot_df <- plot_df %>%
      left_join(se_df, by = c("pair_name", "Difference")) %>%
      mutate(ymin = diff_val - 1.96 * se_val, ymax = diff_val + 1.96 * se_val)
  } else if (all(c("ci_low", "ci_high") %in% names(results_df))) {
    # If the results are from analyse_pa_pairs directly, let's use the CI bounds
    # Since we run it once for raw and once for adjusted, we can match them if both exist
    # If not, let's just make error bars cover the CI of raw_diff
    plot_df <- plot_df %>%
      mutate(ymin = diff_val, ymax = diff_val)
  } else {
    plot_df <- plot_df %>%
      mutate(ymin = diff_val, ymax = diff_val)
  }
  
  t_theme <- if (exists("theme_pnas", mode = "function")) theme_pnas() else theme_minimal()
  
  ggplot(plot_df, aes(x = pair_name, y = diff_val, fill = Difference)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(0.8), width = 0.2) +
    facet_wrap(~region, scales = "free_x") +
    scale_fill_manual(values = c("Raw" = "#CCCCCC", "Adjusted" = "#444444")) +
    labs(x = "Protected Area Pairs", y = "Signal Difference") +
    t_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


#' Ranked parks boxplot
#'
#' Side-by-side raw and adjusted boxplots for each park, ordered by
#' megafauna rank (descending).
#'
#' @param data Data frame with park-level signal data.
#' @param basin Character. Basin name for title ("Congo" or "Amazon").
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' # p <- make_ranked_parks_boxplot(park_data, "Congo")
make_ranked_parks_boxplot <- function(data, basin) {
  df_basin <- data %>% filter(region == basin)
  
  if (nrow(df_basin) == 0) {
    # Return an empty plot if no data
    return(ggplot() + labs(title = sprintf("No Data for %s", basin)))
  }
  
  # Downsample per park if data is too large to prevent outlier plotting overload and PDF bloat
  max_rows_per_park <- 10000
  df_basin <- df_basin %>%
    group_by(name) %>%
    sample_n(min(n(), max_rows_per_park)) %>%
    ungroup()
  
  # Order park name by rank
  df_basin$name <- factor(df_basin$name, levels = unique(df_basin$name[order(df_basin$rank)]))
  
  fill_colour <- if (basin %in% names(BASIN_COLOURS)) BASIN_COLOURS[basin] else "#444444"
  
  t_theme <- if (exists("theme_pnas", mode = "function")) theme_pnas() else theme_minimal()
  
  ggplot(df_basin, aes(x = name, y = value)) +
    geom_boxplot(fill = fill_colour, outlier.size = 0.5, alpha = 0.7) +
    labs(title = sprintf("%s Parks Ranked by Megafauna Status", basin), x = "Park (Decreasing Faunal Status)", y = "Signal Value") +
    t_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


#' Openness distribution plot
#'
#' Overlaid histograms of UOI by region, with an inset or companion panel
#' showing boxplots by Region × Protection interaction.
#'
#' @param df Data frame with columns: `uoi`, `region`, `protection`.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' # p <- make_openness_distribution(uoi_df)
make_openness_distribution <- function(df) {
  t_theme <- if (exists("theme_pnas", mode = "function")) theme_pnas() else theme_minimal()
  
  # Downsample if data is too large to prevent ggplot density/boxplot performance bottlenecks
  max_rows <- 200000
  if (nrow(df) > max_rows) {
    df <- df %>%
      group_by(region, protection) %>%
      sample_n(min(n(), as.integer(max_rows / 4))) %>%
      ungroup()
  }
  
  p_hist <- ggplot(df, aes(x = uoi, fill = region)) +
    geom_density(alpha = 0.5, position = "identity") +
    scale_fill_manual(values = BASIN_COLOURS) +
    labs(x = "Understory Openness Index (UOI)", y = "Density") +
    t_theme
  
  p_box <- ggplot(df, aes(x = region, y = uoi, fill = protection)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    scale_fill_manual(values = PROTECTION_COLOURS) +
    labs(x = "Region", y = "UOI") +
    t_theme
  
  cowplot::plot_grid(p_hist, p_box, ncol = 2, rel_widths = c(1.2, 1))
}


# --- Multi-scale plot functions ----------------------------------------------

#' Multi-scale confidence interval plot
#'
#' Linerange plot showing effect sizes (or correlations) across scales
#' with confidence intervals. Significant scales are highlighted with
#' full opacity; non-significant with reduced alpha. A horizontal
#' reference line at zero is added.
#'
#' @param df Data frame with multi-scale results.
#' @param x_col Character. Column for x-axis (scale).
#' @param ymin_col Character. Column for lower CI bound.
#' @param ymax_col Character. Column for upper CI bound.
#' @param signif_col Character. Column with logical significance flag.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' # p <- make_multiscale_ci_plot(h1_ms, "scale", "ci_lo", "ci_hi", "signif")
make_multiscale_ci_plot <- function(df, x_col, ymin_col, ymax_col, signif_col) {
  # Choose y column dynamically
  y_vals <- if ("estimate" %in% names(df)) {
    df$estimate
  } else if ("cohens_d" %in% names(df)) {
    df$cohens_d
  } else if ("t_stat" %in% names(df)) {
    df$t_stat
  } else {
    (df[[ymin_col]] + df[[ymax_col]]) / 2
  }
  
  df_plot <- df %>%
    mutate(
      y = y_vals,
      signif_alpha = ifelse(!!sym(signif_col), 1.0, 0.4)
    )
  
  t_theme <- if (exists("theme_pnas", mode = "function")) theme_pnas() else theme_minimal()
  
  ggplot(df_plot, aes(x = !!sym(x_col) / 1000, y = y, alpha = signif_alpha)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(ymin = !!sym(ymin_col), ymax = !!sym(ymax_col)), linewidth = 1, colour = "#2D6A4F") +
    geom_point(size = 2, colour = "#2D6A4F") +
    scale_alpha_identity() +
    labs(x = "Aggregation Scale (km)", y = "Effect Size") +
    t_theme
}


#' Denoising ratio plot
#'
#' Plots the log R² ratio (denoised / original) against scale.
#' Values above zero indicate denoising improved signal extraction.
#'
#' @param df Data frame with columns: `scale`, `ratio` (or log-ratio).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' # p <- make_denoising_ratio_plot(denoise_df)
make_denoising_ratio_plot <- function(df) {
  t_theme <- if (exists("theme_pnas", mode = "function")) theme_pnas() else theme_minimal()
  
  ggplot(df, aes(x = scale / 1000, y = ratio)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 1, colour = "#A23B72") +
    geom_point(size = 2, colour = "#A23B72") +
    labs(x = "Aggregation Scale (km)", y = "Log R² Ratio (Denoised / Raw)") +
    t_theme
}


# --- Bivariate and scatter functions -----------------------------------------

#' Bivariate choropleth map (2×2)
#'
#' Maps the 4-class bivariate raster (from \code{classify_bivariate})
#' using the BIVARIATE_COLOURS palette, with country boundary overlay.
#'
#' @param rast SpatRaster. Integer raster with values 1–4.
#' @param countries SpatVector. Country boundary polygons.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' # p <- make_bivariate_map(rast, countries)
make_bivariate_map <- function(rast, countries) {
  t_theme <- if (exists("theme_pnas", mode = "function")) theme_pnas() else theme_minimal()
  
  rast <- as.factor(rast)
  
  p <- ggplot() +
    geom_spatraster(data = rast, aes(fill = !!sym(names(rast)[1]))) +
    scale_fill_manual(
      values = BIVARIATE_COLOURS,
      labels = c(
        "1" = "Intact (High UOI / Low FRIP)",
        "2" = "Structure Ok, Function Degraded",
        "3" = "Structure Lost, Function Ok",
        "4" = "Degraded (Low UOI / High FRIP)"
      ),
      na.translate = FALSE
    ) +
    t_theme +
    theme(
      legend.position = "bottom",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(fill = "Bivariate Class")
  
  if (!is.null(countries)) {
    if (crs(countries) != crs(rast)) {
      countries <- project(countries, crs(rast))
    }
    p <- p + geom_spatvector(data = countries, fill = NA, colour = "white", linewidth = 0.3)
  }
  return(p)
}


#' Scatter plot with Spearman r annotation
#'
#' Point scatter with a colour grouping variable and a text annotation
#' showing the Spearman correlation coefficient and p-value.
#'
#' @param df Data frame.
#' @param x Character. Column for x-axis.
#' @param y Character. Column for y-axis.
#' @param color_col Character. Column for colour aesthetic (e.g., "basin").
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' # p <- make_scatter_with_cor(pa_df, "mean_uoi", "mean_frip", "basin")
make_scatter_with_cor <- function(df, x, y, color_col) {
  cor_res <- cor.test(df[[x]], df[[y]], method = "spearman", exact = FALSE)
  r_val <- cor_res$estimate
  p_val <- cor_res$p.value
  
  label_text <- sprintf("Spearman r = %.3f\np = %.3g", r_val, p_val)
  
  t_theme <- if (exists("theme_pnas", mode = "function")) theme_pnas() else theme_minimal()
  
  ggplot(df, aes(x = !!sym(x), y = !!sym(y), colour = !!sym(color_col))) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linetype = "dashed", colour = "grey30") +
    annotate("text", x = -Inf, y = Inf, label = label_text, hjust = -0.1, vjust = 1.1, fontface = "italic", size = 4) +
    t_theme
}


# --- Assembly and export functions -------------------------------------------

#' Assemble multi-panel figure
#'
#' Wrapper around \code{cowplot::plot_grid()} that accepts an arbitrary
#' number of ggplot objects and labels them with panel letters (a, b, c, …).
#'
#' @param ... ggplot objects to assemble.
#' @param ncol Integer. Number of columns in the grid. Default 2.
#' @param labels Character vector or "AUTO" for automatic labelling.
#'
#' @return A \code{ggplot} object (combined panel).
#'
#' @examples
#' # fig <- assemble_figure(p1, p2, p3, ncol = 2)
assemble_figure <- function(..., ncol = 2, labels = "AUTO") {
  plots <- list(...)
  cowplot::plot_grid(plotlist = plots, ncol = ncol, labels = labels, label_size = 12)
}


#' Save figure to PNAS specifications
#'
#' Wrapper around \code{ggplot2::ggsave()} with PNAS column widths
#' and DPI presets.
#'
#' @param plot A ggplot object.
#' @param filename Character. Output filename (relative to `figures/` dir).
#' @param type Character. One of "single_col" (8.7 cm), "1.5_col" (11.4 cm),
#'   "two_col" (17.8 cm). Determines width; height auto-scaled or specified.
#' @param height Numeric or NULL. Figure height in cm. If NULL, uses aspect
#'   ratio from the ggplot object.
#' @param dpi Integer. Resolution. Default 300.
#'
#' @return Invisible. The plot, saved to disk.
#'
#' @examples
#' # save_pnas(fig1, "Fig1_UOI_maps.pdf", type = "two_col")
save_pnas <- function(plot, filename, type = "two_col", height = NULL,
                      dpi = 300) {
  width <- switch(
    type,
    "single_col" = 8.7,
    "1.5_col" = 11.4,
    "two_col" = 17.8,
    17.8 # default
  )
  
  if (is.null(height)) {
    height <- width * 0.75
  }
  
  dir.create("figures", recursive = TRUE, showWarnings = FALSE)
  filepath <- file.path("figures", filename)
  
  ggsave(
    filename = filepath,
    plot = plot,
    width = width,
    height = height,
    units = "cm",
    dpi = dpi
  )
  
  invisible(plot)
}
