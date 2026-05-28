# =============================================================================
# theme_pnas.R
#
# PNAS-style ggplot2 theme, colour palettes, and figure export helpers.
# Designed for the Defaunation-from-Space manuscript.
# =============================================================================

library(ggplot2)
library(scales)

# -----------------------------------------------------------------------------
# Font registration (Conditional)
# -----------------------------------------------------------------------------
has_showtext <- suppressWarnings(requireNamespace("showtext", quietly = TRUE)) &&
                suppressWarnings(requireNamespace("sysfonts", quietly = TRUE))

base_family <- "sans"
if (has_showtext) {
  library(showtext)
  library(sysfonts)
  tryCatch({
    font_add_google("Roboto Condensed", "Roboto Condensed")
    showtext_auto()
    base_family <- "Roboto Condensed"
    cat("Roboto Condensed registered successfully via showtext.\n")
  }, error = function(e) {
    cat("Warning: Failed to load Google Font Roboto Condensed, falling back to sans.\n")
  })
} else {
  cat("Warning: showtext package not found. Using system default sans font.\n")
}

# =============================================================================
# THEME
# =============================================================================

#' PNAS-style ggplot2 theme
#'
#' A clean, compact theme following PNAS formatting guidelines.
#' Uses Roboto Condensed via showtext, removes gridlines, and sets
#' compact margins suitable for multi-panel figures.
#'
#' @param base_size Numeric. Base font size in points (default 8, PNAS minimum).
#'
#' @return A \code{ggplot2::theme} object.
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_pnas()
#'
#' @export
theme_pnas <- function(base_size = 8) {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      line = element_line(colour = "black", linewidth = 0.25, linetype = 1, lineend = "butt"),
      rect = element_rect(fill = "white", colour = "black", linewidth = 0.25, linetype = 1),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      text = element_text(family = base_family, face = "plain", colour = "black", size = base_size, lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug = FALSE),
      axis.line = element_line(colour = "black", linewidth = 0.25),
      axis.ticks = element_line(colour = "black", linewidth = 0.25),
      axis.ticks.length = unit(0.08, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.25),
      legend.position = "bottom",
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      legend.key.size = unit(0.3, "cm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0),
      plot.caption = element_text(size = base_size - 3, hjust = 1),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
      strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.25),
      strip.text = element_text(size = base_size, face = "bold", margin = margin(t = 2, r = 2, b = 2, l = 2))
    )
}

# =============================================================================
# COLOUR PALETTES
# =============================================================================

#' Basin palette
#'
#' Named colour vector for the two study basins.
#' @export
pal_basin <- c(
  Congo  = "#1B5E20",
  Amazon = "#E65100",
  SE_Asia = "#0D47A1"
)

#' Signal palette
#'
#' Named colour vector for the two primary remote-sensing signals.
#' @export
pal_signal <- c(
  UOI  = "#1B5E20",
  FRIP = "#BF360C"
)

#' Trend palette
#'
#' Named colour vector for Mann-Kendall trend categories.
#' @export
pal_trend <- c(
  Decreasing = "#388E3C",
  Increasing = "#F57C00",
  None       = "grey70"
)

#' Protection status palette
#'
#' Named colour vector for protected / unprotected classification.
#' @export
pal_protect <- c(
  Protected   = "#1565C0",
  Unprotected = "#B71C1C"
)

# =============================================================================
# DIVERGING FILL SCALES
# =============================================================================

#' Diverging fill scale for FRIP (blue → grey95 → red)
#'
#' Continuous fill scale centred on zero using \code{scales::muted()} endpoints.
#' Suitable for FRIP correlation values spanning negative to positive.
#'
#' @param ... Additional arguments passed to
#'   \code{ggplot2::scale_fill_gradient2}.
#'
#' @return A ggplot2 scale object.
#' @export
scale_fill_frip <- function(...) {
  scale_fill_gradient2(
    low = scales::muted("blue"),
    mid = "grey95",
    high = scales::muted("red"),
    midpoint = 0,
    ...
  )
}

#' Diverging fill scale for Mann-Kendall tau (green → grey95 → orange)
#'
#' Continuous fill scale centred on zero.
#' Suitable for tau values indicating temporal trends.
#'
#' @param ... Additional arguments passed to
#'   \code{ggplot2::scale_fill_gradient2}.
#'
#' @return A ggplot2 scale object.
#' @export
scale_fill_tau <- function(...) {
  scale_fill_gradient2(
    low = "#388E3C",
    mid = "grey95",
    high = "#F57C00",
    midpoint = 0,
    ...
  )
}

#' Sequential fill scale for UOI (viridis option D)
#'
#' Continuous fill scale using the viridis-D colour map.
#' Suitable for Understory Openness Index values.
#'
#' @param ... Additional arguments passed to
#'   \code{ggplot2::scale_fill_viridis_c}.
#'
#' @return A ggplot2 scale object.
#' @export
scale_fill_uoi <- function(...) {
  scale_fill_viridis_c(option = "plasma", ...)
}

# =============================================================================
# PNAS FIGURE DIMENSIONS & EXPORT
# =============================================================================

#' PNAS single-column width (cm)
#' @export
PNAS_SINGLE_W_CM <- 8.7

#' PNAS double-column width (cm)
#' @export
PNAS_DOUBLE_W_CM <- 17.8

#' PNAS maximum figure height (cm)
#' @export
PNAS_MAX_H_CM <- 23.0

#' PNAS minimum resolution (dpi)
#' @export
PNAS_DPI <- 600

#' Save a figure to PNAS specifications
#'
#' Wrapper around \code{ggplot2::ggsave} that enforces PNAS column widths,
#' resolution, and file format defaults.
#'
#' @param plot A ggplot2 plot object.
#' @param filename Character. Output file path (e.g. \code{"figures/fig1.pdf"}).
#' @param type Character. One of \code{"single"} (8.7 cm wide) or
#'   \code{"double"} (17.8 cm wide). Default \code{"single"}.
#' @param height_cm Numeric. Figure height in centimetres.
#'   Defaults to width (square) if \code{NULL}.
#' @param ... Additional arguments passed to \code{ggplot2::ggsave}.
#'
#' @return Invisibly returns \code{filename} (side-effect: writes file).
#'. @export
save_pnas <- function(plot, filename, type = c("single", "double"),
                      height_cm = NULL, ...) {
  type <- match.arg(type)
  width_cm <- if (type == "single") PNAS_SINGLE_W_CM else PNAS_DOUBLE_W_CM
  if (is.null(height_cm)) {
    height_cm <- width_cm
  }
  
  # Limit height to max
  if (height_cm > PNAS_MAX_H_CM) {
    warning("Requested height exceeds PNAS maximum height. Capping at ", PNAS_MAX_H_CM, " cm.")
    height_cm <- PNAS_MAX_H_CM
  }
  
  # Convert cm to inches
  w_in <- width_cm / 2.54
  h_in <- height_cm / 2.54
  
  ggsave(
    filename = filename,
    plot = plot,
    width = w_in,
    height = h_in,
    dpi = PNAS_DPI,
    bg = "white",
    ...
  )
  
  invisible(filename)
}

#' Plot Model Selection Bar Chart
#'
#' A unified, highly DRY plotting function to render model selection bar charts
#' for both Framework 1 and Framework 2 (LOBO and AIC variants).
#'
#' @param plot_df Data frame containing the models to display.
#' @param x_var Character. Name of the column to map to the x-axis (e.g., "R2", "Full_DevExpl", "delta_AIC").
#' @param fill_var Character. Name of the column to map to the fill color (e.g., "delta_AIC", "OOS_MAE_log", "Full_DevExpl").
#' @param fill_label Character. Title for the fill color legend.
#' @param x_label Character. Title for the x-axis.
#' @param plot_title Character. Title for the panel.
#' @param colors Character vector. Colors for the gradient (default is the standard deep blue to red palette).
#' @param parsed_labels Parsed plotmath expressions for the y-axis labels.
#'
#' @return A ggplot2 plot object.
#' @export
plot_model_selection_bars <- function(plot_df, x_var, fill_var, fill_label, x_label, plot_title,
                                      colors = c("#0D47A1", "#1976D2", "#64B5F6", "#FFA726", "#F57C00", "#D84315"),
                                      parsed_labels = NULL) {
  p <- ggplot(plot_df, aes(x = .data[[x_var]], y = CleanName, fill = .data[[fill_var]])) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.2)
  
  if (!is.null(parsed_labels)) {
    p <- p + scale_y_discrete(labels = parsed_labels)
  }
  
  p <- p +
    scale_fill_gradientn(colors = colors, name = fill_label) +
    labs(
      title = plot_title,
      x = x_label,
      y = "Model Formulation"
    ) +
    theme_pnas(base_size = 7.5) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 8.5, margin = margin(b = 6, t = 4)),
      axis.text.y = element_text(size = 6.5),
      axis.title.x = element_text(margin = margin(t = 4)),
      axis.title.y = element_text(margin = margin(r = 4)),
      legend.title = element_text(size = 6.0, face = "bold"),
      legend.text = element_text(size = 5.0),
      legend.key.width = unit(0.15, "cm"),
      legend.key.height = unit(0.25, "cm"),
      legend.margin = margin(l = 2, r = 2, unit = "pt"),
      plot.margin = margin(t = 6, r = 4, b = 6, l = 4, unit = "pt")
    )
  
  return(p)
}

