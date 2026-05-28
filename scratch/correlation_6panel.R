# =============================================================================
# correlation_6panel.R
#
# Generates a 6-panel multipanel figure showing the relationship between GEDI
# UOI and standing biomass sub-components (Total, >100kg, >1000kg) on both
# log1p (top row) and raw/linear (bottom row) scales.
# Performs correlation tests for each panel and adds stats to subtitles.
# =============================================================================

source("code/functions/calibration_helpers.R")
library(ggplot2)
library(dplyr)
library(cowplot)

# 1. Extract scale data
cat("Extracting scale data...\n")
joined_data <- extract_scale_data(5000)

# Helper function to create individual panels
create_panel <- function(df, y_var, y_label, log_scale = TRUE, title_prefix = "") {
  y_vals <- df[[y_var]]
  
  # Calculate correlation statistics
  if (log_scale) {
    pearson_test <- cor.test(df$uoi, log1p(y_vals), method = "pearson")
    spearman_test <- cor.test(df$uoi, log1p(y_vals), method = "spearman", exact = FALSE)
  } else {
    pearson_test <- cor.test(df$uoi, y_vals, method = "pearson")
    spearman_test <- cor.test(df$uoi, y_vals, method = "spearman", exact = FALSE)
  }
  
  # Title and subtitle strings
  panel_title <- sprintf("%s: %s", title_prefix, y_label)
  panel_subtitle <- sprintf("Pearson r = %.3f (p = %.3g)\nSpearman rho = %.3f (p = %.3g)",
                            pearson_test$estimate, pearson_test$p.value,
                            spearman_test$estimate, spearman_test$p.value)
  
  # Create plot
  p <- ggplot(df, aes(x = uoi, y = .data[[y_var]])) +
    geom_point(aes(fill = elephant_present_strict, size = trap_days, alpha = w_temp_cluster, shape = basin),
               color = "black", stroke = 0.2) +
    geom_smooth(method = "lm", color = "#2E7D32", fill = "#2E7D32", alpha = 0.12, linewidth = 0.6) +
    scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22)) +
    scale_fill_manual(values = c("Absent" = "#E06666", "Present" = "#2E7D32")) +
    scale_size_continuous(range = c(1.0, 3.5)) +
    scale_alpha_continuous(range = c(0.25, 1.0)) +
    labs(
      title = panel_title,
      subtitle = panel_subtitle,
      x = "Understory Openness Index (UOI)",
      y = ifelse(log_scale, paste(y_label, "(log1p scale)"), paste(y_label, "(raw scale)"))
    ) +
    theme_bw(base_size = 7.5) +
    theme(
      plot.title = element_text(face = "bold", size = 8),
      plot.subtitle = element_text(size = 6.5, face = "italic", color = "#333333"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.margin = margin(6, 6, 6, 6, "pt")
    )
  
  if (log_scale) {
    p <- p + scale_y_continuous(trans = "log1p")
  }
  
  return(p)
}

# 2. Build individual panels
cat("Building panels...\n")
# Top Row: log1p scale
p1 <- create_panel(joined_data, "B_H_index", "Total Biomass", log_scale = TRUE, title_prefix = "A")
p2 <- create_panel(joined_data, "B_H_gt100", "Biomass >100 kg", log_scale = TRUE, title_prefix = "B")
p3 <- create_panel(joined_data, "B_H_gt1000", "Biomass >1000 kg", log_scale = TRUE, title_prefix = "C")

# Bottom Row: raw/linear scale
p4 <- create_panel(joined_data, "B_H_index", "Total Biomass", log_scale = FALSE, title_prefix = "D")
p5 <- create_panel(joined_data, "B_H_gt100", "Biomass >100 kg", log_scale = FALSE, title_prefix = "E")
p6 <- create_panel(joined_data, "B_H_gt1000", "Biomass >1000 kg", log_scale = FALSE, title_prefix = "F")

# 3. Create a shared legend
p_legend_obj <- ggplot(joined_data) +
  geom_point(aes(x = uoi, y = B_H_index, fill = elephant_present_strict, size = trap_days, alpha = w_temp_cluster, shape = basin), color = "black", stroke = 0.3) +
  scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22), name = "Basin:") +
  scale_size_continuous(name = "Effort (Trap-days):", breaks = c(100, 1000, 5000, 15000), range = c(1.2, 3.5)) +
  scale_alpha_continuous(name = "Temporal Weight:", range = c(0.25, 1.0)) +
  scale_fill_manual(values = c("Absent" = "#E06666", "Present" = "#2E7D32"), name = "Elephant:") +
  theme_bw(base_size = 7.5) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6.5)
  )
shared_legend <- get_legend(p_legend_obj)

# 4. Assemble the grid
cat("Assembling final layout...\n")
grid_plots <- plot_grid(
  p1, p2, p3,
  p4, p5, p6,
  ncol = 3,
  nrow = 2,
  align = "hv"
)

final_layout <- plot_grid(
  grid_plots,
  shared_legend,
  ncol = 1,
  rel_heights = c(1.0, 0.08)
)

# 5. Save final outputs
output_png <- "scratch/correlation_6panel.png"
ggsave(output_png, plot = final_layout, width = 9.5, height = 7.0, dpi = 300)
cat(sprintf("✓ Saved 6-panel plot to %s\n", output_png))

output_pdf <- "scratch/correlation_6panel.pdf"
ggsave(output_pdf, plot = final_layout, width = 9.5, height = 7.0)
cat(sprintf("✓ Saved 6-panel plot to %s\n", output_pdf))

# Copy to brain artifact folder for the user
brain_artifact_dir <- "/home/j/.gemini/antigravity/brain/8f51df52-4604-48e0-9ce8-1c52d1cb241c"
if (file.exists(brain_artifact_dir)) {
  file.copy(output_png, file.path(brain_artifact_dir, "correlation_6panel.png"), overwrite = TRUE)
  file.copy(output_pdf, file.path(brain_artifact_dir, "correlation_6panel.pdf"), overwrite = TRUE)
  cat("✓ Copied outputs to brain artifact folder.\n")
}
