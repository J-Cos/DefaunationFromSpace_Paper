# =============================================================================
# correlation_uoi_biomass.R
#
# Creates a simple correlation scatter plot and performs statistical correlation
# tests between GEDI Understory Openness Index (UOI) and Standing Mammal Biomass.
# =============================================================================

source("code/functions/calibration_helpers.R")
library(ggplot2)
library(dplyr)

# 1. Extract scale data
cat("Extracting scale data...\n")
joined_data <- extract_scale_data(5000)

# 2. Perform correlation tests
cat("\n--- Statistical Correlation Tests ---\n")
# Pearson correlation on raw scale
pearson_raw <- cor.test(joined_data$uoi, joined_data$B_H_index, method = "pearson")
cat(sprintf("Pearson (raw scale): r = %.3f, p-value = %e\n", pearson_raw$estimate, pearson_raw$p.value))

# Spearman correlation (rank-based, robust to skewness)
spearman_raw <- cor.test(joined_data$uoi, joined_data$B_H_index, method = "spearman", exact = FALSE)
cat(sprintf("Spearman (rank-based): rho = %.3f, p-value = %e\n", spearman_raw$estimate, spearman_raw$p.value))

# Pearson correlation on log1p scale
pearson_log <- cor.test(joined_data$uoi, log1p(joined_data$B_H_index), method = "pearson")
cat(sprintf("Pearson (log1p scale): r = %.3f, p-value = %e\n", pearson_log$estimate, pearson_log$p.value))

# 3. Create correlation scatter plot
cat("\nGenerating scatter plot...\n")
p <- ggplot(joined_data, aes(x = uoi, y = B_H_index)) +
  geom_point(aes(fill = elephant_present_strict, size = trap_days, alpha = w_temp_cluster, shape = basin),
             color = "black", stroke = 0.3) +
  geom_smooth(method = "lm", color = "#2E7D32", fill = "#2E7D32", alpha = 0.15, linewidth = 0.8) +
  scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22), name = "Basin") +
  scale_fill_manual(values = c("Absent" = "#E06666", "Present" = "#2E7D32"), name = "Elephant Presence") +
  scale_size_continuous(name = "Effort (Trap-days)", range = c(1.5, 4.5)) +
  scale_alpha_continuous(name = "Temporal Weight", range = c(0.25, 1.0)) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 5000, 15000)) +
  labs(
    title = "Correlation: GEDI Openness (UOI) vs. Standing Mammal Biomass",
    subtitle = sprintf("Spearman rho = %.3f (p < 0.001) | Pearson r (log1p) = %.3f (p < 0.001)", 
                       spearman_raw$estimate, pearson_log$estimate),
    x = "Understory Openness Index (UOI)",
    y = "Standing Mammal Biomass Index (log1p scale)"
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    plot.subtitle = element_text(size = 8.5, hjust = 0.5, face = "italic"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Create figures and scratch output directories if they don't exist
dir.create("figures", showWarnings = FALSE)
dir.create("scratch", showWarnings = FALSE)

# Save plot
output_png <- "scratch/correlation_uoi_biomass.png"
ggsave(output_png, plot = p, width = 7, height = 5.5, dpi = 300)
cat(sprintf("✓ Saved scatter plot to %s\n", output_png))

# Save PDF version as well
output_pdf <- "scratch/correlation_uoi_biomass.pdf"
ggsave(output_pdf, plot = p, width = 7, height = 5.5)
cat(sprintf("✓ Saved scatter plot to %s\n", output_pdf))
