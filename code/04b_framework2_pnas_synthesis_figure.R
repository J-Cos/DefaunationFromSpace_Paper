# =============================================================================
# 04b_framework2_pnas_synthesis_figure.R
#
# Generates a unified, publication-quality PNAS multi-panel figure
# combining the key insights of Framework 2:
#   - Panel A: Square-aspect scatter plot comparing Standing Mammal Biomass (B_H)
#     against GEDI Understory Openness (UOI), showing prediction lines from the
#     univariate UOI-only model (LOBO) and the multivariate interaction model (AICc).
#   - Panel B: Model selection bar plot for standard AICc (from Figure 4C).
#   - Panel C: Model selection bar plot for generalizability LOBO (from Figure S2C).
# =============================================================================

library(ggplot2)
library(dplyr)
library(mgcv)
library(scales)
library(patchwork)

# --- 1. Load Data, Models, and Custom Theme ---
source("code/functions/calibration_helpers.R")
source("code/functions/theme_pnas.R")

# Load 5km scale calibrated cluster data
joined_data <- extract_scale_data(5000)

# Load selected best models
model_lobo <- readRDS("outputs/framework2_best_model.RDS")
model_aic <- readRDS("outputs/framework2_best_model_aic.RDS")

# Load covariate selection results CSV
df_selection <- read.csv("outputs/framework2_covariate_model_selection.csv")

# Ensure output directories exist
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("scratch", recursive = TRUE, showWarnings = FALSE)

# Posterior simulation helper to calculate robust confidence intervals from parameter covariance
get_sim_predictions <- function(model_obj, newdata, n_sim = 1000) {
  beta <- coef(model_obj)
  Vp <- vcov(model_obj)
  
  # Predict design matrix (lpmatrix)
  X <- predict(model_obj, newdata = newdata, type = "lpmatrix")
  
  # Draw from posterior parameter distribution using mgcv::rmvn
  set.seed(42)
  beta_sim <- mgcv::rmvn(n_sim, beta, Vp)
  
  # Project to response scale (Tweedie default link is log, so we use exp)
  link_preds <- X %*% t(beta_sim)
  resp_preds <- exp(link_preds)
  
  # Compute fit (median) and 95% confidence intervals (quantiles)
  newdata$fit <- apply(resp_preds, 1, median)
  newdata$lower <- apply(resp_preds, 1, quantile, probs = 0.025)
  newdata$upper <- apply(resp_preds, 1, quantile, probs = 0.975)
  return(newdata)
}

# --- 2. Predict Model Fits & 95% Confidence Intervals over UOI Range ---
uoi_seq <- seq(from = 0.918, to = 0.970, length.out = 300)
covs_to_fill <- c("elevation", "slope", "hnd", "precip", "clay", "forest_fraction")

# Predict LOBO (UOI-only) with posterior simulation CIs
pred_df_lobo <- data.frame(uoi = uoi_seq)
for (cv in covs_to_fill) pred_df_lobo[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
pred_df_lobo <- get_sim_predictions(model_lobo, pred_df_lobo)

# Predict AIC (UOI * Elephant Strict) - Absent
pred_df_aic_absent <- data.frame(uoi = uoi_seq)
pred_df_aic_absent$elephant_present_strict <- factor("Absent", levels = c("Absent", "Present"))
for (cv in covs_to_fill) pred_df_aic_absent[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
pred_df_aic_absent <- get_sim_predictions(model_aic, pred_df_aic_absent)

# Predict AIC (UOI * Elephant Strict) - Present
pred_df_aic_present <- data.frame(uoi = uoi_seq)
pred_df_aic_present$elephant_present_strict <- factor("Present", levels = c("Absent", "Present"))
for (cv in covs_to_fill) pred_df_aic_present[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
pred_df_aic_present <- get_sim_predictions(model_aic, pred_df_aic_present)

# Define custom prediction colors & linetypes
custom_colors <- c(
  "UOI Only (LOBO)" = "black",
  "Elephants Absent (AICc)" = "#D55E00",
  "Elephants Present (AICc)" = "#0072B2"
)

custom_linetypes <- c(
  "UOI Only (LOBO)" = "dashed",
  "Elephants Absent (AICc)" = "solid",
  "Elephants Present (AICc)" = "solid"
)

# --- 3. Panel A: Scatter Plot with Dual-Model Predictions & CIs ---
p_a <- ggplot() +
  # Points colored by elephant_present_strict
  geom_point(data = joined_data, 
             aes(x = uoi, y = B_H_index, fill = elephant_present_strict, size = trap_days, alpha = w_temp_cluster, shape = basin),
             color = "black", stroke = 0.3) +
  # Ribbons plotted behind lines/points for clean layering
  geom_ribbon(data = pred_df_lobo, aes(x = uoi, ymin = lower, ymax = upper), fill = "black", alpha = 0.08) +
  geom_ribbon(data = pred_df_aic_absent, aes(x = uoi, ymin = lower, ymax = upper), fill = "#D55E00", alpha = 0.08) +
  geom_ribbon(data = pred_df_aic_present, aes(x = uoi, ymin = lower, ymax = upper), fill = "#0072B2", alpha = 0.08) +
  
  geom_line(data = pred_df_lobo, aes(x = uoi, y = fit, color = "UOI Only (LOBO)", linetype = "UOI Only (LOBO)"), linewidth = 0.75) +
  geom_line(data = pred_df_aic_absent, aes(x = uoi, y = fit, color = "Elephants Absent (AICc)", linetype = "Elephants Absent (AICc)"), linewidth = 0.75) +
  geom_line(data = pred_df_aic_present, aes(x = uoi, y = fit, color = "Elephants Present (AICc)", linetype = "Elephants Present (AICc)"), linewidth = 0.75) +
  
  scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22), name = "Region:") +
  scale_size_continuous(name = "Effort (Trap-days):", range = c(1.2, 4.0), breaks = c(100, 1000, 5000, 15000)) +
  scale_alpha_continuous(name = "Temporal Alignment:", range = c(0.25, 1.0), breaks = c(0.1, 0.5, 1.0), labels = c("Hist.", "Interm.", "Contemp.")) +
  scale_fill_manual(values = pal_elephant_binary, name = "Elephant Status:") +
  scale_color_manual(values = custom_colors, name = "Model Predictions:") +
  scale_linetype_manual(values = custom_linetypes, name = "Model Predictions:") +
  scale_x_continuous(breaks = seq(0.92, 0.97, by = 0.01)) +
  scale_y_continuous(trans = "log1p", labels = comma_format(), breaks = c(0, 10, 100, 1000, 3000)) +
  coord_cartesian(xlim = c(0.918, 0.970), ylim = c(0, 5000)) +
  labs(
    title = "A",
    x = "GEDI Understory Openness Index (UOI)",
    y = "Total Mammal Biomass Index (log1p scale)"
  ) +
  guides(
    size = "none",
    alpha = "none",
    fill = "none",
    shape = guide_legend(override.aes = list(size = 2, fill = "gray50")),
    color = guide_legend(override.aes = list(linewidth = 0.8)),
    linetype = guide_legend()
  ) +
  theme_pnas(base_size = 7.5) +
  theme(
    legend.position = "right",
    aspect.ratio = 1.0,
    plot.title = element_text(face = "bold", size = 8.5, margin = margin(b = 6, t = 4)),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4)),
    plot.margin = margin(t = 1.5, r = 4, b = 1.5, l = 4, unit = "pt"),
    legend.title = element_text(size = 7.0, face = "bold"),
    legend.text = element_text(size = 6.5),
    legend.key.width = unit(0.4, "cm")
  )

# --- 4. Refit Candidate Models to obtain summaries for plotmath formatting ---
formulas_list <- list(
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

cat("Refitting candidate models to extract p-value significance indicators...\n")
full_models <- list()
for (m_name in names(formulas_list)) {
  full_models[[m_name]] <- gam(formulas_list[[m_name]], family = tw(),
                               weights = w_combined_norm, data = joined_data, method = "ML")
}

# Formatting function for y-axis labels
format_model_label_lobo <- function(model_name, model_obj) {
  clean_name <- sub("^M2\\.[0-9\\.]+[a-z_]*: ", "", model_name)
  tokens <- strsplit(clean_name, "\\s+")[[1]]
  tokens <- tokens[tokens != ""]
  p_table <- summary(model_obj)$p.table
  
  var_map <- list(
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

# --- 5. Panel B: Standard AICc Model Selection (Top 15) ---
plot_df_aic <- df_selection %>%
  arrange(Full_AICc) %>%
  head(15)

plot_df_aic$plotmath_label <- sapply(1:nrow(plot_df_aic), function(i) {
  format_model_label_lobo(plot_df_aic$Model[i], full_models[[plot_df_aic$Model[i]]])
})

plot_sel_df_aic <- plot_df_aic %>%
  mutate(
    CleanName = sub("^M2\\.[0-9\\.]+[a-z_]*: ", "", Model),
    CleanName = factor(CleanName, levels = rev(CleanName))
  )

ordered_exprs_aic <- plot_sel_df_aic$plotmath_label[match(levels(plot_sel_df_aic$CleanName), plot_sel_df_aic$CleanName)]
parsed_labels_aic <- parse(text = ordered_exprs_aic)

plot_sel_df_aic$Full_DevExpl_Pct <- plot_sel_df_aic$Full_DevExpl * 100

p_b <- plot_model_selection_bars(
  plot_df = plot_sel_df_aic,
  x_var = "Full_DevExpl_Pct",
  fill_var = "delta_AICc",
  fill_label = "Delta AICc",
  x_label = "Model Deviance Explained (%)",
  plot_title = "B",
  parsed_labels = parsed_labels_aic
) + labs(y = NULL) + theme(aspect.ratio = 0.25, plot.margin = margin(t = 1.5, r = 4, b = 1.5, l = 4, unit = "pt"))

# --- 6. Panel C: Generalizability Model Selection (Top 15) ---
plot_df_lobo <- df_selection %>%
  head(15)

plot_df_lobo$plotmath_label <- sapply(1:nrow(plot_df_lobo), function(i) {
  format_model_label_lobo(plot_df_lobo$Model[i], full_models[[plot_df_lobo$Model[i]]])
})

plot_sel_df_lobo <- plot_df_lobo %>%
  mutate(
    CleanName = sub("^M2\\.[0-9\\.]+[a-z_]*: ", "", Model),
    CleanName = factor(CleanName, levels = rev(CleanName))
  )

ordered_exprs_lobo <- plot_sel_df_lobo$plotmath_label[match(levels(plot_sel_df_lobo$CleanName), plot_sel_df_lobo$CleanName)]
parsed_labels_lobo <- parse(text = ordered_exprs_lobo)

plot_sel_df_lobo$Full_DevExpl_Pct <- plot_sel_df_lobo$Full_DevExpl * 100

p_c <- plot_model_selection_bars(
  plot_df = plot_sel_df_lobo,
  x_var = "Full_DevExpl_Pct",
  fill_var = "OOS_MAE_log",
  fill_label = "OOS MAE (log1p)",
  x_label = "Model Deviance Explained (%)",
  plot_title = "C",
  parsed_labels = parsed_labels_lobo
) + labs(y = NULL) + theme(aspect.ratio = 0.25, plot.margin = margin(t = 1.5, r = 4, b = 1.5, l = 4, unit = "pt"))

# --- 8. Assemble Multipanel Layout via Patchwork ---
# Stack Panels A, B, C vertically. Patchwork automatically aligns their plotting areas.
fig_final <- p_a + p_b + p_c + plot_layout(ncol = 1)

# --- 9. Save Figures ---
output_png <- "figures/framework2_pnas_synthesis_figure.png"
output_pdf <- "figures/framework2_pnas_synthesis_figure.pdf"

save_pnas(plot = fig_final, filename = output_png, type = "double", height_cm = 20.7)
save_pnas(plot = fig_final, filename = output_pdf, type = "double", height_cm = 20.7)

# Also save to scratch for backward compatibility
save_pnas(plot = fig_final, filename = "scratch/framework2_pnas_synthesis_figure.png", type = "double", height_cm = 20.7)
save_pnas(plot = fig_final, filename = "scratch/framework2_pnas_synthesis_figure.pdf", type = "double", height_cm = 20.7)

cat(sprintf("✓ Successfully saved unified synthesis figure to:\n  - %s\n  - %s\n  - scratch/framework2_pnas_synthesis_figure.png\n  - scratch/framework2_pnas_synthesis_figure.pdf\n", output_png, output_pdf))
