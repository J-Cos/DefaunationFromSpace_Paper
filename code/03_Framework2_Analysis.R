# =============================================================================
# code/03_Framework2_Analysis.R
#
# Performs covariate model selection for Framework 2 (Spaceborne Biomass Prediction)
# using Tweedie GLMs (via mgcv::gam) across 10 candidate models, evaluates them via AIC,
# and generates a publication-quality three-panel PNAS-style figure.
#
# All logic is encapsulated in a clean, unit-testable function.
# =============================================================================

library(terra)
library(dplyr)
library(readr)
library(mgcv)
library(ggplot2)
library(scales)
library(cowplot)

#' Run Framework 2 Analysis
#'
#' Fits a series of Tweedie GLMs predicting standing mammal biomass from GEDI understory
#' openness, conducts model selection via AIC, and outputs diagnostic plots.
#'
#' @param scale_m Numeric. Spatial resolution in meters (default: 5000)
#  @param outputs_dir Character. Directory to save model outputs (default: "outputs")
#' @param figures_dir Character. Directory to save figures (default: "figures")
#'
#' @return A list containing the model selection results data frame and the best fitted model object.
#' @export
run_framework2_analysis <- function(scale_m = 5000, outputs_dir = "outputs", figures_dir = "figures") {
  cat(sprintf("=== Starting Framework 2 Model Selection & Plotting (%d m scale) ===\n\n", scale_m))
  
  # --- 1. Ingest and Calibrate Scale-Specific Cluster Data ---------------------
  source("code/functions/calibration_helpers.R")
  joined_data <- extract_scale_data(scale_m)
  
  cat("✓ Merged and calibrated data successfully. N =", nrow(joined_data), "clusters.\n")
  
  # --- 2. Tweedie GLM Model Selection (10 Candidate Covariate Models) -----------
  cat("\nRunning Tweedie GLM Model Selection across 10 formulations...\n")
  
  models_list <- list(
    # --- Base Models ---
    "M2.1: UOI Only"                                        = gam(B_H_index ~ uoi, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.2: Basin Only"                                      = gam(B_H_index ~ basin, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    
    # --- UOI + Environmental Covariate Set ---
    "M2.3: UOI + Elev"                                      = gam(B_H_index ~ uoi + elevation, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.4: UOI + Slope"                                     = gam(B_H_index ~ uoi + slope, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.5: UOI + HAND"                                      = gam(B_H_index ~ uoi + hnd, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.6: UOI + Precip"                                    = gam(B_H_index ~ uoi + precip, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.7: UOI + Clay"                                      = gam(B_H_index ~ uoi + clay, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.8: UOI + Forest"                                    = gam(B_H_index ~ uoi + forest_fraction, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),

    # --- Main Effect Backbone Set ---
    "M2.9: UOI + Basin"                                     = gam(B_H_index ~ uoi + basin, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.10: UOI + Basin + Elev"                             = gam(B_H_index ~ uoi + basin + elevation, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.11: UOI + Basin + Slope"                            = gam(B_H_index ~ uoi + basin + slope, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.12: UOI + Basin + HAND"                             = gam(B_H_index ~ uoi + basin + hnd, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.13: UOI + Basin + Precip"                           = gam(B_H_index ~ uoi + basin + precip, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.14: UOI + Basin + Clay"                             = gam(B_H_index ~ uoi + basin + clay, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.15: UOI + Basin + Forest"                           = gam(B_H_index ~ uoi + basin + forest_fraction, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    
    # --- Interaction Effect Backbone Set (Decoupled / Compliance) ---
    "M2.16: UOI * Basin"                                    = gam(B_H_index ~ uoi * basin, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.17: UOI * Basin + Elev"                             = gam(B_H_index ~ uoi * basin + elevation, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.18: UOI * Basin + Slope"                            = gam(B_H_index ~ uoi * basin + slope, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.19: UOI * Basin + HAND"                             = gam(B_H_index ~ uoi * basin + hnd, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.20: UOI * Basin + Precip"                           = gam(B_H_index ~ uoi * basin + precip, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.21: UOI * Basin + Clay"                             = gam(B_H_index ~ uoi * basin + clay, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2.22: UOI * Basin + Forest"                           = gam(B_H_index ~ uoi * basin + forest_fraction, family = tw(), weights = w_combined_norm, data = joined_data, method = "REML")
  )
  
  # Compile results table
  results_df <- data.frame(
    Model = names(models_list),
    AIC = sapply(models_list, AIC),
    BIC = sapply(models_list, BIC),
    LogLik = sapply(models_list, function(m) as.numeric(logLik(m))),
    edf = sapply(models_list, function(m) sum(m$edf)),
    dev_expl = sapply(models_list, function(m) summary(m)$dev.expl),
    IsSignificant = sapply(models_list, function(m) {
      p_table <- summary(m)$p.table
      non_intercept_rows <- which(rownames(p_table) != "(Intercept)")
      if (length(non_intercept_rows) == 0) return(FALSE)
      p_vals <- p_table[non_intercept_rows, ncol(p_table), drop = TRUE]
      any(p_vals < 0.05)
    }),
    stringsAsFactors = FALSE
  )
  
  results_df <- results_df %>%
    mutate(delta_AIC = AIC - min(AIC)) %>%
    arrange(AIC)
  
  print(results_df)
  
  dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  write_csv(results_df, file.path(outputs_dir, "framework2_covariate_model_selection.csv"))
  
  # Identify the top best model on AIC
  best_model_name <- results_df$Model[1]
  best_model <- models_list[[best_model_name]]
  cat(sprintf("\n★ Selected Best-Fitting Model: %s (AIC: %.2f, d_AIC: 0.00)\n\n", best_model_name, AIC(best_model)))
  
  # Save best model RDS objects
  saveRDS(formula(best_model), file.path(outputs_dir, "framework2_best_formula.RDS"))
  saveRDS(best_model, file.path(outputs_dir, "framework2_best_model.RDS"))
  
  # --- 3. Generate Panels for Figure 3 -----------------------------------------
  cat("Generating PNAS-styled figure panels...\n")
  source("code/functions/theme_pnas.R")
  
  pal_basin <- c("Amazon" = "#E65100", "Congo" = "#1B5E20")
  
  uoi_seq <- seq(from = 0.918, to = 0.970, length.out = 300)
  
  # --- Panel A: Total standing biomass vs. GEDI UOI ---
  pred_df_amazon <- data.frame(uoi = uoi_seq, basin = factor("Amazon", levels = c("Amazon", "Congo")))
  pred_df_congo <- data.frame(uoi = uoi_seq, basin = factor("Congo", levels = c("Amazon", "Congo")))
  
  covs_to_fill <- c("elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  for (cv in covs_to_fill) {
    pred_df_amazon[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
    pred_df_congo[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
  }
  
  pred_df_amazon$fit <- predict(best_model, newdata = pred_df_amazon, type = "response")
  pred_df_congo$fit <- predict(best_model, newdata = pred_df_congo, type = "response")
  
  pred_total_plot <- rbind(pred_df_amazon, pred_df_congo)
  
  p_a <- ggplot() +
    geom_point(data = joined_data, aes(x = uoi, y = B_H_index, fill = w_temp_cluster, size = trap_days, alpha = homogeneity, shape = basin),
               color = "black", stroke = 0.3) +
    geom_line(data = pred_total_plot, aes(x = uoi, y = fit, color = basin), linewidth = 0.75) +
    
    scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21), name = "Basin/Continent") +
    scale_color_manual(values = pal_basin, name = "Basin/Continent") +
    scale_size_continuous(name = "Effort (Trap-days)", range = c(1.2, 4.0), breaks = c(100, 1000, 5000, 15000)) +
    scale_alpha_continuous(name = "Spatial Homogeneity", range = c(0.35, 1.0), breaks = c(0, 0.5, 1.0), labels = c("Low", "Medium", "High")) +
    scale_fill_gradientn(
      colors = c("#D32F2F", "#F57C00", "#FBC02D", "#388E3C"),
      values = c(0, 0.25, 0.5, 1.0),
      limits = c(0.1, 1.0),
      name = "Temporal Weight (W_temp)"
    ) +
    scale_x_continuous(breaks = seq(0.92, 0.97, by = 0.01), limits = c(0.918, 0.970)) +
    scale_y_continuous(trans = "log1p", labels = comma_format(), breaks = c(0, 10, 100, 1000, 3000), limits = c(0, 5000)) +
    labs(
      title = "A. Standing Mammal Biomass vs. GEDI Openness",
      subtitle = sprintf("Best Fit: %s", best_model_name),
      x = "GEDI Understory Openness Index (UOI)",
      y = "Total Mammal Biomass Index (log1p scale)"
    ) +
    theme_pnas(base_size = 7.5) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 8.5, margin = margin(b = 6, t = 4)),
      axis.title.x = element_text(margin = margin(t = 4)),
      axis.title.y = element_text(margin = margin(r = 4)),
      plot.margin = margin(t = 6, r = 4, b = 6, l = 4, unit = "pt")
    )
  
  # --- Panel B: Residual stability and temporal independence ---
  joined_data$residuals <- residuals(best_model, type = "deviance")
  
  p_b <- ggplot(joined_data, aes(x = w_temp_cluster, y = residuals)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#555555", linewidth = 0.4) +
    geom_point(aes(fill = w_temp_cluster, size = trap_days, alpha = homogeneity, shape = basin), color = "black", stroke = 0.3) +
    geom_smooth(method = "lm", formula = y ~ x, color = "#2E7D32", linewidth = 0.6, se = TRUE, alpha = 0.1) +
    
    scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21), name = "Basin/Continent") +
    scale_size_continuous(name = "Effort (Trap-days)", range = c(1.2, 4.0), breaks = c(100, 1000, 5000, 15000)) +
    scale_alpha_continuous(name = "Spatial Homogeneity", range = c(0.35, 1.0), breaks = c(0, 0.5, 1.0), labels = c("Low", "Medium", "High")) +
    scale_fill_gradientn(
      colors = c("#D32F2F", "#F57C00", "#FBC02D", "#388E3C"),
      values = c(0, 0.25, 0.5, 1.0),
      limits = c(0.1, 1.0),
      name = "Temporal Weight (W_temp)"
    ) +
    scale_x_continuous(breaks = seq(0.1, 1.0, by = 0.2), limits = c(0.08, 1.02)) +
    scale_y_continuous(breaks = seq(-10, 8, by = 2), limits = c(-10.2, 8.2)) +
    labs(
      title = "B. Residual Independence & Temporal Stability",
      subtitle = sprintf("Deviance Residuals from %s", best_model_name),
      x = "Cluster Temporal Alignment Weight (W_temp)",
      y = "Best Model Deviance Residuals"
    ) +
    theme_pnas(base_size = 7.5) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 8.5, margin = margin(b = 6, t = 4)),
      axis.title.x = element_text(margin = margin(t = 4)),
      axis.title.y = element_text(margin = margin(r = 4)),
      plot.margin = margin(t = 6, r = 4, b = 6, l = 4, unit = "pt")
    )
  
  # --- Panel C: Model Selection Bar Plot ---
  # Helper function for dynamic plotmath bolding of significant variables
  format_model_label_f2 <- function(model_name, model_obj) {
    clean_name <- sub("^M2\\.[0-9a-z]+: ", "", model_name)
    clean_name <- gsub(" (Shared)", "", clean_name, fixed = TRUE)
    
    tokens <- strsplit(clean_name, "\\s+")[[1]]
    tokens <- tokens[tokens != ""]
    
    p_table <- summary(model_obj)$p.table
    
    var_map <- list(
      "Biomass"       = "B_H_index",
      "Megafauna"     = "B_H_gt100",
      "Basin"         = "basinCongo",
      "UOI"           = "uoi",
      "Elevation"     = "elevation",
      "Slope"         = "slope",
      "HAND"          = "hnd",
      "Precipitation" = "precip",
      "Clay"          = "clay",
      "Forest"        = "forest_fraction",
      "UOI:Basin"     = "uoi:basinCongo",
      "Biomass:Basin" = "B_H_index:basinCongo"
    )
    
    plotmath_tokens <- sapply(tokens, function(tok) {
      if (tok %in% c("+", "*", ":")) return(sprintf("plain(\" %s \")", tok))
      if (tok == "Only") return(sprintf("plain(\" %s\")", tok))
      
      matched_term <- var_map[[tok]]
      if (!is.null(matched_term)) {
        is_sig <- FALSE
        if (matched_term %in% rownames(p_table)) {
          p_val <- p_table[matched_term, ncol(p_table)]
          is_sig <- !is.na(p_val) && p_val < 0.05
        }
        return(ifelse(is_sig, sprintf("bold(\"%s\")", tok), sprintf("plain(\"%s\")", tok)))
      } else {
        return(sprintf("plain(\"%s\")", tok))
      }
    })
    
    paste(plotmath_tokens, collapse = " * ")
  }
  
  results_df$plotmath_label <- sapply(1:nrow(results_df), function(i) {
    format_model_label_f2(results_df$Model[i], models_list[[results_df$Model[i]]])
  })
  
  plot_sel_df <- results_df %>%
    mutate(
      CleanName = sub("^M2\\.[0-9a-z]+: ", "", Model),
      CleanName = factor(CleanName, levels = rev(CleanName))
    )
  
  ordered_exprs <- plot_sel_df$plotmath_label[match(levels(plot_sel_df$CleanName), plot_sel_df$CleanName)]
  parsed_labels <- parse(text = ordered_exprs)
  
  p_c <- ggplot(plot_sel_df, aes(x = dev_expl * 100, y = CleanName, fill = delta_AIC)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.2) +
    scale_y_discrete(labels = parsed_labels) +
    scale_fill_gradientn(
      colors = c("#0D47A1", "#29B6F6", "#E0E0E0", "#FF7043", "#D84315"),
      name = "Delta AIC"
    ) +
    labs(
      title = "C. Spaceborne Model Selection (Tweedie GLMs)",
      x = "Model Deviance Explained (%)",
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
      legend.text = element_text(size = 5.5),
      legend.key.width = unit(0.15, "cm"),
      legend.key.height = unit(0.35, "cm"),
      legend.margin = margin(l = 2, r = 2, unit = "pt"),
      plot.margin = margin(t = 6, r = 4, b = 6, l = 4, unit = "pt")
    )
  
  # --- 4. Construct Clean Shared Legend ---
  p_legend_obj <- ggplot(joined_data) +
    geom_point(aes(x = uoi, y = B_H_index, fill = w_temp_cluster, size = trap_days, alpha = homogeneity, shape = basin), color = "black", stroke = 0.3) +
    scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21), name = "Basin:") +
    scale_size_continuous(name = "Effort (Trap-days):", breaks = c(100, 1000, 5000, 15000), range = c(1.2, 4.0)) +
    scale_alpha_continuous(name = "Spatial Homogeneity:", range = c(0.35, 1.0), breaks = c(0, 0.5, 1.0), labels = c("Low", "Medium", "High")) +
    scale_fill_gradientn(
      colors = c("#D32F2F", "#F57C00", "#FBC02D", "#388E3C"),
      values = c(0, 0.25, 0.5, 1.0),
      limits = c(0.1, 1.0),
      name = "Temporal Alignment Weight (W_temp):"
    ) +
    theme_pnas(base_size = 7.5) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 7.0, face = "bold"),
      legend.text = element_text(size = 6.5)
    )
  
  shared_legend <- cowplot::get_legend(p_legend_obj)
  
  # --- 5. Assemble and Save Multipanel Figure ---
  cat("Assembling panels into 3-panel PNAS-style layout...\n")
  
  row1 <- cowplot::plot_grid(
    p_a, p_b,
    ncol = 2,
    align = "h",
    axis = "tb",
    rel_widths = c(1.0, 1.0)
  )
  
  fig_final <- cowplot::plot_grid(
    row1,
    p_c,
    shared_legend,
    ncol = 1,
    rel_heights = c(1.0, 0.9, 0.12)
  )
  
  fig_png_path <- file.path(figures_dir, "figure4.png")
  save_pnas(
    plot = fig_final,
    filename = fig_png_path,
    type = "double",
    height_cm = 12.5
  )
  
  # Copy to brain artifact directory
  brain_artifact_dir <- "/home/j/.gemini/antigravity/brain/913e5cea-7c99-4b21-8124-ea8455da8457"
  if (file.exists(brain_artifact_dir)) {
    file.copy(fig_png_path,
              file.path(brain_artifact_dir, "figure4.png"),
              overwrite = TRUE)
    # Also copy pdf version
    fig_pdf_path <- sub("\\.png$", ".pdf", fig_png_path)
    if (file.exists(fig_pdf_path)) {
      file.copy(fig_pdf_path,
                file.path(brain_artifact_dir, "figure4.pdf"),
                overwrite = TRUE)
    }
    cat("✓ Copied integrated figure to brain artifacts folder.\n")
  }
  
  cat(sprintf("✓ Saved figure to %s\n", fig_png_path))
  cat("=== Framework 2 Integrated Analysis Completed Successfully ===\n")
  
  return(list(
    results_df = results_df,
    best_model = best_model
  ))
}

# --- Execute directly if called from terminal ---
run_framework2_analysis()
