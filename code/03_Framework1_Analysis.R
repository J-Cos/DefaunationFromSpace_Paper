# =============================================================================
# code/03_Framework1_Analysis.R
#
# Performs covariate model selection for Framework 1 (GEDI UOI as Response) using
# Beta Regression (via mgcv::gam) across 10 candidate models, evaluates them via AIC,
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

#' Run Framework 1 Analysis
#'
#' Runs GEDI understory openness (UOI) beta regressions, performs model selection
#' on AIC, generates diagnostic panels, and saves the output figure and model data.
#'
#' @param scale_m Numeric. Spatial resolution in meters (default: 5000)
#  @param outputs_dir Character. Directory to save model outputs (default: "outputs")
#' @param figures_dir Character. Directory to save figures (default: "figures")
#'
#' @return A list containing the model selection results data frame and the best fitted model object.
#' @export
run_framework1_analysis <- function(scale_m = 5000, outputs_dir = "outputs", figures_dir = "figures") {
  cat(sprintf("=== Starting Framework 1 Model Selection & Plotting (%d m scale) ===\n\n", scale_m))
  
  # --- 1. Ingest and Calibrate Scale-Specific Cluster Data ---------------------
  source("code/functions/calibration_helpers.R")
  joined_data <- extract_scale_data(scale_m)
  
  cat("✓ Merged and calibrated data successfully. N =", nrow(joined_data), "clusters.\n")
  
  # --- 2. Covariate Model Selection using Beta Regression -----------------------
  cat("\nRunning Beta Regression Covariate Model Selection...\n")
  
  models_list <- list(
    # --- Base Models ---
    "M1: Biomass Only"              = gam(uoi ~ B_H_index, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M2: Megafauna Only"            = gam(uoi ~ B_H_gt100, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M3: Basin Only"                = gam(uoi ~ basin, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    
    # --- Biomass + Environmental Covariate Set ---
    "M4: Biomass + Elev"            = gam(uoi ~ B_H_index + elevation, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M5: Biomass + Slope"           = gam(uoi ~ B_H_index + slope, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M6: Biomass + HAND"            = gam(uoi ~ B_H_index + hnd, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M7: Biomass + Precip"          = gam(uoi ~ B_H_index + precip, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M8: Biomass + Clay"            = gam(uoi ~ B_H_index + clay, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M9: Biomass + Forest"          = gam(uoi ~ B_H_index + forest_fraction, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),

    # --- Main Effect Backbone Set ---
    "M10: Biomass + Basin"           = gam(uoi ~ B_H_index + basin, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M11: Biomass + Basin + Elev"    = gam(uoi ~ B_H_index + basin + elevation, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M12: Biomass + Basin + Slope"   = gam(uoi ~ B_H_index + basin + slope, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M13: Biomass + Basin + HAND"    = gam(uoi ~ B_H_index + basin + hnd, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M14: Biomass + Basin + Precip"  = gam(uoi ~ B_H_index + basin + precip, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M15: Biomass + Basin + Clay"    = gam(uoi ~ B_H_index + basin + clay, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M16: Biomass + Basin + Forest"  = gam(uoi ~ B_H_index + basin + forest_fraction, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    
    # --- Interaction Effect Backbone Set (Decoupled) ---
    "M17: Biomass * Basin"           = gam(uoi ~ B_H_index * basin, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M18: Biomass * Basin + Elev"    = gam(uoi ~ B_H_index * basin + elevation, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M19: Biomass * Basin + Slope"   = gam(uoi ~ B_H_index * basin + slope, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M20: Biomass * Basin + HAND"    = gam(uoi ~ B_H_index * basin + hnd, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M21: Biomass * Basin + Precip"  = gam(uoi ~ B_H_index * basin + precip, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M22: Biomass * Basin + Clay"    = gam(uoi ~ B_H_index * basin + clay, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML"),
    "M23: Biomass * Basin + Forest"  = gam(uoi ~ B_H_index * basin + forest_fraction, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML")
  )
  
  # Compile results table
  results_df <- data.frame(
    Model = names(models_list),
    AIC = sapply(models_list, AIC),
    LogLik = sapply(models_list, function(m) as.numeric(logLik(m))),
    edf = sapply(models_list, function(m) sum(m$edf)),
    R2 = sapply(models_list, function(m) {
      r2 <- summary(m)$r.sq
      if (is.null(r2) || is.na(r2)) return(summary(m)$dev.expl)
      return(r2)
    }),
    dev_expl = sapply(models_list, function(m) summary(m)$dev.expl),
    IsSignificant = sapply(models_list, function(m) {
      p_table <- summary(m)$p.table
      non_intercept_rows <- which(rownames(p_table) != "(Intercept)")
      if (length(non_intercept_rows) == 0) return(FALSE)
      p_vals <- p_table[non_intercept_rows, ncol(p_table), drop = TRUE]
      all(p_vals < 0.05)
    }),
    stringsAsFactors = FALSE
  )
  
  results_df <- results_df %>%
    mutate(delta_AIC = AIC - min(AIC)) %>%
    arrange(AIC)
  
  print(results_df)
  
  dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  write_csv(results_df, file.path(outputs_dir, "framework1_covariate_model_selection.csv"))
  
  # Identify the best model
  best_model_name <- results_df$Model[1]
  best_model <- models_list[[best_model_name]]
  cat(sprintf("\n★ Selected Best-Fitting Model: %s (AIC: %.2f, d_AIC: 0.00)\n\n", best_model_name, AIC(best_model)))
  
  # Save best model RDS objects to outputs
  saveRDS(formula(best_model), file.path(outputs_dir, "framework1_best_formula.RDS"))
  saveRDS(best_model, file.path(outputs_dir, "framework1_best_model.RDS"))
  
  # --- 3. Generate Panels for PNAS Double-Column Figure -------------------------
  cat("Generating PNAS figure panels...\n")
  source("code/functions/theme_pnas.R")
  
  # --- Panel A: Scatter Plot with Best Model Fit ---
  biomass_seq <- seq(0, 5000, length.out = 300)
  pred_df_amazon <- data.frame(B_H_index = biomass_seq, basin = factor("Amazon", levels = c("Amazon", "Congo")))
  pred_df_congo <- data.frame(B_H_index = biomass_seq, basin = factor("Congo", levels = c("Amazon", "Congo")))
  
  covs_to_fill <- c("B_H_gt100", "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  for (cv in covs_to_fill) {
    pred_df_amazon[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
    pred_df_congo[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
  }
  
  pred_df_amazon$fit_link <- predict(best_model, newdata = pred_df_amazon, type = "link")
  pred_df_amazon$fit <- exp(pred_df_amazon$fit_link) / (1 + exp(pred_df_amazon$fit_link))
  
  pred_df_congo$fit_link <- predict(best_model, newdata = pred_df_congo, type = "link")
  pred_df_congo$fit <- exp(pred_df_congo$fit_link) / (1 + exp(pred_df_congo$fit_link))
  
  pred_plot <- rbind(pred_df_amazon, pred_df_congo)
  
  p_a <- ggplot() +
    geom_point(data = joined_data, aes(x = B_H_index, y = uoi, fill = w_temp_cluster, size = trap_days, alpha = homogeneity, shape = basin),
               color = "black", stroke = 0.3) +
    geom_line(data = pred_plot, aes(x = B_H_index, y = fit, color = basin), linewidth = 0.75) +
    
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
    scale_x_continuous(trans = "log1p", labels = comma_format(), breaks = c(0, 10, 100, 1000, 3000), limits = c(0, 5000)) +
    scale_y_continuous(breaks = seq(0.92, 0.97, by = 0.01), limits = c(0.918, 0.972)) +
    labs(
      title = "A. GEDI Openness vs. Standing Biomass",
      x = "Mammal Standing Biomass Index (log1p scale)",
      y = "GEDI Understory Openness Index (UOI)"
    ) +
    theme_pnas(base_size = 7.5) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 8.5, margin = margin(b = 6, t = 4)),
      axis.title.x = element_text(margin = margin(t = 4)),
      axis.title.y = element_text(margin = margin(r = 4)),
      plot.margin = margin(t = 6, r = 4, b = 6, l = 4, unit = "pt")
    )
  
  # --- Panel B: Model Selection Bar Plot ---
  # Helper function for dynamic plotmath bolding of significant variables
  format_model_label <- function(model_name, model_obj) {
    clean_name <- gsub("^M[0-9\\.]+[a-z]*: ", "", model_name)
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
      "Elev"          = "elevation",
      "Slope"         = "slope",
      "HAND"          = "hnd",
      "Precipitation" = "precip",
      "Precip"        = "precip",
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
    format_model_label(results_df$Model[i], models_list[[results_df$Model[i]]])
  })
  
  plot_sel_df <- results_df %>%
    mutate(
      CleanName = gsub("M[0-9]+: ", "", Model),
      CleanName = factor(CleanName, levels = rev(CleanName))
    )
  
  ordered_exprs <- plot_sel_df$plotmath_label[match(levels(plot_sel_df$CleanName), plot_sel_df$CleanName)]
  parsed_labels <- parse(text = ordered_exprs)
  
  p_b <- ggplot(plot_sel_df, aes(x = R2, y = CleanName, fill = delta_AIC)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.2) +
    scale_y_discrete(labels = parsed_labels) +
    scale_fill_gradientn(
      colors = c("#0D47A1", "#29B6F6", "#E0E0E0", "#FF7043", "#D84315"),
      name = "Delta AIC"
    ) +
    labs(
      title = "C. Covariate Model Selection (Beta Regression)",
      x = "Adjusted Pseudo-R² (Goodness of Fit)",
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
  
  # --- Panel C: Residual Diagnostic vs Temporal Offset ---
  joined_data$residuals <- residuals(best_model, type = "deviance")
  
  p_c <- ggplot(joined_data, aes(x = w_temp_cluster, y = residuals)) +
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
    scale_y_continuous(breaks = seq(-3, 3, by = 1), limits = c(-2.8, 2.8)) +
    labs(
      title = "B. Residual Independence & Temporal Stability",
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
  
  # --- 4. Construct Clean Shared Legend ---
  p_legend_obj <- ggplot(joined_data) +
    geom_point(aes(x = B_H_index, y = uoi, fill = w_temp_cluster, size = trap_days, alpha = homogeneity, shape = basin), color = "black", stroke = 0.3) +
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
  cat("Assembling panels into two-row PNAS-style layout...\n")
  
  row1 <- cowplot::plot_grid(
    p_a, p_c,
    ncol = 2,
    align = "h",
    axis = "tb",
    rel_widths = c(1.0, 1.0)
  )
  
  fig_final <- cowplot::plot_grid(
    row1,
    p_b,
    shared_legend,
    ncol = 1,
    rel_heights = c(1.0, 0.9, 0.12)
  )
  
  # Save to the designated figures directory
  fig_png_path <- file.path(figures_dir, "figure3.png")
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
              file.path(brain_artifact_dir, "figure3.png"),
              overwrite = TRUE)
    # Also copy pdf version
    fig_pdf_path <- sub("\\.png$", ".pdf", fig_png_path)
    if (file.exists(fig_pdf_path)) {
      file.copy(fig_pdf_path,
                file.path(brain_artifact_dir, "figure3.pdf"),
                overwrite = TRUE)
    }
    cat("✓ Copied integrated figure to brain artifacts folder.\n")
  }
  
  cat(sprintf("✓ Saved figure to %s\n", fig_png_path))
  cat("=== Framework 1 Integrated Analysis Completed Successfully ===\n")
  
  return(list(
    results_df = results_df,
    best_model = best_model
  ))
}

# --- Execute directly if called from terminal ---
run_framework1_analysis()
