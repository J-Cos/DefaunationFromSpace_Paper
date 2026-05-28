# =============================================================================
# code/03_Framework1_Analysis.R
#
# Performs covariate model selection for Framework 1 (GEDI UOI as Response) using
# Beta Regression (via mgcv::gam) across 24 candidate models (excluding elephant
# and megafauna history models), evaluates them via AIC, and generates a
# publication-quality three-panel PNAS-style figure.
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
#' @param outputs_dir Character. Directory to save model outputs (default: "outputs")
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
  cat("\nRunning Beta Regression Covariate Model Selection (Excluding Elephant & MegaHx)...\n")
  base_formulas <- list(
    # --- Base Models ---
    "M1: Biomass Only"                       = uoi ~ B_H_index,
    "M2: Megafauna Only"                     = uoi ~ B_H_gt100,
    
    # --- Biomass + Environmental Covariate Set ---
    "M4: Biomass + Elev"                     = uoi ~ B_H_index + elevation,
    "M5: Biomass + Slope"                    = uoi ~ B_H_index + slope,
    "M6: Biomass + HAND"                     = uoi ~ B_H_index + hnd,
    "M7: Biomass + Precip"                   = uoi ~ B_H_index + precip,
    "M8: Biomass + Clay"                     = uoi ~ B_H_index + clay,
    "M9: Biomass + Forest"                   = uoi ~ B_H_index + forest_fraction,

    # --- Basin Base and Main Effect Set ---
    "M24: Basin Only"                        = uoi ~ basin,
    "M25: Biomass + Basin"                   = uoi ~ B_H_index + basin,
    "M26: Biomass + Basin + Elev"            = uoi ~ B_H_index + basin + elevation,
    "M27: Biomass + Basin + Slope"           = uoi ~ B_H_index + basin + slope,
    "M28: Biomass + Basin + HAND"            = uoi ~ B_H_index + basin + hnd,
    "M29: Biomass + Basin + Precip"          = uoi ~ B_H_index + basin + precip,
    "M30: Biomass + Basin + Clay"            = uoi ~ B_H_index + basin + clay,
    "M31: Biomass + Basin + Forest"          = uoi ~ B_H_index + basin + forest_fraction,
    
    # --- Basin Interaction Effect Set ---
    "M32: Biomass * Basin"                   = uoi ~ B_H_index * basin,
    "M33: Biomass * Basin + Elev"            = uoi ~ B_H_index * basin + elevation,
    "M34: Biomass * Basin + Slope"           = uoi ~ B_H_index * basin + slope,
    "M35: Biomass * Basin + HAND"            = uoi ~ B_H_index * basin + hnd,
    "M36: Biomass * Basin + Precip"          = uoi ~ B_H_index * basin + precip,
    "M37: Biomass * Basin + Clay"            = uoi ~ B_H_index * basin + clay,
    "M38: Biomass * Basin + Forest"          = uoi ~ B_H_index * basin + forest_fraction
  )

  # Dynamically construct full set of base formulas, adding elephant alternates for any basin models
  expanded_base_formulas <- list()
  for (name in names(base_formulas)) {
    f <- base_formulas[[name]]
    expanded_base_formulas[[name]] <- f
    
    # If the formula contains 'basin', create the ElephantPossible and ElephantStrict alternates
    if ("basin" %in% all.vars(f)) {
      # 1. Elephant Possible
      name_possible <- gsub("Basin", "ElephantPossible", name)
      f_str_possible <- deparse(f)
      f_str_possible <- gsub("basin", "elephant_present_possible", f_str_possible)
      expanded_base_formulas[[name_possible]] <- as.formula(f_str_possible)
      
      # 2. Elephant Strict
      name_strict <- gsub("Basin", "ElephantStrict", name)
      f_str_strict <- deparse(f)
      f_str_strict <- gsub("basin", "elephant_present_strict", f_str_strict)
      expanded_base_formulas[[name_strict]] <- as.formula(f_str_strict)
    }
  }

  # Dynamically construct models_list with two alternatives added for each expanded base model
  models_list <- list()
  for (name in names(expanded_base_formulas)) {
    f <- expanded_base_formulas[[name]]
    
    # Base model
    models_list[[name]] <- gam(f, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML")
    
    # Alternative 1: current + biomass >100kg
    if (!("B_H_gt100" %in% all.vars(f))) {
      name_gt100 <- paste0(name, " + Megafauna")
      f_gt100 <- as.formula(paste(deparse(f), "+ B_H_gt100"))
      models_list[[name_gt100]] <- gam(f_gt100, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML")
    }
    
    # Alternative 2: current + biomass >1000kg
    if (!("B_H_gt1000" %in% all.vars(f))) {
      name_gt1000 <- paste0(name, " + Megafauna1000")
      f_gt1000 <- as.formula(paste(deparse(f), "+ B_H_gt1000"))
      models_list[[name_gt1000]] <- gam(f_gt1000, family = betar(link = "logit"), weights = w_combined_norm, data = joined_data, method = "REML")
    }
  }
  
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
  
  covs_to_fill <- c("B_H_gt100", "B_H_gt1000", "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
  
  # Detect whether best model uses megafaunaHistory, basin, or elephant as the grouping variable
  best_formula_vars <- all.vars(formula(best_model))
  uses_megahx <- "megafaunaHistory" %in% best_formula_vars
  uses_basin <- "basin" %in% best_formula_vars
  
  if (uses_megahx) {
    # --- Megafauna History-based fit lines (2 lines: NewWorld, OldWorld) ---
    pal_mega <- c("NewWorld" = "#E65100", "OldWorld" = "#00695C")
    pred_list <- lapply(levels(joined_data$megafaunaHistory), function(m) {
      nd <- data.frame(B_H_index = biomass_seq, megafaunaHistory = factor(m, levels = levels(joined_data$megafaunaHistory)))
      nd$basin <- factor(ifelse(m == "NewWorld", "Amazon", "Congo"), levels = levels(joined_data$basin))
      nd$elephant_present <- factor(ifelse(m == "NewWorld", "Absent", "Present"), levels = levels(joined_data$elephant_present))
      for (cv in covs_to_fill) nd[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
      nd$fit_link <- predict(best_model, newdata = nd, type = "link")
      nd$fit <- plogis(nd$fit_link)
      nd
    })
    pred_plot <- do.call(rbind, pred_list)
    
    p_a <- ggplot() +
      geom_point(data = joined_data, aes(x = B_H_index, y = uoi, fill = megafaunaHistory, size = trap_days, alpha = w_temp_cluster, shape = basin),
                 color = "black", stroke = 0.3) +
      geom_line(data = pred_plot, aes(x = B_H_index, y = fit, color = megafaunaHistory), linewidth = 0.75) +
      
      scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22), name = "Basin/Continent") +
      scale_color_manual(values = pal_mega, name = "Megafauna History") +
      scale_size_continuous(name = "Effort (Trap-days)", range = c(1.2, 4.0), breaks = c(100, 1000, 5000, 15000)) +
      scale_alpha_continuous(name = "Temporal Alignment Weight", range = c(0.25, 1.0), breaks = c(0.1, 0.5, 1.0), labels = c("Historical", "Intermediate", "Contemp.")) +
      scale_fill_manual(values = pal_mega, name = "Megafauna History") +
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
      
  } else if (uses_basin) {
    # --- Basin-based fit lines (3 lines: Amazon, Congo, SE_Asia) ---
    pal_basin <- c("Amazon" = "#E65100", "Congo" = "#1B5E20", "SE_Asia" = "#0D47A1")
    pred_list <- lapply(levels(joined_data$basin), function(b) {
      nd <- data.frame(B_H_index = biomass_seq, basin = factor(b, levels = levels(joined_data$basin)))
      nd$elephant_present <- factor(ifelse(b == "Amazon", "Absent", "Present"), levels = c("Absent", "Present"))
      nd$megafaunaHistory <- factor(ifelse(b %in% c("Congo", "SE_Asia"), "OldWorld", "NewWorld"), levels = levels(joined_data$megafaunaHistory))
      for (cv in covs_to_fill) nd[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
      nd$fit_link <- predict(best_model, newdata = nd, type = "link")
      nd$fit <- plogis(nd$fit_link)
      nd
    })
    pred_plot <- do.call(rbind, pred_list)
    
    p_a <- ggplot() +
      geom_point(data = joined_data, aes(x = B_H_index, y = uoi, fill = basin, size = trap_days, alpha = w_temp_cluster, shape = basin),
                 color = "black", stroke = 0.3) +
      geom_line(data = pred_plot, aes(x = B_H_index, y = fit, color = basin), linewidth = 0.75) +
      
      scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22), name = "Basin/Continent") +
      scale_color_manual(values = pal_basin, name = "Basin/Continent") +
      scale_size_continuous(name = "Effort (Trap-days)", range = c(1.2, 4.0), breaks = c(100, 1000, 5000, 15000)) +
      scale_alpha_continuous(name = "Temporal Alignment Weight", range = c(0.25, 1.0), breaks = c(0.1, 0.5, 1.0), labels = c("Historical", "Intermediate", "Contemp.")) +
      scale_fill_manual(values = pal_basin, name = "Basin/Continent") +
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
    
  } else {
    # Detect which elephant variable is in the formula
    uses_ele_strict <- "elephant_present_strict" %in% best_formula_vars
    uses_ele_possible <- "elephant_present_possible" %in% best_formula_vars
    
    ele_col <- if (uses_ele_strict) {
      "elephant_present_strict"
    } else if (uses_ele_possible) {
      "elephant_present_possible"
    } else {
      "elephant_present"
    }
    
    # --- Elephant-based fit lines (2 lines: Absent, Present) ---
    pred_df_absent <- data.frame(B_H_index = biomass_seq)
    pred_df_absent[[ele_col]] <- factor("Absent", levels = c("Absent", "Present"))
    pred_df_absent$basin <- factor("Amazon", levels = levels(joined_data$basin))
    pred_df_absent$megafaunaHistory <- factor("NewWorld", levels = levels(joined_data$megafaunaHistory))
    for (cv in covs_to_fill) {
      pred_df_absent[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
    }
    pred_df_absent$fit_link <- predict(best_model, newdata = pred_df_absent, type = "link")
    pred_df_absent$fit <- plogis(pred_df_absent$fit_link)
    
    pred_df_present <- data.frame(B_H_index = biomass_seq)
    pred_df_present[[ele_col]] <- factor("Present", levels = c("Absent", "Present"))
    pred_df_present$basin <- factor("Congo", levels = levels(joined_data$basin))
    pred_df_present$megafaunaHistory <- factor("OldWorld", levels = levels(joined_data$megafaunaHistory))
    for (cv in covs_to_fill) {
      pred_df_present[[cv]] <- median(joined_data[[cv]], na.rm = TRUE)
    }
    pred_df_present$fit_link <- predict(best_model, newdata = pred_df_present, type = "link")
    pred_df_present$fit <- plogis(pred_df_present$fit_link)
    
    pred_plot <- rbind(pred_df_absent, pred_df_present)
    
    p_a <- ggplot() +
      geom_point(data = joined_data, aes(x = B_H_index, y = uoi, fill = .data[[ele_col]], size = trap_days, alpha = w_temp_cluster, shape = basin),
                 color = "black", stroke = 0.3) +
      geom_line(data = pred_plot, aes(x = B_H_index, y = fit, color = .data[[ele_col]]), linewidth = 0.75) +
      
      scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22), name = "Basin/Continent") +
      scale_color_manual(values = c("Absent" = "#E06666", "Present" = "#2E7D32"), name = "Elephant Presence") +
      scale_size_continuous(name = "Effort (Trap-days)", range = c(1.2, 4.0), breaks = c(100, 1000, 5000, 15000)) +
      scale_alpha_continuous(name = "Temporal Alignment Weight", range = c(0.25, 1.0), breaks = c(0.1, 0.5, 1.0), labels = c("Historical", "Intermediate", "Contemp.")) +
      scale_fill_manual(values = c("Absent" = "#E06666", "Present" = "#2E7D32"), name = "Elephant Presence") +
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
  }
  
  # --- Panel B: Model Selection Bar Plot ---
  # Helper function for dynamic plotmath bolding of significant variables
  format_model_label <- function(model_name, model_obj) {
    clean_name <- gsub("^M[0-9\\.]+[a-z_]*: ", "", model_name)
    clean_name <- gsub(" (Shared)", "", clean_name, fixed = TRUE)
    
    tokens <- strsplit(clean_name, "\\s+")[[1]]
    tokens <- tokens[tokens != ""]
    
    p_table <- summary(model_obj)$p.table
    
    var_map <- list(
      "Biomass"          = "B_H_index",
      "Megafauna"        = "B_H_gt100",
      "Megafauna1000"    = "B_H_gt1000",
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
  
  results_df$plotmath_label <- sapply(1:nrow(results_df), function(i) {
    format_model_label(results_df$Model[i], models_list[[results_df$Model[i]]])
  })
  
  plot_sel_df <- results_df %>%
    head(20) %>%
    mutate(
      CleanName = gsub("^M[0-9\\.]+[a-z_]*: ", "", Model),
      CleanName = factor(CleanName, levels = rev(CleanName))
    )
  
  ordered_exprs <- plot_sel_df$plotmath_label[match(levels(plot_sel_df$CleanName), plot_sel_df$CleanName)]
  parsed_labels <- parse(text = ordered_exprs)
  
  p_b <- plot_model_selection_bars(
    plot_df = plot_sel_df,
    x_var = "R2",
    fill_var = "delta_AIC",
    fill_label = "Delta AIC",
    x_label = "Adjusted Pseudo-R² (Goodness of Fit)",
    plot_title = "C. Covariate Model Selection (Beta Regression)",
    parsed_labels = parsed_labels
  )
  
  # --- Panel C: Residual Diagnostic vs Temporal Offset ---
  joined_data$residuals <- residuals(best_model, type = "deviance")
  
  p_c <- ggplot(joined_data, aes(x = w_temp_cluster, y = residuals)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#555555", linewidth = 0.4) +
    geom_point(aes(fill = basin, size = trap_days, alpha = w_temp_cluster, shape = basin), color = "black", stroke = 0.3) +
    geom_smooth(method = "lm", aes(weight = w_combined_norm), formula = y ~ x, color = "#2E7D32", linewidth = 0.6, se = TRUE, alpha = 0.1) +
    
    scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22), name = "Basin/Continent") +
    scale_size_continuous(name = "Effort (Trap-days)", range = c(1.2, 4.0), breaks = c(100, 1000, 5000, 15000)) +
    scale_alpha_continuous(name = "Temporal Alignment Weight", range = c(0.25, 1.0), breaks = c(0.1, 0.5, 1.0), labels = c("Historical", "Intermediate", "Contemporaneous")) +
    scale_fill_manual(values = c("Amazon" = "#E65100", "Congo" = "#1B5E20", "SE_Asia" = "#0D47A1"), name = "Basin/Continent") +
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
    geom_point(aes(x = B_H_index, y = uoi, fill = basin, size = trap_days, alpha = w_temp_cluster, shape = basin), color = "black", stroke = 0.3) +
    scale_shape_manual(values = c("Amazon" = 24, "Congo" = 21, "SE_Asia" = 22), name = "Basin:") +
    scale_size_continuous(name = "Effort (Trap-days):", breaks = c(100, 1000, 5000, 15000), range = c(1.2, 4.0)) +
    scale_alpha_continuous(name = "Temporal Alignment:", range = c(0.25, 1.0), breaks = c(0.1, 0.5, 1.0), labels = c("Hist.", "Interm.", "Contemp.")) +
    scale_fill_manual(values = c("Amazon" = "#E65100", "Congo" = "#1B5E20", "SE_Asia" = "#0D47A1"), name = "Basin:") +
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
  brain_artifact_dir <- "/home/j/.gemini/antigravity/brain/8f51df52-4604-48e0-9ce8-1c52d1cb241c"
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
