# =============================================================================
# code/framework1_integrated_analysis.R
#
# Performs covariate model selection for Framework 1 (GEDI UOI as Response) using
# Beta Regression (via mgcv::gam) across 10 candidate models, evaluates them via AIC,
# and generates a publication-quality three-panel PNAS-style figure:
#   Panel A: Scatter plot of GEDI UOI vs. Biomass, colored by temporal offset weighting,
#            with best-fitting Beta Regression curves.
#   Panel B: Model selection comparison bar plot showing delta AIC.
#   Panel C: Model residuals vs. temporal offset weighting to validate stability.
# =============================================================================

library(terra)
library(dplyr)
library(readr)
library(mgcv)
library(ggplot2)
library(scales)
library(cowplot)

cat("=== Starting Integrated Framework 1 Model Selection & Plotting ===\n\n")

# --- 1. Load Protected Area / MCP Polygons and Raster Stacks -----------------
geojson_path <- "outputs/camera_traps_robust_buffered_mcps.geojson"
r_congo_path <- "outputs/EOdata/analysis_stack_5000_Congo.tif"
r_amazon_path <- "outputs/EOdata/analysis_stack_5000_Amazon.tif"

mcps <- terra::vect(geojson_path)
mcps_congo <- mcps[mcps$region == "Congo", ]
mcps_amazon <- mcps[mcps$region == "Amazon", ]

r_congo <- rast(r_congo_path)
r_amazon <- rast(r_amazon_path)

names(r_congo) <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                    "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
names(r_amazon) <- c("frip", "frip_mk_tau", "uoi", "uoi_sd", "rh98", "gedi_n",
                     "elevation", "slope", "hnd", "precip", "clay", "forest_fraction")

# --- 2. Extract Raster Pixel Values inside MCP Polygons ---------------------
extracted_congo <- terra::extract(r_congo, mcps_congo, df = TRUE)
mcp_congo_df <- as.data.frame(mcps_congo)
mcp_congo_df$ID <- 1:nrow(mcp_congo_df)
pixel_congo <- merge(extracted_congo, mcp_congo_df, by = "ID") %>%
  filter(!is.na(uoi) & !is.na(frip)) %>%
  select(-ID) %>%
  mutate(basin = "Congo")

extracted_amazon <- terra::extract(r_amazon, mcps_amazon, df = TRUE)
mcp_amazon_df <- as.data.frame(mcps_amazon)
mcp_amazon_df$ID <- 1:nrow(mcp_amazon_df)
pixel_amazon <- merge(extracted_amazon, mcp_amazon_df, by = "ID") %>%
  filter(!is.na(uoi) & !is.na(frip)) %>%
  select(-ID) %>%
  mutate(basin = "Amazon")

pixel_data <- rbind(pixel_congo, pixel_amazon)

# Aggregate to cluster level
joined_data <- pixel_data %>%
  group_by(cluster_id, region, basin, trap_days, n_species, B_H_index, M_H_index, B_H_gt50, B_H_gt100, megafauna_fraction) %>%
  summarise(
    n_pixels = n(),
    uoi_sd = ifelse(is.na(sd(uoi, na.rm = TRUE)), 0, sd(uoi, na.rm = TRUE)),
    uoi = mean(uoi, na.rm = TRUE),
    elevation = mean(elevation, na.rm = TRUE),
    slope = mean(slope, na.rm = TRUE),
    hnd = mean(hnd, na.rm = TRUE),
    precip = mean(precip, na.rm = TRUE),
    clay = mean(clay, na.rm = TRUE),
    forest_fraction = mean(forest_fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>% filter(trap_days >= 10)

# --- 3. Compute Temporal Weights per Cluster --------------------------------
cat("Loading camera trap detections database to compute temporal weights...\n")
det_all <- read_csv("outputs/camera_traps_joint_detections.csv", show_col_types = FALSE)

# Safe date parsing
det_all$start_date <- as.Date(det_all$start_date)
det_all$end_date <- as.Date(det_all$end_date)

# Haversine distance single-linkage clustering to map deployments to spatial clusters
haversine_dist <- function(lon1, lat1, lon2, lat2) {
  r <- 6371.0
  rad <- pi / 180
  dlon <- (lon2 - lon1) * rad
  dlat <- (lat2 - lat1) * rad
  lat1 <- lat1 * rad
  lat2 <- lat2 * rad
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * asin(sqrt(a))
  return(r * c)
}

coords_df <- det_all %>% 
  select(region, longitude, latitude) %>% 
  distinct() %>% 
  mutate(cluster_id_geo = "")

for (reg in unique(coords_df$region)) {
  sub_indices <- which(coords_df$region == reg)
  sub <- coords_df[sub_indices, ]
  n <- nrow(sub)
  if (n == 0) next
  if (n > 1) {
    dist_mat <- matrix(0, nrow=n, ncol=n)
    for (i in 1:n) {
      for (j in 1:n) {
        dist_mat[i,j] <- haversine_dist(sub$longitude[i], sub$latitude[i], sub$longitude[j], sub$latitude[j])
      }
    }
    hc <- hclust(as.dist(dist_mat), method="single")
    labels <- cutree(hc, h=11.1)
  } else {
    labels <- 1
  }
  coords_df$cluster_id_geo[sub_indices] <- paste0(reg, "_", sprintf("%02d", labels))
}

det_all <- det_all %>%
  left_join(coords_df, by = c("region", "longitude", "latitude"))

# Deployments dates and sampling effort
deployments <- det_all %>%
  select(region, cluster_id_geo, deployment_id, start_date, end_date, trap_days) %>%
  distinct()

# GEDI Launch baseline (April 17, 2019)
gedi_start <- as.Date("2019-04-17")

# Compute deployment temporal weights
deployments <- deployments %>%
  mutate(
    years_before_gedi = as.numeric(gedi_start - start_date) / 365.25,
    w_temp = case_when(
      years_before_gedi <= 1.0  ~ 1.0,
      years_before_gedi <= 6.0  ~ 0.5,
      years_before_gedi <= 11.0 ~ 0.25,
      TRUE                      ~ 0.1
    )
  )

# Aggregate to cluster level
cluster_temp_metrics <- deployments %>%
  group_by(cluster_id_geo) %>%
  summarise(
    w_temp_cluster = sum(trap_days * w_temp) / sum(trap_days),
    .groups = "drop"
  )

# Join temporal weights back to spatial dataset (matching names)
joined_data <- joined_data %>%
  left_join(cluster_temp_metrics, by = c("cluster_id" = "cluster_id_geo"))

# Fallback: if any NAs, set to median (or 1.0)
joined_data$w_temp_cluster[is.na(joined_data$w_temp_cluster)] <- median(joined_data$w_temp_cluster, na.rm = TRUE)

# Weights formulation
joined_data$uoi_se <- joined_data$uoi_sd / sqrt(joined_data$n_pixels)
reg_uoi <- median(joined_data$uoi_se[joined_data$uoi_se > 0])
joined_data$w_uoi <- log10(joined_data$trap_days) / (joined_data$uoi_se + reg_uoi)
joined_data$w_uoi_norm <- joined_data$w_uoi / mean(joined_data$w_uoi)
joined_data$basin <- factor(joined_data$basin, levels = c("Amazon", "Congo"))

# Spatial Homogeneity definition (inverse of standard error, normalized to 0-1)
raw_homo <- 1 / (joined_data$uoi_se + reg_uoi)
joined_data$homogeneity <- (raw_homo - min(raw_homo)) / (max(raw_homo) - min(raw_homo))

cat("✓ Merged data successfully. N =", nrow(joined_data), "clusters.\n")

# --- 4. Covariate Model Selection using Beta Regression -----------------------
cat("\nRunning Beta Regression Covariate Model Selection...\n")

models_list <- list(
  "M1: Biomass Only"           = gam(uoi ~ B_H_index, family = betar(link = "logit"), weights = w_uoi_norm, data = joined_data, method = "REML"),
  "M2: Megafauna Only"         = gam(uoi ~ B_H_gt100, family = betar(link = "logit"), weights = w_uoi_norm, data = joined_data, method = "REML"),
  "M3: Biomass + Basin"        = gam(uoi ~ B_H_index + basin, family = betar(link = "logit"), weights = w_uoi_norm, data = joined_data, method = "REML"),
  "M4: Biomass * Basin"        = gam(uoi ~ B_H_index * basin, family = betar(link = "logit"), weights = w_uoi_norm, data = joined_data, method = "REML"),
  "M5: Biomass + Basin + Elev"  = gam(uoi ~ B_H_index + basin + elevation, family = betar(link = "logit"), weights = w_uoi_norm, data = joined_data, method = "REML"),
  "M6: Biomass + Basin + Slope" = gam(uoi ~ B_H_index + basin + slope, family = betar(link = "logit"), weights = w_uoi_norm, data = joined_data, method = "REML"),
  "M7: Biomass + Basin + HAND"  = gam(uoi ~ B_H_index + basin + hnd, family = betar(link = "logit"), weights = w_uoi_norm, data = joined_data, method = "REML"),
  "M8: Biomass + Basin + Precip" = gam(uoi ~ B_H_index + basin + precip, family = betar(link = "logit"), weights = w_uoi_norm, data = joined_data, method = "REML"),
  "M9: Biomass + Basin + Clay"  = gam(uoi ~ B_H_index + basin + clay, family = betar(link = "logit"), weights = w_uoi_norm, data = joined_data, method = "REML"),
  "M10: Biomass + Basin + Forest" = gam(uoi ~ B_H_index + basin + forest_fraction, family = betar(link = "logit"), weights = w_uoi_norm, data = joined_data, method = "REML")
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
    p_vals <- p_table[non_intercept_rows, "Pr(>|z|)", drop = TRUE]
    all(p_vals < 0.05)
  }),
  stringsAsFactors = FALSE
)


results_df <- results_df %>%
  mutate(delta_AIC = AIC - min(AIC)) %>%
  arrange(AIC)

print(results_df)

# Save the model selection table to a text summary
dir.create("outputs", recursive = TRUE, showWarnings = FALSE)
write_csv(results_df, "outputs/framework1_covariate_model_selection.csv")

# Identify the best model
best_model_name <- results_df$Model[1]
best_model <- models_list[[best_model_name]]
cat(sprintf("\n★ Selected Best-Fitting Model: %s (AIC: %.2f, d_AIC: 0.00)\n\n", best_model_name, AIC(best_model)))

# --- 5. Generate Panels for PNAS Double-Column Figure -------------------------
cat("Generating PNAS figure panels...\n")
source("code/functions/theme_pnas.R")

# --- Panel A: Scatter Plot with Best Model Fit ---
biomass_seq <- seq(0, 5000, length.out = 300)
pred_df_amazon <- data.frame(B_H_index = biomass_seq, basin = factor("Amazon", levels = c("Amazon", "Congo")))
pred_df_congo <- data.frame(B_H_index = biomass_seq, basin = factor("Congo", levels = c("Amazon", "Congo")))

# Populate environmental covariates with median values if present in best model formula
covs_to_fill <- c("elevation", "slope", "hnd", "precip", "clay", "forest_fraction")
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
  # Empirical points
  geom_point(data = joined_data, aes(x = B_H_index, y = uoi, fill = w_temp_cluster, size = trap_days, alpha = homogeneity, shape = basin),
             color = "black", stroke = 0.3) +
  # Best model curves
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
plot_sel_df <- results_df %>%
  mutate(
    CleanName = gsub("M[0-9]+: ", "", Model),
    CleanName = factor(CleanName, levels = rev(CleanName))
  )

# Programmatically build plotmath labels (bold for fully significant models, plain otherwise)
y_levels <- levels(plot_sel_df$CleanName)
matched_sig <- plot_sel_df$IsSignificant[match(y_levels, plot_sel_df$CleanName)]
math_labels <- ifelse(matched_sig,
                      paste0("bold(\"", y_levels, "\")"),
                      paste0("plain(\"", y_levels, "\")"))
parsed_labels <- parse(text = math_labels)

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

# --- 6. Construct Clean Shared Legend ---
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

# --- 7. Assemble and Save Multipanel Figure ---
cat("Assembling panels into two-row PNAS-style layout...\n")

# Row 1: A and B side-by-side
row1 <- cowplot::plot_grid(
  p_a, p_c,
  ncol = 2,
  align = "h",
  axis = "tb",
  rel_widths = c(1.0, 1.0)
)

# Row 2: C (which is p_b) at full width
row2 <- p_b

fig_final <- cowplot::plot_grid(
  row1,
  row2,
  shared_legend,
  ncol = 1,
  rel_heights = c(1.0, 0.9, 0.12)
)

# Save as PNAS double-column figure (17.8 cm wide) with a generous height for 2 rows
save_pnas(
  plot = fig_final,
  filename = "outputs/framework1_integrated_pnas_figure.png",
  type = "double",
  height_cm = 12.5
)

# Copy to brain artifact directory
brain_artifact_dir <- "/home/j/.gemini/antigravity/brain/913e5cea-7c99-4b21-8124-ea8455da8457"
if (file.exists(brain_artifact_dir)) {
  file.copy("outputs/framework1_integrated_pnas_figure.png",
            file.path(brain_artifact_dir, "framework1_integrated_pnas_figure.png"),
            overwrite = TRUE)
  cat("✓ Copied integrated figure to brain artifacts folder.\n")
}

cat("✓ Saved figure to outputs/framework1_integrated_pnas_figure.png\n")
cat("=== Integrated Framework 1 Analysis Completed Successfully ===\n")
