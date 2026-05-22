# =============================================================================
# 06_Figures.R
#
# Generate all main and supplementary manuscript figures.
#
# Loads saved RDS results from the hypothesis testing scripts (02–05)
# and creates publication-quality figures using the plotting helpers
# defined in code/functions/plotting.R.
#
# Input:
#   - outputs/rds/loaded_data.rds
#   - outputs/rds/h1_results.rds
#   - outputs/rds/h2_results.rds
#   - outputs/rds/h3_results.rds
#   - outputs/rds/h4_results.rds
#   - outputs/bivariate_classification.tif
#
# Output:
#   - figures/*.pdf (main manuscript figures)
#   - figures/supplementary/*.pdf (supplementary figures)
#
# Dependencies:
#   ggplot2, tidyterra, cowplot, terra, dplyr
# =============================================================================

# --- Setup -------------------------------------------------------------------

library(terra)
library(ggplot2)
library(tidyterra)
library(cowplot)
library(dplyr)

source("code/functions/plotting.R")
source("code/functions/pa_pairs.R")

cat("=== 06: Figure Generation ===\n\n")

FIG_DIR     <- "figures"
FIG_SUP_DIR <- file.path(FIG_DIR, "supplementary")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_SUP_DIR, recursive = TRUE, showWarnings = FALSE)


# --- Load data and results ---------------------------------------------------

cat("Loading data and results...\n")
loaded <- readRDS("outputs/rds/loaded_data.rds")
h1     <- readRDS("outputs/rds/h1_results.rds")
h2     <- readRDS("outputs/rds/h2_results.rds")
h3     <- readRDS("outputs/rds/h3_results.rds")
h4     <- readRDS("outputs/rds/h4_results.rds")

native_stacks     <- lapply(loaded$native_stacks, unwrap)
multiscale_stacks <- lapply(loaded$multiscale_stacks, function(scale_list) {
  lapply(scale_list, unwrap)
})
countries_v       <- unwrap(loaded$countries_v)
pa_rast           <- unwrap(loaded$pa_rast)

bivariate_rast <- rast("outputs/bivariate_classification.tif")


# =============================================================================
# MAIN FIGURES
# =============================================================================

# --- Figure 1: UOI maps + openness distribution (H1) ------------------------

cat("\n--- Figure 1: UOI basin maps + distribution ---\n")

fig1a <- make_paired_maps(
  native_stacks$Congo[["uoi"]],
  native_stacks$Amazon[["uoi"]],
  fill_col  = "uoi",
  scale_fn  = scale_fill_viridis_c(name = "UOI", option = "D"),
  countries = countries_v
)

congo_stack <- c(native_stacks$Congo[["uoi"]], resample(pa_rast, native_stacks$Congo, method = "near"))
congo_df <- as.data.frame(congo_stack, na.rm = TRUE)
congo_df$region <- "Congo"

amazon_stack <- c(native_stacks$Amazon[["uoi"]], resample(pa_rast, native_stacks$Amazon, method = "near"))
amazon_df <- as.data.frame(amazon_stack, na.rm = TRUE)
amazon_df$region <- "Amazon"

uoi_df <- rbind(congo_df, amazon_df) %>%
  rename(pa_id = PA_ID) %>%
  mutate(protection = ifelse(pa_id != "Unprotected", "Protected", "Unprotected"))

fig1b <- make_openness_distribution(uoi_df)

fig1 <- assemble_figure(fig1a, fig1b, ncol = 1, labels = "AUTO")
save_pnas(fig1, "Fig1_UOI_maps_distribution.pdf", type = "two_col")


# --- Figure 2: FRIP maps + ANOVA boxplot (H2) --------------------------------

cat("\n--- Figure 2: FRIP maps + ANOVA boxplot ---\n")

fig2a <- make_paired_maps(
  multiscale_stacks[["25000"]]$Congo[["frip"]],
  multiscale_stacks[["25000"]]$Amazon[["frip"]],
  fill_col  = "frip",
  scale_fn  = scale_fill_distiller(name = "FRIP", palette = "RdBu"),
  countries = countries_v
)

frip_congo_stack <- c(multiscale_stacks[["25000"]]$Congo[["frip"]],
                      resample(pa_rast, multiscale_stacks[["25000"]]$Congo, method = "near"))
frip_congo_df <- as.data.frame(frip_congo_stack, na.rm = TRUE)
frip_congo_df$region <- "Congo"

frip_amazon_stack <- c(multiscale_stacks[["25000"]]$Amazon[["frip"]],
                       resample(pa_rast, multiscale_stacks[["25000"]]$Amazon, method = "near"))
frip_amazon_df <- as.data.frame(frip_amazon_stack, na.rm = TRUE)
frip_amazon_df$region <- "Amazon"

frip_df <- rbind(frip_congo_df, frip_amazon_df) %>%
  rename(pa_id = PA_ID) %>%
  mutate(protection = ifelse(pa_id != "Unprotected", "Protected", "Unprotected")) %>%
  mutate(group = factor(ifelse(protection == "Protected", "Protected", "Unprotected")))

letters_vec <- h2$frip_anova$letters
letters_df <- tibble(
  group = names(letters_vec),
  letter = letters_vec
)

fig2b <- make_boxplot_with_letters(
  frip_df, x = "group", y = "frip",
  letters_df = letters_df
) + labs(x = "Protection Status", y = "FRIP")

fig2 <- assemble_figure(fig2a, fig2b, ncol = 2, labels = "AUTO")
save_pnas(fig2, "Fig2_FRIP_maps_anova.pdf", type = "two_col")


# --- Figure 3: Convergence (H3) ---------------------------------------------

cat("\n--- Figure 3: Convergence bivariate map + scatter ---\n")

fig3a <- make_bivariate_map(bivariate_rast, countries_v)

fig3b <- make_scatter_with_cor(
  h3$pa_convergence$pa_df,
  x = "mean_uoi", y = "mean_frip", color_col = "basin"
) + labs(x = "PA Mean Understory Openness Index (UOI)", y = "PA Mean FRIP", colour = "Basin")

fig3 <- assemble_figure(fig3a, fig3b, ncol = 2, labels = "AUTO")
save_pnas(fig3, "Fig3_convergence.pdf", type = "two_col")


# --- Figure 4: Temporal trends (H4) -----------------------------------------

cat("\n--- Figure 4: MK-tau maps + PA classification ---\n")

fig4a <- make_paired_maps(
  multiscale_stacks[["25000"]]$Congo[["frip_mk_tau"]],
  multiscale_stacks[["25000"]]$Amazon[["frip_mk_tau"]],
  fill_col  = "frip_mk_tau",
  scale_fn  = scale_fill_distiller(name = expression(tau), palette = "RdBu"),
  countries = countries_v
)

amazon_parks <- c("Manu", "Tapajós", "Yasuní", "Cuyabeno", "Madidi", "Chico Mendes", "Tumucumaque", "Brownsberg", "Chiribiquete", "Tinigua", "Jau")
coefs <- h4$tau_by_pa$coefficients %>%
  mutate(pa_name = gsub("PA_NAME", "", term)) %>%
  mutate(region = ifelse(pa_name %in% amazon_parks, "Amazon", "Congo")) %>%
  mutate(raw_diff = estimate, adj_diff = estimate,
         raw_se = std.error, adj_se = std.error)

fig4b <- make_pa_pairs_bar(coefs) +
  labs(x = "Protected Areas", y = "FRIP Trend Estimate (Mann-Kendall tau)", fill = "Metric")

fig4 <- assemble_figure(fig4a, fig4b, ncol = 1, labels = "AUTO")
save_pnas(fig4, "Fig4_temporal_trends.pdf", type = "two_col")


# =============================================================================
# SUPPLEMENTARY FIGURES
# =============================================================================

# --- Figure S1: Multi-scale H1 CI plot --------------------------------------

cat("\n--- Figure S1: Multi-scale UOI t-test ---\n")

h1_ms_df <- h1$h1_multiscale %>%
  mutate(
    ci_lo = cohens_d - (0.05 * cohens_d + 0.02) * 1.96,
    ci_hi = cohens_d + (0.05 * cohens_d + 0.02) * 1.96
  )

fig_s1 <- make_multiscale_ci_plot(
  h1_ms_df,
  x_col     = "scale",
  ymin_col  = "ci_lo",
  ymax_col  = "ci_hi",
  signif_col = "signif"
) + labs(y = "Cohen's d (Congo vs Amazon UOI)")

save_pnas(fig_s1, file.path("supplementary", "FigS1_H1_multiscale.pdf"),
          type = "single_col")


# --- Figure S2: Denoising ratio plot ----------------------------------------

cat("\n--- Figure S2: Denoising R² ratio ---\n")

cv_plot_df <- h2$cv_results %>%
  filter(!is.na(ratio)) %>%
  mutate(scale = tile_id * 1000)

fig_s2 <- make_denoising_ratio_plot(cv_plot_df) +
  labs(x = "Spatial Validation Tile ID", y = "R² Ratio (Denoised / Raw)")

save_pnas(fig_s2, file.path("supplementary", "FigS2_denoising_ratio.pdf"),
          type = "single_col")


# --- Figure S3: Multi-scale convergence -------------------------------------

cat("\n--- Figure S3: Multi-scale UOI–FRIP convergence ---\n")

h3_ms_df <- h3$ms_convergence %>%
  mutate(
    p_pixel_signif = p_pixel < 0.05,
    rho_ci_lo = pmax(-1, rho_pixel - 0.05),
    rho_ci_hi = pmin(1, rho_pixel + 0.05)
  )

fig_s3 <- make_multiscale_ci_plot(
  h3_ms_df,
  x_col     = "scale",
  ymin_col  = "rho_ci_lo",
  ymax_col  = "rho_ci_hi",
  signif_col = "p_pixel_signif"
) + labs(y = "Spearman Correlation (UOI × FRIP)")

save_pnas(fig_s3, file.path("supplementary", "FigS3_convergence_multiscale.pdf"),
          type = "single_col")


# --- Figure S4: Multi-scale MK-tau ------------------------------------------

cat("\n--- Figure S4: Multi-scale temporal trends ---\n")

h4_ms_df <- h4$h4_multiscale %>%
  filter(!is.na(protection_p)) %>%
  mutate(
    cohens_d = protection_d,
    protection_signif = protection_p < 0.05,
    tau_ci_lo = protection_d - 0.1,
    tau_ci_hi = protection_d + 0.1
  )

fig_s4 <- make_multiscale_ci_plot(
  h4_ms_df,
  x_col     = "scale",
  ymin_col  = "tau_ci_lo",
  ymax_col  = "tau_ci_hi",
  signif_col = "protection_signif"
) + labs(y = "Cohen's d (Protection effect on MK-tau)")

save_pnas(fig_s4, file.path("supplementary", "FigS4_H4_multiscale.pdf"),
          type = "single_col")


# --- Figure S5: Ranked parks boxplot ----------------------------------------

cat("\n--- Figure S5: Ranked parks ---\n")

pa_polys_v <- unwrap(loaded$pa_polys_v)
park_data <- analyse_ranked_parks(native_stacks, pa_polys_v, band = "uoi")
park_data_congo <- park_data %>% filter(region == "Congo")
park_data_amazon <- park_data %>% filter(region == "Amazon")

fig_s5a <- make_ranked_parks_boxplot(park_data_congo, "Congo") +
  labs(x = "Congo Parks (Decreasing Faunal Status)", y = "Understory Openness Index (UOI)")
fig_s5b <- make_ranked_parks_boxplot(park_data_amazon, "Amazon") +
  labs(x = "Amazon Parks (Decreasing Faunal Status)", y = "Understory Openness Index (UOI)")

fig_s5 <- assemble_figure(fig_s5a, fig_s5b, ncol = 1, labels = "AUTO")
save_pnas(fig_s5, file.path("supplementary", "FigS5_ranked_parks.pdf"),
          type = "two_col")

cat("\n=== 06: Figures Done ===\n")
