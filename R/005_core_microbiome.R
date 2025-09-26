#####################################################################
# Corn Interrow/Row Analysis
#
# This script analyzes the core microbiome across different corn plots
# and locations, identifying ASVs that are consistently present across
# different combinations of samples.
#
# Author: Jaejin Lee
# Modified by: Bolívar Aponte Rolón
# Date: 2025-05-05
#####################################################################

# Load required libraries
source("R/utils/000_setup.R")


#--------------------------------------------------------
# Subset data by corn plant, plots, and locations
#--------------------------------------------------------

# Subset to include only Corn samples
corn_physeq <- subset_samples(sabr_2023_physeq, plant == "Corn")

corn_physeq_list <- get_prevalent_rare(
  corn_physeq,
  thresholds = c(90, 80, 70, 60, 30),
  detection = 0 / 100,
  include.lowest = FALSE
)

corn_prevalence_list <- list(
  full_community = sabr_2023_physeq,
  prevalent_90 = corn_physeq_list$prevalent_90,
  prevalent_80 = corn_physeq_list$prevalent_80,
  prevalent_70 = corn_physeq_list$prevalent_70,
  prevalent_60 = corn_physeq_list$prevalent_60,
  prevalent_30 = corn_physeq_list$prevalent_30
)

# ------------------------------------------------------------
# Plot relative abundance of ASVs vs gnha across prevalence levels
# ------------------------------------------------------------

# Build combined long data frame
corn_combined_rel_abund <- imap_dfr(
  corn_prevalence_list,
  ~ relab_prevalence_long(.x, prevalence_label = .y)
)

# Select top ASVs per location within each prevalence threshold
# Top 20 per location within each prevalence
per_loc_top <- corn_combined_rel_abund %>%
  group_by(prevalence_level, sampling_location, asv) %>%
  summarise(
    mean_rel = mean(rel_abundance, na.rm = TRUE),
    .groups = "drop_last"
  ) %>%
  slice_max(order_by = mean_rel, n = 20, with_ties = FALSE) %>%
  ungroup()

# Collapse to per-prevalence unique ASVs, rank by overall mean and cap at 20
corn_top20_pairs_union <- per_loc_top %>%
  group_by(prevalence_level, asv) %>%
  summarise(
    mean_rel_overall = mean(mean_rel, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(prevalence_level) %>%
  slice_max(order_by = mean_rel_overall, n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  select(prevalence_level, asv)

corn_combined_top_asv <- corn_combined_rel_abund %>%
  dplyr::inner_join(
    corn_top20_pairs_union,
    by = c("prevalence_level", "asv")
  ) %>%
  mutate(
    prevalence_level = factor(
      prevalence_level,
      levels = c(
        "full_community",
        "prevalent_30",
        "prevalent_60",
        "prevalent_70",
        "prevalent_80",
        "prevalent_90"
      )
    )
  )

#---------------------
# Scatter plot
#---------------------
# Scatter + smooth (lm) faceted by prevalence level (rows) and ASV (cols)
# p_scatter <- ggplot(
#   corn_combined_top_asv %>%
#     filter(!is.na(gnha)),
#   aes(x = gnha, y = rel_abundance, color = sampling_location)
# ) +
#   geom_point(alpha = 0.35, size = 0.9, ) +
#   geom_smooth(
#     aes(color = sampling_location),
#     method = "lm",
#     formula = y ~ x,
#     se = FALSE,
#     #color = "#d7301f",
#     linewidth = 0.5
#   ) +
#   ggpmisc::stat_poly_eq(
#     ggpmisc::use_label(c(
#       "adj.R2",
#       # "f",
#       "p"
#       # "n"
#     )),
#     #label.y = "top",
#     label.x = 0.02,
#     label.y = seq(0.95, by = -0.1, length.out = 2),
#     size = 1
#   ) +
#   facet_grid(
#     prevalence_level ~ fct_reorder(asv, parse_number(asv), .fun = min),
#     scales = "free_y"
#   ) +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
#   scale_x_log10() +
#   labs(
#     title = "CORN: ASV Relative Abundance vs gnha Across Prevalence Thresholds",
#     x = "log10(gnha)",
#     y = "Relative Abundance",
#     subtitle = sprintf(
#       "Top %d ASVs per location within each prevalence threshold selected by mean relative abundance",
#       20
#     )
#   ) +
#   theme_bw(base_size = 10) +
#   theme(
#     strip.text = element_text(size = 4),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.grid.minor = element_blank()
#   )

# p_scatter

# Make one plot per prevalence with facet_grid/wrap over ASV
nested_names <- corn_combined_top_asv %>%
  group_by(prevalence_level) %>%
  arrange(prevalence_level) %>%
  tidyr::nest()

scatter_plots_by_prev <- corn_combined_top_asv %>%
  group_by(prevalence_level) %>%
  arrange(prevalence_level) %>%
  group_split(.keep = TRUE) %>%
  setNames(., nested_names$prevalence_level) %>%
  imap(function(df_prev, prev_lbl) {
    ggplot(
      df_prev,
      aes(x = gnha, y = rel_abundance, color = sampling_location)
    ) +
      geom_point(alpha = 0.35, size = 0.9) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
      ggpmisc::stat_poly_eq(
        ggpmisc::use_label(c("p")),
        label.x = 0.02,
        label.y = c(0.95, 0.85),
        size = 2
      ) +
      facet_grid(
        prevalence_level ~ fct_reorder(asv, parse_number(asv), .fun = min),
        scales = "free_y"
      ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      scale_x_log10(
        breaks = c(0.1, 1, 10, 100, 1000),
        limits = c(0.1, 1000),
        labels = scales::label_number()
      ) +
      labs(
        title = paste("CORN: ASV Relative Abundance:", prev_lbl),
        x = "log10(gnha)",
        y = "Relative Abundance",
        subtitle = sprintf(
          "Top %d ASVs per location within each prevalence threshold selected by mean relative abundance",
          20
        )
      ) +
      theme_bw(base_size = 10) +
      theme(
        strip.text = element_text(size = 5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        panel.grid.minor = element_blank()
      )
  })

# final_fig <- patchwork::wrap_plots(scatter_plots_by_prev, ncol = 1)
# final_fig
scatter_plots_by_prev

safe_name <- function(x) str_replace_all(x, "[^A-Za-z0-9._-]+", "_")
purrr::iwalk(
  scatter_plots_by_prev,
  ~ {
    base <- file.path("data/output/plots/", paste0("scatter_", safe_name(.y)))
    ggsave(
      paste0(base, ".png"),
      .x,
      width = 14,
      height = 8,
      dpi = 600,
      units = "in"
    )
    ggsave(paste0(base, ".pdf"), .x, width = 14, height = 8, units = "in")
  }
)

#------------------
# Binned plots
#------------------

# Aggregate view: Mean relative abundance per gnha bin (global)
# Prepping df
n_bins <- 8
global_breaks <- quantile(
  corn_combined_top_asv$gnha,
  probs = seq(0, 1, length.out = n_bins + 1),
  na.rm = TRUE
) %>%
  unique()

corn_combined_binned <- corn_combined_top_asv %>%
  mutate(
    gnha_bin = cut(
      gnha,
      breaks = global_breaks,
      include.lowest = TRUE,
      dig.lab = 6
    )
  ) %>%
  filter(!is.na(gnha_bin)) %>%
  group_by(prevalence_level, asv, gnha_bin) %>%
  summarise(mean_rel = mean(rel_abundance, na.rm = TRUE), .groups = "drop")

# Plot

p_binned <- ggplot(
  corn_combined_binned,
  aes(x = gnha_bin, y = mean_rel, group = asv, color = asv)
) +
  geom_line(linewidth = 0.4, alpha = 0.85) +
  geom_point(size = 1) +
  facet_grid(~prevalence_level, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title = "CORN: Mean Relative Abundance vs gnha (Global Quantile Bins)",
    x = "gnha (binned)",
    y = "Mean Relative Abundance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_binned


# List and save
corn_top_asv_n2o_plots <- list(
  scatter_facet = p_scatter,
  binned_lines = p_binned
)

purrr::iwalk(
  corn_top_asv_n2o_plots,
  ~ ggsave(
    filename = paste0("data/output/plots/", .y, ".png"),
    plot = .x,
    width = if (.y == "scatter_facet") 14 else 10,
    height = if (.y == "scatter_facet") 8 else 6,
    dpi = 600
  )
)
################################################################
# Scraps
# 4. Format prevalence labels for nicer facet strip titles
# fmt_prev <- function(x) {
#   # e.g. prevalent_90 -> "Prevalent ≥90%"
#   str_replace(x, "^prevalent_?(\\d+)$", "Prevalent ≥\\1%") %>%
#     str_replace("^prevalent_?(\\d+)", "Prevalent ≥\\1%")
# }
# combined_top <- combined_top %>%
#   mutate(prevalence_level_fmt = fmt_prev(prevalence_level))

# Siginificance
# alpha_level <- 0.05

# 1. Compute stats per facet (example: prevalence_level_fmt + asv)
# facet_stats <- corn_combined_top_asv %>%
#   group_by(prevalence_level, sampling_location, asv) %>%
#   summarise(
#     p_value = {
#       fit <- lm(rel_abundance ~ gnha, data = .)
#       anova(fit)[["Pr(>F)"]][1]
#     },
#     slope = {
#       fit <- lm(rel_abundance ~ gnha, data = .)
#       coef(fit)[["gnha"]]
#     },
#     n = n(),
#     .groups = "drop"
#   ) %>%
#   mutate(significant = p_value < alpha_level) %>%
#   group_by(prevalence_level, sampling_location) %>%
#   mutate(
#     p_adj_bh = p.adjust(p_value, method = "BH"),
#     sig_bh = p_adj_bh < 0.05
#   ) %>%
#   ungroup()

# Order levels so they match 'asv' facet order
# asv_levels <- sort(unique(corn_combined_top_asv$asv))
# strip_colors <- ifelse(asv_levels %in% facet_stats, "#d7601bff", "#f0f0f0")
