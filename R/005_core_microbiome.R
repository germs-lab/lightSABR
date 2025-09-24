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
#####################################################################
# ------------------------------------------------------------
# Plot relative abundance of ASVs vs gnha across prevalence levels
# ------------------------------------------------------------

# Build combined long data frame
combined_rel_abund <- imap_dfr(
  corn_prevalence_list,
  ~ relab_prevalence_long(.x, prevalence_label = .y)
)

# Select top ASVs to avoid overplotting (by global mean relative abundance)
top_n_asvs <- 20
top_asvs <- combined_rel_abund %>%
  group_by(asv) %>%
  summarise(mean_rel = mean(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_rel)) %>%
  slice_head(n = top_n_asvs) %>%
  pull(asv)

combined_top <- combined_rel_abund %>%
  filter(asv %in% top_asvs)

# 4. Format prevalence labels for nicer facet strip titles
# fmt_prev <- function(x) {
#   # e.g. prevalent_90 -> "Prevalent ≥90%"
#   str_replace(x, "^prevalent_?(\\d+)$", "Prevalent ≥\\1%") %>%
#     str_replace("^prevalent_?(\\d+)", "Prevalent ≥\\1%")
# }
# combined_top <- combined_top %>%
#   mutate(prevalence_level_fmt = fmt_prev(prevalence_level))

# Scatter + smooth (lm) faceted by prevalence level (rows) and ASV (cols)
p_scatter <- ggplot(
  combined_top,
  aes(x = gnha, y = rel_abundance, color = sampling_location)
) +
  geom_point(alpha = 0.35, size = 0.9, ) +
  geom_smooth(
    aes(color = sampling_location),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    #color = "#d7301f",
    linewidth = 0.5
  ) +
  facet_grid(prevalence_level ~ asv, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_log10() +
  labs(
    title = "CORN: ASV Relative Abundance vs gnha Across Prevalence Thresholds",
    x = "log10(gnha)",
    y = "Relative Abundance",
    subtitle = sprintf(
      "Top %d ASVs by selected mean relative abundance",
      top_n_asvs
    )
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

p_scatter

# Aggregate view: Mean relative abundance per gnha bin (global)
# Prepping df
n_bins <- 8
global_breaks <- quantile(
  combined_top$gnha,
  probs = seq(0, 1, length.out = n_bins + 1),
  na.rm = TRUE
) %>%
  unique()

combined_binned <- combined_top %>%
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
  combined_binned,
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
