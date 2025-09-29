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
per_location_top <- corn_combined_rel_abund %>%
  group_by(prevalence_level, sampling_location, asv) %>%
  summarise(
    mean_rel = mean(rel_abundance, na.rm = TRUE),
    .groups = "drop_last"
  ) %>%
  slice_max(order_by = mean_rel, n = 20, with_ties = FALSE) %>%
  ungroup()

# Collapse to per-prevalence unique ASVs, rank by overall mean and cap at 20
corn_top20_pairs_union <- per_location_top %>%
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

# Data ready for plotting
nested_names <- corn_combined_top_asv %>%
  group_by(prevalence_level) %>%
  arrange(prevalence_level) %>%
  tidyr::nest()

dfs_by_prev <- corn_combined_top_asv %>%
  group_by(prevalence_level) %>%
  arrange(prevalence_level) %>%
  group_split(.keep = TRUE) %>%
  setNames(., nested_names$prevalence_level)

# Prepping df for binned plots
n_bins <- 8
global_breaks <- quantile(
  corn_combined_top_asv$gnha,
  probs = seq(0, 1, 0.25),
  na.rm = TRUE
) %>%
  unique()

dfs_by_prev_binned <- corn_combined_top_asv %>%
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
  summarise(
    mean_rel = mean(rel_abundance, na.rm = TRUE),
    .groups = "drop_last"
  ) %>%
  group_by(prevalence_level) %>%
  arrange(prevalence_level) %>%
  group_split(.keep = TRUE) %>%
  setNames(., nested_names$prevalence_level)

# Theming
axes_themes <- theme(
  plot.title = element_text(face = "bold"),
  strip.text = element_text(size = 5),
  axis.title = element_markdown(face = "bold"),
  axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
  panel.grid.minor = element_blank()
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
scatter_plots_by_prev <- dfs_by_prev %>%
  imap(function(df_prev, prev_label) {
    # Percentages
    pct <- stringr::str_extract(prev_label, "\\d+")

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
      facet_wrap(
        ~ fct_reorder(asv, parse_number(asv), .fun = min),
        scales = "free_y"
      ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      scale_x_log10(
        breaks = c(0.1, 1, 10, 100, 1000),
        limits = c(0.1, 1000),
        labels = scales::label_number()
      ) +
      labs(
        title = sprintf(
          "CORN: ASV Relative Abundance at %s%% prevalence threshold",
          pct
        ),
        x = "log<sub>10</sub>[g N (N<sub>2</sub>O) ha<sup>-1</sup>]", # "g N (N₂O) ha⁻¹"
        y = "Relative Abundance",
        subtitle = sprintf(
          "Top %d ASVs per location within each prevalence threshold selected by mean relative abundance",
          20
        )
      ) +
      theme_bw(base_size = 8) +
      axes_themes +
      theme(legend.position = "bottom")
  })


scatter_plots_by_prev

#------------------
# Binned plots
#------------------

# Aggregate view: Mean relative abundance per gnha bin (global)

binned_plots_by_prev <- dfs_by_prev_binned %>%
  imap(function(df_prev, prev_label) {
    # Percentages
    pct <- stringr::str_extract(prev_label, "\\d+")

    ggplot(
      df_prev,
      aes(x = gnha_bin, y = mean_rel, group = asv, color = asv)
    ) +
      geom_line(linewidth = 0.4, alpha = 0.85) +
      geom_point(size = 1) +
      # facet_grid(
      #   prevalence_level ~ fct_reorder(asv, parse_number(asv), .fun = min),
      #   scales = "free_y"
      # ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      labs(
        title = sprintf(
          "CORN: ASV  Mean Relative Abundance at %s%% prevalence threshold",
          pct
        ),
        x = "log<sub>10</sub>[g N (N<sub>2</sub>O) ha<sup>-1</sup>]", # "g N (N₂O) ha⁻¹"
        y = "Mean Relative Abundance",
        subtitle = sprintf(
          "Top %d ASVs per location within each prevalence threshold selected by mean relative abundance",
          20
        ),
        caption = "Global Quantiles Binned"
      ) +
      theme_bw() +
      axes_themes +
      guides(
        color = guide_legend(ncol = 2),
        fill = guide_legend(ncol = 2),
        shape = guide_legend(ncol = 2),
        linetype = guide_legend(ncol = 2)
      )
  })
binned_plots_by_prev

# # List and save
# corn_top_asv_n2o_plots <- list(
#   scatter_facet = scatter_plots_by_prev,
#   binned_lines = binned_plots_by_prev
# )

# safe_name <- function(x) str_replace_all(x, "[^A-Za-z0-9._-]+", "_")

# # recursive saver: handles ggplot/grob/gtable OR nested lists of them
# save_plot_or_list <- function(obj, base) {
#   if (inherits(obj, c("ggplot", "gg", "grob", "gtable"))) {
#     ggsave(
#       paste0(base, ".png"),
#       obj,
#       width = 14,
#       height = 8,
#       dpi = 600,
#       units = "in",
#       create.dir = TRUE
#     )
#     ggsave(
#       paste0(base, ".pdf"),
#       obj,
#       width = 14,
#       height = 8,
#       units = "in",
#       create.dir = TRUE
#     )
#   } else if (is.list(obj)) {
#     iwalk(
#       obj,
#       function(child, child_name) {
#         child_base <- file.path(
#           dirname(base),
#           paste0(basename(base), "_", safe_name(child_name))
#         )
#         save_plot_or_list(child, child_base)
#       }
#     )
#   } else {
#     warning(
#       "Skipping unsupported object type: ",
#       paste(class(obj), collapse = ", ")
#     )
#   }
# }

# purrr::iwalk(
#   corn_top_asv_n2o_plots,
#   ~ {
#     base <- file.path("data/output/plots", paste0(safe_name(.y)))
#     save_plot_or_list(.x, base)
#   }
# )

#--------------------------------------------------------
# Taxonomy tables
#--------------------------------------------------------

tables_by_prev <- corn_combined_top_asv %>%
  distinct(
    prevalence_level,
    asv,
    kingdom,
    phylum,
    class,
    # order,
    # family,
    genus,
    species
  ) %>%
  mutate(
    asv_num = readr::parse_number(asv),
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
  ) %>%
  arrange(
    prevalence_level,
    asv_num,
    asv
  ) %>%
  select(-asv_num) %>%
  group_split(prevalence_level, .keep = TRUE) %>%
  setNames(., nested_names$prevalence_level)

# tables_by_prev %>%
# group_by(prevalence_level) %>%
#   group_map(
#     ~ gt(.x %>% select(-prevalence_level)) %>%
#       tab_header(title = paste0("Prevalence: ", .y$prevalence_level)),
#     .keep = TRUE
#   ) -> gt_tables

# # Print all gt tables
# lapply(gt_tables, print)

purrr::iwalk(
  tables_by_prev,
  ~ {
    cat("\n\n")
    cat(sprintf("Prevalence: %s\n\n", .y))
    print(knitr::kable(.x %>% select(-prevalence_level)))
  }
)
# Maybe for prettier tables
# library(dplyr)
# library(tidyr)
# library(readr)
# library(purrr)
# library(gt)

# # 1) Choose colors to match your plots
# row_col <- "#1b9e77" # Row color
# interrow_col <- "#d95f02" # interrow color

# # 2) Correlation direction per prevalence × ASV × location
# corr_signs <- corn_combined_top_asv %>%
#   group_by(prevalence_level, asv, sampling_location) %>%
#   summarise(
#     r = suppressWarnings(cor(
#       gnha,
#       rel_abundance,
#       use = "pairwise.complete.obs"
#     )),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     arrow = dplyr::case_when(
#       is.na(r) ~ "–", # no data / undefined
#       r > 0 ~ "↑",
#       r < 0 ~ "↓",
#       TRUE ~ "–"
#     )
#   ) %>%
#   select(prevalence_level, asv, sampling_location, arrow) %>%
#   pivot_wider(
#     names_from = sampling_location,
#     values_from = arrow,
#     names_glue = "arrow_{sampling_location}"
#   )

# # 3) Base taxonomy table (distinct ASV taxonomy per prevalence)
# tax_base <- corn_combined_top_asv %>%
#   distinct(
#     prevalence_level,
#     asv,
#     kingdom,
#     phylum,
#     class,
#     family,
#     genus,
#     species
#   ) %>%
#   mutate(
#     asv_num = readr::parse_number(asv),
#     prevalence_level = factor(
#       prevalence_level,
#       levels = c(
#         "full_community",
#         "prevalent_30",
#         "prevalent_60",
#         "prevalent_70",
#         "prevalent_80",
#         "prevalent_90"
#       )
#     )
#   )

# # 4) Join arrows and order rows
# tax_with_corr <- tax_base %>%
#   left_join(corr_signs, by = c("prevalence_level", "asv")) %>%
#   arrange(prevalence_level, asv_num, asv)

# # 5) Split into one data frame per prevalence; name from the data (avoids mismatches)
# tables_by_prev <- tax_with_corr %>%
#   group_split(prevalence_level, .keep = TRUE)

# names(tables_by_prev) <- tax_with_corr %>%
#   distinct(prevalence_level) %>%
#   pull(prevalence_level) %>%
#   as.character()

# # 6) Build nicely formatted gt tables with colored arrows
# gt_tables_by_prev <- imap(
#   tables_by_prev,
#   function(df_prev, prev_label) {
#     # Ensure arrow columns exist even if one location is absent
#     if (!"arrow_Row" %in% names(df_prev)) {
#       df_prev$arrow_Row <- "–"
#     }
#     if (!"arrow_interrow" %in% names(df_prev)) {
#       df_prev$arrow_interrow <- "–"
#     }

#     df_prev %>%
#       select(
#         prevalence_level,
#         asv,
#         kingdom,
#         phylum,
#         class,
#         family,
#         genus,
#         species,
#         arrow_Row,
#         arrow_interrow
#       ) %>%
#       gt() %>%
#       cols_label(
#         prevalence_level = "Prevalence",
#         asv = "ASV",
#         kingdom = "Kingdom",
#         phylum = "Phylum",
#         class = "Class",
#         family = "Family",
#         genus = "Genus",
#         species = "Species",
#         arrow_Row = "Row",
#         arrow_interrow = "Interrow"
#       ) %>%
#       tab_header(
#         title = paste0("ASV Taxonomy and Correlation Direction — ", prev_label)
#       ) %>%
#       # Color the arrow columns
#       tab_style(
#         style = cell_text(color = row_col, weight = "bold"),
#         locations = cells_body(columns = "arrow_Row")
#       ) %>%
#       tab_style(
#         style = cell_text(color = interrow_col, weight = "bold"),
#         locations = cells_body(columns = "arrow_interrow")
#       ) %>%
#       # Optional: center the arrow columns
#       cols_align(align = "center", columns = c(arrow_Row, arrow_interrow)) %>%
#       # Optional: compact look
#       opt_table_lines("none")
#   }
# )

# # 7) Print all tables (in a Quarto report, this will render each table inline)
# invisible(lapply(gt_tables_by_prev, print))
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
