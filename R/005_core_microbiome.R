#####################################################################
# Corn Interrow/Row Analysis
#
# This script analyzes the core microbiome across different corn plots
# and locations, identifying ASVs that are consistently present across
# different combinations of samples.
#
# Author: Bolívar Aponte Rolón & Jaejin Lee
# Date: 2025-05-05
# Last modified: 2025-10-14
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
        ggpmisc::use_label(c("adj.R2", "p")),
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


scatter_plots_by_prev$full_community

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
# -----------------------
# #List and save
# #------------------------
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

# Pretty tables with gt

# Correlation direction per prevalence × ASV × location
corn_corr_stats <- corn_combined_rel_abund %>%
  filter(prevalence_level == "full_community") %>% # Filter prevalence levels to obtain the community you want
  group_by(prevalence_level, asv, sampling_location) %>%
  summarise(
    n = sum(is.finite(gnha) & is.finite(rel_abundance)),
    r = {
      mask <- is.finite(gnha) & is.finite(rel_abundance)
      x <- gnha[mask]
      y <- rel_abundance[mask]
      if (length(x) >= 3 && stats::var(x) > 0 && stats::var(y) > 0) {
        stats::cor(x, y)
      } else {
        NA_real_
      }
    },
    r2 = r^2,
    r2_adj = ifelse(
      n - 1L - 1 > 0,
      1 - (1 - r2) * (n - 1) / (n - 1L - 1),
      NA_real_
    ), # p=1 or 1L for simple regression/correlation
    p_value = {
      mask <- is.finite(gnha) & is.finite(rel_abundance)
      x <- gnha[mask]
      y <- rel_abundance[mask]
      if (length(x) >= 3 && stats::var(x) > 0 && stats::var(y) > 0) {
        stats::cor.test(
          x,
          y,
          alternative = "two.sided",
          method = "pearson"
        )$p.value
      } else {
        NA_real_
      }
    },
    arrow = dplyr::case_when(
      is.na(r) ~ "–",
      r > 0 ~ "▲",
      r < 0 ~ "▼",
      TRUE ~ "–"
    ),
    .groups = "drop"
  )

# Select significantly correlated ASV (mostly, positive)
corn_corr_significant <- corn_corr_stats %>%
  filter(p_value < 0.05)

# Base taxonomy table (distinct ASV taxonomy per prevalence)
tax_base <- corn_combined_rel_abund %>%
  filter(prevalence_level == "full_community") %>%
  distinct(
    prevalence_level,
    asv,
    kingdom,
    phylum,
    class,
    family,
    genus,
    species
  )

# Join taxonomy table and correlation table
corn_tax_with_corr <- corn_corr_significant %>%
  dplyr::left_join(tax_base, by = c("prevalence_level", "asv")) %>%
  relocate(n:arrow, .after = species, ) %>%
  #select(-prevalence_level) %>%
  arrange(sampling_location, fct_reorder(asv, parse_number(asv), .fun = min), r)

# Split into one data frame per prevalence; name from the data (avoids mismatches)
corn_tables_by_prev <- corn_tax_with_corr %>%
  group_split(prevalence_level, .keep = TRUE)

names(corn_tables_by_prev) <- corn_tax_with_corr %>%
  distinct(prevalence_level) %>%
  pull(prevalence_level) %>%
  as.character()

# Build nicely formatted gt tables with colored arrows
pos_arrow <- "#38d005ff"
neg_arrow <- "#ff0000"


gt_tables_by_prev <- imap(
  corn_tables_by_prev,
  function(df_prev, prev_label) {
    # Ensure arrow columns exist even if one location is absent
    df_prev %>%
      select(
        -prevalence_level,
        r
      ) %>%
      gt::gt(
        groupname_col = "sampling_location",
        row_group_as_column = TRUE
      ) %>%
      cols_label(
        asv = "ASV",
        sampling_location = "Sampling location",
        kingdom = "Kingdom",
        phylum = "Phylum",
        class = "Class",
        family = "Family",
        genus = "Genus",
        species = "Species",
        n = md("*n*"),
        r2 = "R\u00B2",
        r2_adj = "Adj. R\u00B2",
        p_value = "p-value",
        arrow = "Correlation with g N (N\u2082O) ha\u207B\u00B9"
      ) %>%
      fmt_number(columns = c(r2, r2_adj), decimals = 3) %>%
      fmt(
        columns = p_value,
        fns = function(x) {
          ifelse(
            is.na(x),
            NA_character_,
            ifelse(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3))
          )
        }
      ) %>%
      tab_header(
        title = paste0("ASV Taxonomy and Correlation Direction — ", prev_label),
        subtitle = md("sample *n* = 994; ASV *n* = 988")
      ) %>%
      tab_style(
        style = cell_text(weight = "bolder"),
        locations = cells_title(groups = c("title", "subtitle"))
      ) %>%
      tab_stubhead(label = "Sampling Location") %>%
      tab_style(
        style = cell_text(weight = "bold"),
        locations = list(
          cells_column_labels(columns = everything()),
          cells_stubhead()
        )
      ) %>%
      # Stub formatting
      tab_style(
        style = cell_text(color = "black", weight = "bolder"),
        locations = cells_row_groups()
      ) %>%
      tab_style(
        style = list(
          cell_fill(color = "#f78f8fff"),
          cell_text(color = "black")
        ),
        locations = cells_row_groups(matches("Interrow"))
      ) %>%
      tab_style(
        style = list(
          cell_fill(color = "#31878fff"),
          cell_text(color = "black", )
        ),
        locations = cells_row_groups(starts_with("Row"))
      ) %>%

      # Color the arrow columns
      tab_style(
        style = cell_text(size = px(26), weight = "bold"),
        locations = cells_body(columns = "arrow")
      ) %>%
      tab_style(
        style = cell_text(color = pos_arrow),
        locations = cells_body(columns = "arrow", rows = arrow == "▲")
      ) %>%
      tab_style(
        style = cell_text(color = neg_arrow),
        locations = cells_body(columns = "arrow", rows = arrow == "▼")
      ) %>%
      # Optional: center the arrow columns
      cols_align(align = "center", columns = c(arrow)) %>%
      # Optional: compact look
      opt_table_lines("none")
  }
)

# Print all tables
invisible(lapply(gt_tables_by_prev, print))


gtsave(
  gt_tables_by_prev$full_community,
  filename = "data/output/significant_asvs.html"
)

#-------------------------------------------------------------
# Select significantly correlated AVS for downtream BLASTing
#-------------------------------------------------------------

duplicates <- corn_tax_with_corr |>
  group_by(asv) |>
  filter(n() > 1) |>
  ungroup()

corn_signif_asvs_strings <- corn_tax_with_corr |>
  distinct(asv) |>
  pull(asv)

# Subset to phyloseq to FASTA file

corn_signif_physeq_fcm <- corn_prevalence_list$full_community %>%
  prune_taxa(taxa_names(.) %in% corn_signif_asvs_strings, .)


physeq2fasta <- function(physeq, seq_col) {
  physeq_df <- physeq %>%
    phyloseq::tax_table(.) %>%
    phyloseq::as.data.frame(.)

  asv_headers <- paste0(">", rownames(physeq_df)) # ASV names must be rownames

  asv_seqs <- physeq_df[[seq_col]]

  asv_fasta <- c(rbind(asv_headers, asv_seqs))

  return(list(
    asv_fasta = asv_fasta,
    asv_headers = asv_headers
  ))
}


corn_signif_asvs_fasta <- physeq2fasta(
  corn_signif_physeq_fcm,
  seq_col = "sequence"
)


write(
  corn_signif_asvs_fasta$asv_fasta,
  file.path(
    "data/output/processed/sequences/sabr_2023_corn_signif_asv_fcm.fa"
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
