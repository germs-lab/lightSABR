-
  #####################################################################
  # Beta Diversity Analyses
  #
  # This script explores general community diversity patterns in SABR data set.
  # We explore, richness, alpha and beta diversity using NMDS, PCoA, PERMANOVAS , etc.

  # Author: Bolívar Aponte Rolón
  # Last modified: 2025-09-19
  #####################################################################

  # Load required libraries
  source("R/utils/000_setup.R")

#--------------------------------------------------------
# Section 1: Calculate Bray-Curtis distance and perform NMDS analysis
#--------------------------------------------------------

# Calculate Bray-Curtis distance matrix

main_raw_nmds <- BRCore::brc_nmds(
  asv_matrix = t(as(phyloseq::otu_table(sabr_2023_physeq), "matrix")), # Samples as rows
  physeq = sabr_2023_physeq,
  ncores = 1,
  k = 2,
  trymax = 500
)

# Run 111 stress 0.1878301
# ... Procrustes: rmse 5.135828e-05  max resid 0.0006460752
# ... Similar to previous best
# *** Best solution repeated 1 times

BRCore::brc_gg_ordi(.data = main_raw_nmds$nmds_df, .color = plant, "NMDS") +
  geom_point(
    aes(color = plant, shape = sampling_location),
    stroke = 1,
    alpha = 1,
    na.rm = TRUE,
    size = 2.5
  ) +
  facet_wrap(~sampling_date) +
  theme_minimal()


#--------------------------------------------------------
# Generate NMDS from phyloseq
# object at different thresholds
#--------------------------------------------------------

prevalence_list <- list(
  prevalent_100 = sabr_2023_physeq,
  prevalent_90 = main_physeq_list$original_phyloseq_lists$prevalent_90,
  prevalent_80 = main_physeq_list$original_phyloseq_lists$prevalent_80,
  prevalent_70 = main_physeq_list$original_phyloseq_lists$prevalent_70,
  prevalent_60 = main_physeq_list$original_phyloseq_lists$prevalent_60,
  prevalent_30 = main_physeq_list$original_phyloseq_lists$prevalent_30
)

# plot theme
nmds_theme <- list(
  geom_point(
    aes(color = plant, shape = sampling_location),
    stroke = 1,
    alpha = 1,
    na.rm = TRUE,
    size = 2.5
  ),
  facet_wrap(~sampling_date),
  theme_minimal()
)


prevalence_nmds_map <- purrr::imap(
  prevalence_list,
  function(ps_obj, prev_label) {
    mat <- as(phyloseq::otu_table(ps_obj), "matrix")
    # Ensure samples are rows (adjust if brc_nmds expects otherwise)
    if (phyloseq::taxa_are_rows(ps_obj)) {
      mat <- t(mat)
    }
    result <- BRCore::brc_nmds(
      asv_matrix = mat,
      physeq = ps_obj,
      ncores = 1,
      k = 2,
      trymax = 500
    )
    # Attach label
    attr(result, "prevalence_level") <- prev_label
    result
  }
)

# Verifying that all calculations are based on
# their corresponding matrices.
# Use to double check things. If all is good, skip.
# purrr::imap(prevalence_nmds_map, function(res, nm) {
#   cat("\n---", nm, "---\n")
#   if (is.null(res)) {
#     cat("Result is NULL\n")
#   } else {
#     print(head(res$nmds_scores))
#     cat("Stress:", res$ordi$stress, "\n")
#   }
#   invisible(NULL)
# })

# sapply(prevalence_list, function(ps) {
#   m <- as(phyloseq::otu_table(ps), "matrix")
#   if (phyloseq::taxa_are_rows(ps)) {
#     m <- t(m)
#   }
#   message(
#     ": dim=",
#     paste(dim(m), collapse = "x"),
#     " hash=",
#     digest::digest(m)
#   )
# })

# Plots should all be different

# Test plot
BRCore::brc_gg_ordi(
  .data = prevalence_nmds_map$prevalent_30$nmds_df,
  .color = plant,
  "NMDS"
) +
  nmds_theme
########################

nmds_df_list <- list(
  prevalent_100 = prevalence_nmds_map$prevalent_100$nmds_df,
  prevalent_90 = prevalence_nmds_map$prevalent_90$nmds_df,
  prevalent_80 = prevalence_nmds_map$prevalent_80$nmds_df,
  prevalent_70 = prevalence_nmds_map$prevalent_70$nmds_df,
  prevalent_60 = prevalence_nmds_map$prevalent_60$nmds_df,
  prevalent_30 = prevalence_nmds_map$prevalent_30$nmds_df
)

#--------------------------------------------------------
# Visualize NMDS plots
#--------------------------------------------------------

beta_nmds_plots <-
  imap(nmds_df_list, function(nmds_df, prev_label) {
    nmds_df <- nmds_df %>% mutate(prevalence_level = prev_label)
    # Nice title
    # Map "prevalent30" or "prevalent_30" -> "prevalent_30" to access the full result
    prev_key <- stringr::str_replace(
      prev_label,
      "^prevalent_?(\\d+).*",
      "prevalent_\\1"
    )

    # Percentages
    pct <- stringr::str_extract(prev_label, "\\d+")

    # ASVs used in the NMDS = number of species scores rows
    asv_count <- nrow(vegan::scores(
      prevalence_nmds_map[[prev_key]]$ordi,
      display = "species"
    ))

    title_label <- bquote(
      "ASVs Prevalent at " *
        .(pct) *
        "% occupancy (" *
        italic(n) ==
        .(asv_count) * ")"
    )

    p <- BRCore::brc_gg_ordi(
      .data = nmds_df,
      .color = plant,
      "NMDS"
    ) +
      nmds_theme +
      labs(title = title_label)

    # Attach label
    attr(p, "prevalence_level") <- prev_label
    p
  })

beta_nmds_plots

#--------------------------------------------------------
# Section 2: dbRDAs and PERMANOVAs
#--------------------------------------------------------
#--------------------------------------------------------
# dbRDA BC_CORE COMMUNITY ANALYSIS
#--------------------------------------------------------
# Transform data using Hellinger transformation
# Hellinger transformation

asv_matrices_list <- list(
  asvs_100pct = as(phyloseq::otu_table(sabr_2023_physeq), "matrix"),
  asvs_90pct = as(
    phyloseq::otu_table(main_physeq_list$original_phyloseq_lists$prevalent_90),
    "matrix"
  ),
  asvs_80pct = as(
    phyloseq::otu_table(main_physeq_list$original_phyloseq_lists$prevalent_80),
    "matrix"
  ),
  asvs_70pct = as(
    phyloseq::otu_table(main_physeq_list$original_phyloseq_lists$prevalent_70),
    "matrix"
  ),
  asvs_60pct = as(
    phyloseq::otu_table(main_physeq_list$original_phyloseq_lists$prevalent_60),
    "matrix"
  ),
  asvs_30pct = as(
    phyloseq::otu_table(main_physeq_list$original_phyloseq_lists$prevalent_30),
    "matrix"
  )
)

hellinger_matrices_list <- purrr::map(asv_matrices_list, function(asv_matrix) {
  hell_matrix <- decostand(asv_matrix, MARGIN = 1, method = "hellinger") #Rows are taxa

  hell_matrix_clean <- hell_matrix %>%
    as.data.frame() %>%
    select(where(~ is.numeric(.) && sum(.) > 0)) %>%
    # rownames_to_column(., var = "id") %>%
    # mutate(sample_id = id) %>%
    # column_to_rownames(., var = "id") %>%
    # relocate(., sample_id) %>%
    as.matrix()

  return(hell_matrix_clean)
})


# Prepare environmental data
dbrda_traits <- sabr_2023_metadata_clean %>%
  select(., c(plot:sampling_location, sampling_date, gnha:gwc_g_g)) %>%
  dplyr::filter(
    rownames(.) %in% colnames(hellinger_matrices_list$asvs_100pct)
  ) %>%
  rownames_to_column(., var = "id") %>%
  mutate(sample_id = id) %>%
  column_to_rownames(., var = "id") %>%
  relocate(., sample_id)

# Build dbRDA models with different constraints
dbrda_00 <- capscale(
  hellinger_matrices_list$asvs_90pct ~ 1,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 2,
  na.action = na.omit
)
str(dbrda_00)

dbrda_01_bc_core <- dbrda(
  t(hell_matrix) ~ .,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 8,
  na.action = na.omit
)

dbrda_02 <- dbrda(
  t(hellinger_matrices_list$asvs_90pct) ~ plot + plant,
  distance = "bray",
  dfun = vegdist,
  data = dbrda_traits,
  parallel = 2,
  na.action = na.omit
)
BRCore::brc_flex_ordi(
  dbrda_02,
  dbrda_traits,
  color_var = "plant",
  sample_id_col = "sample_id",
  biplot = TRUE,
  biplot_n = 18
)
