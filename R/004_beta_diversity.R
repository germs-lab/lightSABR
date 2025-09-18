#####################################################################
# Beta Diversity Analyses
#
# This script explores general community diversity patterns in SABR data set.
# We explore, richness, alpha and beta diversity using NMDS, PCoA, PERMANOVAS , etc.

# Author: Jaejin Lee & Bolívar Aponte Rolón
# Last modified: 2025-08-21
#####################################################################

# Load required libraries
source("R/utils/000_setup.R")


#--------------------------------------------------------
# Calculate Bray-Curtis distance and perform NMDS analysis
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


# Needs fine tuning

# Perform NMDS analysis (k=2 reduces to 2 dimensions)
# trymax=100 ensures convergence by trying up to 100 different random starts
nmds <- ordinate(
  sabr_2023_physeq,
  method = "NMDS",
  distance = "bray",
  trymax = 500
)

#--------------------------------------------------------
# Generate and display NMDS plot
#--------------------------------------------------------

# # Create NMDS plot with samples colored by Plant, shaped by Location, and faceted by Date
# nmds_plot <- plot_ordination(
#   sabr_2023_physeq,
#   nmds,
#   color = "plant",
#   shape = "sampling_location"
# ) +
#   geom_point(size = 3, alpha = 0.8) +
#   facet_wrap(~sampling_date, ncol = 2) + # Separate panels by Date
#   labs(
#     title = "NMDS of Microbial Communities by Date",
#     x = "NMDS1",
#     y = "NMDS2"
#   ) +
#   theme_minimal() +
#   theme(
#     legend.position = "right",
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )

# # Display the plot
# nmds_plot

# # Optionally, save the plot
# ggsave(
#   "data/output/plots/nmds_plot.png",
#   nmds_plot,
#   width = 10,
#   height = 8,
#   dpi = 300
# )

#--------------------------------------------------------
# Create a phyloseq object with only high-prevalence ASVs
#--------------------------------------------------------

test_list <- list(
  prevalent_90 = main_physeq_list$original_phyloseq_lists$prevalent_90,
  prevalent_80 = main_physeq_list$original_phyloseq_lists$prevalent_80,
  prevalent_70 = main_physeq_list$original_phyloseq_lists$prevalent_70,
  prevalent_60 = main_physeq_list$original_phyloseq_lists$prevalent_60,
  prevalent_30 = main_physeq_list$original_phyloseq_lists$prevalent_30
)


test_nmds <- BRCore::brc_nmds(
  asv_matrix = t(as(phyloseq::otu_table(sabr_2023_physeq), "matrix")),
  physeq = sabr_2023_physeq,
  ncores = 1,
  k = 2,
  trymax = 500
)


test_map <- purrr::imap(
  test_list,
  function(ps_obj, prev_label) {
    # Guard against empty or trivial objects
    if (phyloseq::ntaxa(ps_obj) < 2 || phyloseq::nsamples(ps_obj) < 2) {
      return(NULL)
    }
    mat <- as(phyloseq::otu_table(sabr_2023_physeq), "matrix")
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

# Suppose each element in test_map has something like res$points or res$ordination
coords_tbl <- purrr::imap_dfr(
  test_map,
  function(res, prev_label) {
    if (is.null(res)) {
      return(dplyr::tibble())
    }
    # Adjust slot access to whatever brc_nmds returns (example: res$scores or res$points)
    df <- as.data.frame(res$points)
    df$sample_id <- rownames(df)
    df$prevalence <- prev_label
    df
  }
)

# # Extract OTU table for high-prevalence ASVs
# asv_high_prev <- asv_table_data[high_prevalence_asvs, ]
# asv_high_prev <- otu_table(asv_high_prev, taxa_are_rows = TRUE)

# # Create a new phyloseq object with only high-prevalence ASVs
# ps_high_prev <- merge_phyloseq(ps_rel, asv_high_prev)

# #--------------------------------------------------------
# # NMDS analysis for high-prevalence ASVs
# #--------------------------------------------------------

# # Calculate Bray-Curtis distance for the high-prevalence dataset
# bc_dist_high_prev <- phyloseq::distance(ps_high_prev, method = "bray")

# # Perform NMDS for high-prevalence ASVs
# nmds_high_prev <- ordinate(
#   ps_high_prev,
#   method = "NMDS",
#   distance = "bray",
#   trymax = 100
# )

# # Create and display NMDS plot for highly prevalent ASVs
# nmds_plot_high_prev <- plot_ordination(
#   ps_high_prev,
#   nmds_high_prev,
#   color = "Plant",
#   shape = "Location"
# ) +
#   geom_point(size = 3, alpha = 0.8) +
#   facet_wrap(~Date, ncol = 2) +
#   labs(
#     title = "NMDS of Highly Prevalent OTUs",
#     x = "NMDS1",
#     y = "NMDS2"
#   ) +
#   theme_minimal() +
#   theme(
#     legend.position = "right",
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )

# print(nmds_plot_high_prev)

# #--------------------------------------------------------
# # Identify and analyze low-prevalence ASVs (<70% samples)
# #--------------------------------------------------------

# # Define threshold for low prevalence (70% of samples)
# threshold_u70 <- total_samples * 0.7
# low_prevalence_asvs_u70 <- names(asv_sample_counts[
#   asv_sample_counts < threshold_u70
# ])

# # Extract taxonomy for low-prevalence ASVs
# low_prevalence_taxa_u70 <- tax_table_data[low_prevalence_asvs_u70, ]

# # Display information about low-prevalence ASVs
# cat("ASVs found in <70% samples:", length(low_prevalence_asvs_u70), "\n")
# print(low_prevalence_taxa_u70)

# # Create a new phyloseq object with only low-prevalence ASVs
# asv_low_prev_u70 <- otu_table_data[, low_prevalence_asvs_u70]
# asv_low_prev_u70 <- otu_table(asv_low_prev_u70, taxa_are_rows = FALSE)
# ps_low_prev_u70 <- merge_phyloseq(ps_rel, asv_low_prev_u70)

# # NMDS for low-prevalence ASVs (<70%)
# bc_dist_low_prev_u70 <- phyloseq::distance(ps_low_prev_u70, method = "bray")
# nmds_low_prev_u70 <- ordinate(
#   ps_low_prev_u70,
#   method = "NMDS",
#   distance = "bray",
#   trymax = 100
# )

# # Create and display NMDS plot for low-prevalence ASVs (<70%)
# nmds_plot_low_prev_u70 <- plot_ordination(
#   ps_low_prev_u70,
#   nmds_low_prev_u70,
#   color = "Plant",
#   shape = "Location"
# ) +
#   geom_point(size = 3, alpha = 0.8) +
#   facet_wrap(~Date, ncol = 2) +
#   labs(
#     title = "NMDS of Low Prevalence OTUs (Present in <70% of Samples)",
#     x = "NMDS1",
#     y = "NMDS2"
#   ) +
#   theme_minimal() +
#   theme(
#     legend.position = "right",
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )

# print(nmds_plot_low_prev_u70)

# #--------------------------------------------------------
# # Identify and analyze very low-prevalence ASVs (<30% samples)
# #--------------------------------------------------------

# # Define threshold for very low prevalence (30% of samples)
# threshold_u30 <- total_samples * 0.3
# low_prevalence_asvs_u30 <- names(asv_sample_counts[
#   asv_sample_counts < threshold_u30
# ])

# # Extract taxonomy for very low-prevalence ASVs
# low_prevalence_taxa_u30 <- tax_table_data[low_prevalence_asvs_u30, ]

# # Display information about very low-prevalence ASVs
# cat("ASVs found in <30% samples:", length(low_prevalence_asvs_u30), "\n")
# print(low_prevalence_taxa_u30)

# # Create a new phyloseq object with only very low-prevalence ASVs
# asv_low_prev_u30 <- otu_table_data[, low_prevalence_asvs_u30]
# asv_low_prev_u30 <- otu_table(asv_low_prev_u30, taxa_are_rows = FALSE)
# ps_low_prev_u30 <- merge_phyloseq(ps_rel, asv_low_prev_u30)

# # NMDS for very low-prevalence ASVs (<30%)
# bc_dist_low_prev_u30 <- phyloseq::distance(ps_low_prev_u30, method = "bray")
# nmds_low_prev_u30 <- ordinate(
#   ps_low_prev_u30,
#   method = "NMDS",
#   distance = "bray",
#   trymax = 100
# )

# # Create and display NMDS plot for very low-prevalence ASVs (<30%)
# nmds_plot_low_prev_u30 <- plot_ordination(
#   ps_low_prev_u30,
#   nmds_low_prev_u30,
#   color = "Plant",
#   shape = "Location"
# ) +
#   geom_point(size = 3, alpha = 0.8) +
#   facet_wrap(~Date, ncol = 2) +
#   labs(
#     title = "NMDS of Low Prevalence OTUs (Present in <30% of Samples)",
#     x = "NMDS1",
#     y = "NMDS2"
#   ) +
#   theme_minimal() +
#   theme(
#     legend.position = "right",
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )

# print(nmds_plot_low_prev_u30)
