#####################################################################
# Community Diversity Analyses
#
# This script explores general community diversity patterns in SABR data set.
# We explore, richness, alpha and beta diversity using NMDS, PCoA, PERMANOVAS , etc.
# Author: Jaejin Lee
# Modified by: Bolívar Aponte Rolón
# Date: 2025-05-05
#####################################################################

# Load required libraries
source("R/utils/000_setup.R")

#------------------------------------------------------
# Community Diversity
#--------------------------------------------------------

#plot_bar(physeq, x = "samples", fill = "phylum")

## Richness
plot_richness(
  mtr_physeq,
  x = "plot",
  measures = c("Observed", "Shannon", "Simpson"),
  color = "plant",
  shape = "sampling_location",
  scales = "free"
)

# Reshape the data for ggplot (long format for diversity indices)
mtr_rrfy_df_long <- mtr_rrfy_df %>%
  tidyr::pivot_longer(
    cols = c(observed, shannon, simpson, invsimpson),
    names_to = "Diversity_Index",
    values_to = "Value"
  )

ggplot(mtr_rrfy_df_long, aes(x = sampling_location, y = Value)) +
  geom_boxplot(alpha = 0, width = 0.6) +
  geom_jitter(aes(color = sampling_location), width = 0.2, alpha = 0.03) +
  facet_wrap(~Diversity_Index, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Diversity Indices Across Sampling Locations",
    x = "Sampling Location",
    y = "Diversity Index Value",
    color = "Sampling Location"
  )
#--------------------------------------------------------
# Calculate Bray-Curtis distance and perform NMDS analysis
#--------------------------------------------------------

# Calculate Bray-Curtis distance matrix
bc_dist <- phyloseq::distance(mtr_physeq, method = "bray")

# Perform NMDS analysis (k=2 reduces to 2 dimensions)
# trymax=100 ensures convergence by trying up to 100 different random starts
nmds <- ordinate(mtr_physeq, method = "NMDS", distance = "bray", trymax = 100)

#--------------------------------------------------------
# Generate and display NMDS plot
#--------------------------------------------------------

# Create NMDS plot with samples colored by Plant, shaped by Location, and faceted by Date
nmds_plot <- plot_ordination(
  physeq,
  nmds,
  color = "plant",
  shape = "sampling_location"
) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~sampling_date, ncol = 2) + # Separate panels by Date
  labs(
    title = "NMDS of Microbial Communities by Date",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Display the plot
nmds_plot

# Optionally, save the plot
ggsave(
  "data/output/plots/nmds_plot.png",
  nmds_plot,
  width = 10,
  height = 8,
  dpi = 300
)


# Much like phyloseq::transform_sample_counts()
tse <- transformAssay(ps_tse, method = "relabundance")
phy2 <- convertToPhyloseq(tse, assay.type = "relabundance")
#--------------------------------------------------------
# Create a phyloseq object with only high-prevalence ASVs
#--------------------------------------------------------

# Extract OTU table for high-prevalence ASVs
asv_high_prev <- asv_table_data[high_prevalence_asvs, ]
asv_high_prev <- otu_table(asv_high_prev, taxa_are_rows = TRUE)

# Create a new phyloseq object with only high-prevalence ASVs
ps_high_prev <- merge_phyloseq(ps_rel, asv_high_prev)

#--------------------------------------------------------
# NMDS analysis for high-prevalence ASVs
#--------------------------------------------------------

# Calculate Bray-Curtis distance for the high-prevalence dataset
bc_dist_high_prev <- phyloseq::distance(ps_high_prev, method = "bray")

# Perform NMDS for high-prevalence ASVs
nmds_high_prev <- ordinate(
  ps_high_prev,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

# Create and display NMDS plot for highly prevalent ASVs
nmds_plot_high_prev <- plot_ordination(
  ps_high_prev,
  nmds_high_prev,
  color = "Plant",
  shape = "Location"
) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~Date, ncol = 2) +
  labs(
    title = "NMDS of Highly Prevalent OTUs",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_high_prev)

#--------------------------------------------------------
# Identify and analyze low-prevalence ASVs (<70% samples)
#--------------------------------------------------------

# Define threshold for low prevalence (70% of samples)
threshold_u70 <- total_samples * 0.7
low_prevalence_asvs_u70 <- names(asv_sample_counts[
  asv_sample_counts < threshold_u70
])

# Extract taxonomy for low-prevalence ASVs
low_prevalence_taxa_u70 <- tax_table_data[low_prevalence_asvs_u70, ]

# Display information about low-prevalence ASVs
cat("ASVs found in <70% samples:", length(low_prevalence_asvs_u70), "\n")
print(low_prevalence_taxa_u70)

# Create a new phyloseq object with only low-prevalence ASVs
asv_low_prev_u70 <- otu_table_data[, low_prevalence_asvs_u70]
asv_low_prev_u70 <- otu_table(asv_low_prev_u70, taxa_are_rows = FALSE)
ps_low_prev_u70 <- merge_phyloseq(ps_rel, asv_low_prev_u70)

# NMDS for low-prevalence ASVs (<70%)
bc_dist_low_prev_u70 <- phyloseq::distance(ps_low_prev_u70, method = "bray")
nmds_low_prev_u70 <- ordinate(
  ps_low_prev_u70,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

# Create and display NMDS plot for low-prevalence ASVs (<70%)
nmds_plot_low_prev_u70 <- plot_ordination(
  ps_low_prev_u70,
  nmds_low_prev_u70,
  color = "Plant",
  shape = "Location"
) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~Date, ncol = 2) +
  labs(
    title = "NMDS of Low Prevalence OTUs (Present in <70% of Samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_low_prev_u70)

#--------------------------------------------------------
# Identify and analyze very low-prevalence ASVs (<30% samples)
#--------------------------------------------------------

# Define threshold for very low prevalence (30% of samples)
threshold_u30 <- total_samples * 0.3
low_prevalence_asvs_u30 <- names(asv_sample_counts[
  asv_sample_counts < threshold_u30
])

# Extract taxonomy for very low-prevalence ASVs
low_prevalence_taxa_u30 <- tax_table_data[low_prevalence_asvs_u30, ]

# Display information about very low-prevalence ASVs
cat("ASVs found in <30% samples:", length(low_prevalence_asvs_u30), "\n")
print(low_prevalence_taxa_u30)

# Create a new phyloseq object with only very low-prevalence ASVs
asv_low_prev_u30 <- otu_table_data[, low_prevalence_asvs_u30]
asv_low_prev_u30 <- otu_table(asv_low_prev_u30, taxa_are_rows = FALSE)
ps_low_prev_u30 <- merge_phyloseq(ps_rel, asv_low_prev_u30)

# NMDS for very low-prevalence ASVs (<30%)
bc_dist_low_prev_u30 <- phyloseq::distance(ps_low_prev_u30, method = "bray")
nmds_low_prev_u30 <- ordinate(
  ps_low_prev_u30,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

# Create and display NMDS plot for very low-prevalence ASVs (<30%)
nmds_plot_low_prev_u30 <- plot_ordination(
  ps_low_prev_u30,
  nmds_low_prev_u30,
  color = "Plant",
  shape = "Location"
) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~Date, ncol = 2) +
  labs(
    title = "NMDS of Low Prevalence OTUs (Present in <30% of Samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_low_prev_u30)
