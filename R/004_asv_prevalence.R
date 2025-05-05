#####################################################################
# ASV Prevalence Analysis
#
# This script analyzes the prevalence of ASVs across samples,
# identifying and examining high-prevalence ASVs that are present
# in a large proportion of samples.
#
# Author: Bolívar Aponte Rolón
# Date: 2025-05-05
#####################################################################

# Load required libraries
library(phyloseq)
library(ggplot2)

# Ensure ps_rel (phyloseq object with relative abundances) is loaded

#--------------------------------------------------------
# Calculate ASV prevalence across samples
#--------------------------------------------------------

# Extract OTU table as data frame
asv_table_data <- as.data.frame(otu_table(ps_rel))

# Count how many samples contain each ASV (presence > 0)
asv_sample_counts <- colSums(asv_table_data > 0)

# Create a data frame for visualization
asv_count_df <- data.frame(
  OTU = names(asv_sample_counts),
  Sample_Counts = asv_sample_counts
)

#--------------------------------------------------------
# Visualize ASV prevalence
#--------------------------------------------------------

# Create a bar plot showing the number of samples each ASV is found in
ggplot(asv_count_df, aes(x = reorder(OTU, -Sample_Counts), y = Sample_Counts)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Number of Samples Each OTU is Found In",
    x = "ASV",
    y = "Number of Samples"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), # Remove x-axis labels as there are too many ASVs
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

#--------------------------------------------------------
# Identify high-prevalence ASVs (>90% samples)
#--------------------------------------------------------

# Get total number of samples
total_samples <- nrow(otu_table_data)
total_samples

# Define threshold for high prevalence (90% of samples)
threshold <- total_samples * 0.9
high_prevalence_asvs <- names(asv_sample_counts[asv_sample_counts >= threshold])

# Display information about high-prevalence ASVs
cat("ASVs found in >90% samples:", length(high_prevalence_asvs), "\n")
cat("List of ASVs found in >90% samples:\n", high_prevalence_asvs, "\n")

#--------------------------------------------------------
# Extract taxonomy for high-prevalence ASVs
#--------------------------------------------------------

# Get taxonomy table
tax_table_data <- tax_table(ps_rel)

# Extract taxonomy for high-prevalence ASVs
high_prevalence_taxa <- tax_table_data[high_prevalence_asvs, ]
print(high_prevalence_taxa)

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
