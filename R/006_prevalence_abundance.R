#####################################################################
# Prevalence vs. Abundance Analysis
#
# This script analyzes the relationship between ASV prevalence
# (proportion of samples an ASV is found in) and abundance
# (relative abundance when present).
#
# Author: Bolívar Aponte Rolón
# Date: 2025-05-05
#####################################################################

# Load required libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(scales)

# Ensure ps_rel (phyloseq object with relative abundances) is loaded

#--------------------------------------------------------
# Extract and export all ASV data
#--------------------------------------------------------

# Extract ASV table
asv_table_data <- as.data.frame(otu_table(ps_rel))

# Save the complete ASV table
write.csv(
  asv_table_data,
  file = "analysis/asv_table_ps_rel_all.csv",
  row.names = TRUE
)

# Calculate total number of samples
total_samples <- nrow(asv_table_data)

#--------------------------------------------------------
# Calculate prevalence and abundance metrics
#--------------------------------------------------------

# Calculate prevalence (percentage of samples where ASV is present)
asv_prevalence <- colSums(asv_table_data > 0) / total_samples * 100

# Calculate mean abundance for each ASV
asv_abundance <- colMeans(asv_table_data)

# Create a data frame for analysis
prevalence_abundance_df <- data.frame(
  OTU = names(asv_prevalence),
  Prevalence = asv_prevalence,
  Abundance = asv_abundance
)

#--------------------------------------------------------
# Create prevalence vs. abundance plot with 90% threshold
#--------------------------------------------------------

# Create scatter plot with highlighted high-prevalence ASVs
ggplot(prevalence_abundance_df, aes(x = Prevalence, y = Abundance)) +
  geom_point(alpha = 0.6, color = "gray60") + # All ASVs (gray)
  geom_point(
    data = prevalence_abundance_df[prevalence_abundance_df$Prevalence >= 90, ],
    aes(x = Prevalence, y = Abundance),
    color = "blue",
    alpha = 0.8
  ) + # ASVs >90% prevalence (blue)
  scale_y_log10(labels = function(y) round(log10(y), 0)) + # Log scale for abundance
  geom_vline(xintercept = 90, linetype = "dashed", color = "red", size = 1) + # 90% line
  labs(
    title = "Prevalence vs. Abundance (All ASVs, 90% Core Highlighted)",
    x = "Prevalence (%)",
    y = "Mean Relative Abundance (Log scale)"
  ) +
  theme_minimal()

#--------------------------------------------------------
# Extract and analyze high-prevalence ASVs
#--------------------------------------------------------

# Identify ASVs with >90% prevalence
blue_asvs <- prevalence_abundance_df$OTU[
  prevalence_abundance_df$Prevalence >= 90
]

cat("ASVs >90% corn samples:", length(blue_asvs), "\n")

# Extract taxonomy data
tax_table_data <- as.data.frame(tax_table(ps_rel))

# Get taxonomy for high-prevalence ASVs
blue_asvs_taxa <- tax_table_data[blue_asvs, ]
print(blue_asvs_taxa)

# Create comprehensive data frame with taxonomy, prevalence, and abundance
blue_asvs_df <- blue_asvs_taxa %>%
  mutate(
    OTU_Sequence = rownames(blue_asvs_taxa), # ASV sequences
    Prevalence = asv_prevalence[blue_asvs], # Prevalence values
    Mean_Relative_Abundance = asv_abundance[blue_asvs] # Mean abundance values
  ) %>%
  select(OTU_Sequence, everything()) # Rearrange columns

print(blue_asvs_df)

# Export high-prevalence ASV data
write.csv(
  blue_asvs_df,
  file = "analysis/blue_asvs_table.csv",
  row.names = FALSE
)

#--------------------------------------------------------
# Advanced filtering: Both high prevalence and high abundance
#--------------------------------------------------------

# Calculate prevalence
asv_prevalence_1 <- colSums(asv_table_data > 0) / total_samples * 100

# Calculate abundance
asv_abundance_1 <- colMeans(asv_table_data)

# Define threshold for top 1% abundance
abundance_threshold_1 <- quantile(asv_abundance_1, 0.99)

# Create data frame
prevalence_abundance_df_1 <- data.frame(
  OTU = names(asv_prevalence_1),
  Prevalence = asv_prevalence_1,
  Abundance = asv_abundance_1
)

# Create plot highlighting ASVs that meet both criteria
ggplot(prevalence_abundance_df_1, aes(x = Prevalence, y = Abundance)) +
  geom_point(alpha = 0.6, color = "gray60") + # All ASVs (gray)
  geom_point(
    data = prevalence_abundance_df_1[
      prevalence_abundance_df_1$Prevalence >= 80 &
        prevalence_abundance_df_1$Abundance >= abundance_threshold_1,
    ],
    aes(x = Prevalence, y = Abundance),
    color = "blue",
    alpha = 0.8
  ) + # ASVs meeting both criteria (blue)
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(.x))
  ) + # Log scale
  geom_vline(xintercept = 80, linetype = "dashed", color = "red", size = 1) + # Prevalence threshold
  geom_hline(
    yintercept = abundance_threshold_1,
    linetype = "dashed",
    color = "green",
    size = 1
  ) + # Abundance threshold
  labs(
    title = "Prevalence vs. Abundance (80% Prevalence & Top 1% Abundance Highlighted)",
    x = "Prevalence (%)",
    y = "Mean Relative Abundance (Log scale)"
  ) +
  theme_minimal()

#--------------------------------------------------------
# Extract ASVs meeting both criteria
#--------------------------------------------------------

# Identify ASVs that meet both criteria
selected_otus_1 <- prevalence_abundance_df_1$OTU[
  prevalence_abundance_df_1$Prevalence >= 80 &
    prevalence_abundance_df_1$Abundance >= abundance_threshold_1
]

# Count the number of selected ASVs
cat(
  "OTUs meeting both criteria (Prevalence ≥ 80% & Top 1% Abundance):",
  length(selected_otus_1),
  "\n"
)

# Get taxonomy for selected ASVs
selected_otus_taxa_1 <- tax_table_data[selected_otus_1, , drop = FALSE]

# Create comprehensive data frame
selected_otus_df_1 <- selected_otus_taxa_1 %>%
  mutate(
    OTU_Sequence = rownames(selected_otus_taxa_1),
    Prevalence = asv_prevalence_1[selected_otus_1],
    Mean_Relative_Abundance = asv_abundance_1[selected_otus_1]
  ) %>%
  select(OTU_Sequence, everything())

# Display the results
print(selected_otus_df_1)

# Export the data
write.csv(
  selected_otus_df_1,
  file = "analysis/selected_otus_table_1.csv",
  row.names = FALSE
)
