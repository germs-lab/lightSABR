#####################################################################
# Core Microbiome Analysis
#
# This script analyzes the core microbiome across different corn plots
# and locations, identifying ASVs that are consistently present across
# different combinations of samples.
#
# Author: Bolívar Aponte Rolón
# Date: 2025-05-05
#####################################################################

# Load required libraries
library(phyloseq)
library(ggplot2)

# Ensure ps_rel (phyloseq object with relative abundances) is loaded

#--------------------------------------------------------
# Subset data by corn plant, plots, and locations
#--------------------------------------------------------

# Subset to include only Corn samples
ps_corn <- subset_samples(ps_rel, Plant == "Corn")
ps_corn
sample_data(ps_corn)
otu_table(ps_corn)
tax_table(ps_corn)

#--------------------------------------------------------
# Create phyloseq objects for each Corn plot
#--------------------------------------------------------

# Subset by plot
ps_corn_P01 <- subset_samples(ps_corn, Plot == "P01")
ps_corn_P08 <- subset_samples(ps_corn, Plot == "P08")
ps_corn_P11 <- subset_samples(ps_corn, Plot == "P11")
ps_corn_P14 <- subset_samples(ps_corn, Plot == "P14")

#--------------------------------------------------------
# Further subset by location (Interrow)
#--------------------------------------------------------

# Subset for Interrow locations
ps_corn_P01_I <- subset_samples(ps_corn_P01, Location == "Interrow")
ps_corn_P08_I <- subset_samples(ps_corn_P08, Location == "Interrow")
ps_corn_P11_I <- subset_samples(ps_corn_P11, Location == "Interrow")
ps_corn_P14_I <- subset_samples(ps_corn_P14, Location == "Interrow")

#--------------------------------------------------------
# Further subset by location (Row)
#--------------------------------------------------------

# Subset for Row locations
ps_corn_P01_R <- subset_samples(ps_corn_P01, Location == "Row")
ps_corn_P08_R <- subset_samples(ps_corn_P08, Location == "Row")
ps_corn_P11_R <- subset_samples(ps_corn_P11, Location == "Row")
ps_corn_P14_R <- subset_samples(ps_corn_P14, Location == "Row")

# Verify created objects
ps_corn_P01_I
ps_corn_P08_I
ps_corn_P11_I
ps_corn_P14_I
ps_corn_P01_R
ps_corn_P08_R
ps_corn_P11_R
ps_corn_P14_R

#--------------------------------------------------------
# Core microbiome analysis for P01 (Interrow)
#--------------------------------------------------------

# Extract ASV table for P01_I
asv_table_data_P01_I <- as.data.frame(otu_table(ps_corn_P01_I))
dim(asv_table_data_P01_I)

# Count number of samples in P01_I
total_samples_P01_I <- nrow(asv_table_data_P01_I)

# Count number of samples each ASV is present in
asv_sample_counts_P01_I <- colSums(asv_table_data_P01_I > 0)

# Identify ASVs present in >80% of P01_I samples
threshold_80_P01_I <- total_samples_P01_I * 0.8
core_asvs_P01_I <- names(asv_sample_counts_P01_I[
  asv_sample_counts_P01_I >= threshold_80_P01_I
])

cat("ASVs found in >80% P01_I samples:", length(core_asvs_P01_I), "\n")
head(asv_table_data_P01_I[, core_asvs_P01_I])

#--------------------------------------------------------
# Core microbiome analysis for P01 + P08 (Interrow)
#--------------------------------------------------------

# Merge P01_I and P08_I
ps_corn_P01_P08_I <- merge_phyloseq(ps_corn_P01_I, ps_corn_P08_I)

# Extract ASV table
asv_table_data_P01_P08_I <- as.data.frame(otu_table(ps_corn_P01_P08_I))

# Count total samples
total_samples_P01_P08_I <- nrow(asv_table_data_P01_P08_I)

# Count samples per ASV
asv_sample_counts_P01_P08_I <- colSums(asv_table_data_P01_P08_I > 0)

# Identify core ASVs (>80%)
threshold_80_P01_P08_I <- total_samples_P01_P08_I * 0.8
core_asvs_P01_P08_I <- names(asv_sample_counts_P01_P08_I[
  asv_sample_counts_P01_P08_I >= threshold_80_P01_P08_I
])

cat(
  "ASVs found in >80% P01_I + P08_I samples:",
  length(core_asvs_P01_P08_I),
  "\n"
)

# NMDS analysis for P01_I + P08_I
nmds_P01_P08_I <- ordinate(
  ps_corn_P01_P08_I,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

# Create NMDS plot
nmds_plot_P01_P08_I <- plot_ordination(
  ps_corn_P01_P08_I,
  nmds_P01_P08_I,
  color = "Plot",
  shape = "Date"
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "NMDS (P01-IR + P08-IR, 13 samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_P01_P08_I)

#--------------------------------------------------------
# Core microbiome analysis for P01 + P08 + P11 (Interrow)
#--------------------------------------------------------

# Merge P01_I, P08_I, and P11_I
ps_corn_P01_P08_P11_I <- merge_phyloseq(
  ps_corn_P01_I,
  ps_corn_P08_I,
  ps_corn_P11_I
)

# Extract ASV table
asv_table_data_P01_P08_P11_I <- as.data.frame(otu_table(ps_corn_P01_P08_P11_I))

# Count total samples
total_samples_P01_P08_P11_I <- nrow(asv_table_data_P01_P08_P11_I)

# Count samples per ASV
asv_sample_counts_P01_P08_P11_I <- colSums(asv_table_data_P01_P08_P11_I > 0)

# Identify core ASVs (>80%)
threshold_80_P01_P08_P11_I <- total_samples_P01_P08_P11_I * 0.8
core_asvs_P01_P08_P11_I <- names(asv_sample_counts_P01_P08_P11_I[
  asv_sample_counts_P01_P08_P11_I >= threshold_80_P01_P08_P11_I
])

cat(
  "ASVs found in >80% P01_I + P08_I + P11_I samples:",
  length(core_asvs_P01_P08_P11_I),
  "\n"
)

# NMDS analysis for P01_I + P08_I + P11_I
nmds_P01_P08_P11_I <- ordinate(
  ps_corn_P01_P08_P11_I,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

# Create NMDS plot
nmds_plot_P01_P08_P11_I <- plot_ordination(
  ps_corn_P01_P08_P11_I,
  nmds_P01_P08_P11_I,
  color = "Plot",
  shape = "Date"
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "NMDS (P01-IR + P08-IR + P11-IR, 20 samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_P01_P08_P11_I)

#--------------------------------------------------------
# Core microbiome analysis for all plots (Interrow)
#--------------------------------------------------------

# Merge all interrow samples
ps_corn_P01_P08_P11_P14_I <- merge_phyloseq(
  ps_corn_P01_I,
  ps_corn_P08_I,
  ps_corn_P11_I,
  ps_corn_P14_I
)

# Extract ASV table
asv_table_data_P01_P08_P11_P14_I <- as.data.frame(otu_table(
  ps_corn_P01_P08_P11_P14_I
))

# Count total samples
total_samples_P01_P08_P11_P14_I <- nrow(asv_table_data_P01_P08_P11_P14_I)

# Count samples per ASV
asv_sample_counts_P01_P08_P11_P14_I <- colSums(
  asv_table_data_P01_P08_P11_P14_I > 0
)

# Identify core ASVs (>80%)
threshold_80_P01_P08_P11_P14_I <- total_samples_P01_P08_P11_P14_I * 0.8
core_asvs_P01_P08_P11_P14_I <- names(asv_sample_counts_P01_P08_P11_P14_I[
  asv_sample_counts_P01_P08_P11_P14_I >= threshold_80_P01_P08_P11_P14_I
])

cat(
  "ASVs found in >80% P01_I + P08_I + P11_I + P14_I samples:",
  length(core_asvs_P01_P08_P11_P14_I),
  "\n"
)

# NMDS analysis for all interrow samples
nmds_P01_P08_P11_P14_I <- ordinate(
  ps_corn_P01_P08_P11_P14_I,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

# Create NMDS plot
nmds_plot_P01_P08_P11_P14_I <- plot_ordination(
  ps_corn_P01_P08_P11_P14_I,
  nmds_P01_P08_P11_P14_I,
  color = "Plot",
  shape = "Date"
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "NMDS (P01-IR + P08-IR + P11-IR + P14-IR, 27 samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_P01_P08_P11_P14_I)

#--------------------------------------------------------
# Core microbiome with Interrow + P01 Row
#--------------------------------------------------------

# Merge all interrow samples and P01 row samples
ps_corn_I_P01_R <- merge_phyloseq(
  ps_corn_P01_I,
  ps_corn_P08_I,
  ps_corn_P11_I,
  ps_corn_P14_I,
  ps_corn_P01_R
)

# Extract ASV table
asv_table_data_corn_I_P01_R <- as.data.frame(otu_table(ps_corn_I_P01_R))

# Count total samples
total_samples_corn_I_P01_R <- nrow(asv_table_data_corn_I_P01_R)

# Count samples per ASV
asv_sample_counts_corn_I_P01_R <- colSums(asv_table_data_corn_I_P01_R > 0)

# Identify core ASVs (>80%)
threshold_80_corn_I_P01_R <- total_samples_corn_I_P01_R * 0.8
core_asvs_corn_I_P01_R <- names(asv_sample_counts_corn_I_P01_R[
  asv_sample_counts_corn_I_P01_R >= threshold_80_corn_I_P01_R
])

cat("ASVs found in >80% corn_I+P01_RI", length(core_asvs_corn_I_P01_R), "\n")

# NMDS analysis
nmds_corn_I_P01_R <- ordinate(
  ps_corn_I_P01_R,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

# Create NMDS plot
nmds_plot_corn_I_P01_R <- plot_ordination(
  ps_corn_I_P01_R,
  nmds_corn_I_P01_R,
  color = "Plot",
  shape = "Date"
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "NMDS (Corn Interrow + P01-R, 33 samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_corn_I_P01_R)

#--------------------------------------------------------
# Expanding to include P08 Row samples
#--------------------------------------------------------

# Merge all interrow samples and P01, P08 row samples
ps_corn_I_P01_P08_R <- merge_phyloseq(
  ps_corn_P01_I,
  ps_corn_P08_I,
  ps_corn_P11_I,
  ps_corn_P14_I,
  ps_corn_P01_R,
  ps_corn_P08_R
)

# Extract ASV table
asv_table_data_corn_I_P01_P08_R <- as.data.frame(otu_table(ps_corn_I_P01_P08_R))

# Count total samples
total_samples_corn_I_P01_P08_R <- nrow(asv_table_data_corn_I_P01_P08_R)

# Count samples per ASV
asv_sample_counts_corn_I_P01_P08_R <- colSums(
  asv_table_data_corn_I_P01_P08_R > 0
)

# Identify core ASVs (>80%)
threshold_80_corn_I_P01_P08_R <- total_samples_corn_I_P01_P08_R * 0.8
core_asvs_corn_I_P01_P08_R <- names(asv_sample_counts_corn_I_P01_P08_R[
  asv_sample_counts_corn_I_P01_P08_R >= threshold_80_corn_I_P01_P08_R
])

cat(
  "ASVs found in >80% corn_Interrow + P01_R + P08_R samples:",
  length(core_asvs_corn_I_P01_P08_R),
  "\n"
)

# NMDS analysis
nmds_corn_I_P01_R_P08_R <- ordinate(
  ps_corn_I_P01_P08_R,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

# Create NMDS plot
nmds_plot_corn_I_P01_R_P08_R <- plot_ordination(
  ps_corn_I_P01_P08_R,
  nmds_corn_I_P01_R_P08_R,
  color = "Plot",
  shape = "Date"
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "NMDS (Corn Interrow + P01-R + P08-R, 40 samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_corn_I_P01_R_P08_R)

#--------------------------------------------------------
# Expanding to include P11 Row samples
#--------------------------------------------------------

# Merge all interrow samples and P01, P08, P11 row samples
ps_corn_I_P01_P08_P11_R <- merge_phyloseq(
  ps_corn_P01_I,
  ps_corn_P08_I,
  ps_corn_P11_I,
  ps_corn_P14_I,
  ps_corn_P01_R,
  ps_corn_P08_R,
  ps_corn_P11_R
)

# Extract ASV table
asv_table_data_corn_I_P01_P08_P11_R <- as.data.frame(otu_table(
  ps_corn_I_P01_P08_P11_R
))

# Count total samples
total_samples_corn_I_P01_P08_P11_R <- nrow(asv_table_data_corn_I_P01_P08_P11_R)

# Count samples per ASV
asv_sample_counts_corn_I_P01_P08_P11_R <- colSums(
  asv_table_data_corn_I_P01_P08_P11_R > 0
)

# Identify core ASVs (>80%)
threshold_80_corn_I_P01_P08_P11_R <- total_samples_corn_I_P01_P08_P11_R * 0.8
core_asvs_corn_I_P01_P08_P11_R <- names(asv_sample_counts_corn_I_P01_P08_P11_R[
  asv_sample_counts_corn_I_P01_P08_P11_R >= threshold_80_corn_I_P01_P08_P11_R
])

cat(
  "ASVs found in >80% corn_I+P01_P08_P11_R samples:",
  length(core_asvs_corn_I_P01_P08_P11_R),
  "\n"
)

# NMDS analysis
nmds_corn_I_P01_R_P08_R_P11_R <- ordinate(
  ps_corn_I_P01_P08_P11_R,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

# Create NMDS plot
nmds_plot_corn_I_P01_R_P08_R_P11_R <- plot_ordination(
  ps_corn_I_P01_P08_P11_R,
  nmds_corn_I_P01_R_P08_R_P11_R,
  color = "Plot",
  shape = "Date"
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "NMDS (Corn Interrow + P01-R + P08-R + P11-R, 47 samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_corn_I_P01_R_P08_R_P11_R)

#--------------------------------------------------------
# Including all corn samples (all interrow + all row)
#--------------------------------------------------------

# Merge all interrow and row samples
ps_corn_I_P01_P08_P11_P14_R <- merge_phyloseq(
  ps_corn_P01_I,
  ps_corn_P08_I,
  ps_corn_P11_I,
  ps_corn_P14_I,
  ps_corn_P01_R,
  ps_corn_P08_R,
  ps_corn_P11_R,
  ps_corn_P14_R
)

# Extract ASV table
asv_table_data_corn_I_P01_P08_P11_P14_R <- as.data.frame(otu_table(
  ps_corn_I_P01_P08_P11_P14_R
))

# Count total samples
total_samples_corn_I_P01_P08_P11_P14_R <- nrow(
  asv_table_data_corn_I_P01_P08_P11_P14_R
)

# Count samples per ASV
asv_sample_counts_corn_I_P01_P08_P11_P14_R <- colSums(
  asv_table_data_corn_I_P01_P08_P11_P14_R > 0
)

# Identify core ASVs (>80%)
threshold_80_corn_I_P01_P08_P11_P14_R <- total_samples_corn_I_P01_P08_P11_P14_R *
  0.8
core_asvs_corn_I_P01_P08_P11_P14_R <- names(asv_sample_counts_corn_I_P01_P08_P11_P14_R[
  asv_sample_counts_corn_I_P01_P08_P11_P14_R >=
    threshold_80_corn_I_P01_P08_P11_P14_R
])

cat(
  "ASVs found in >80% corn_I+R:",
  length(core_asvs_corn_I_P01_P08_P11_P14_R),
  "\n"
)

# NMDS analysis
nmds_corn_I_R <- ordinate(
  ps_corn_I_P01_P08_P11_P14_R,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

# Create NMDS plot
nmds_plot_corn_I_R <- plot_ordination(
  ps_corn_I_P01_P08_P11_P14_R,
  nmds_corn_I_R,
  color = "Plot",
  shape = "Date"
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "NMDS (Corn Interrow + Row, 53 samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_corn_I_R)

#--------------------------------------------------------
# Compare core ASVs with different thresholds (70%)
#--------------------------------------------------------

# Re-define thresholds to 70% for comparison
threshold_80_P01_I <- total_samples_P01_I * 0.7
core_asvs_P01_I <- names(asv_sample_counts_P01_I[
  asv_sample_counts_P01_I >= threshold_80_P01_I
])
cat("ASVs found in >70% P01_I samples:", length(core_asvs_P01_I), "\n")

threshold_80_P01_P08_I <- total_samples_P01_P08_I * 0.7
core_asvs_P01_P08_I <- names(asv_sample_counts_P01_P08_I[
  asv_sample_counts_P01_P08_I >= threshold_80_P01_P08_I
])
cat(
  "ASVs found in >70% P01_I + P08_I samples:",
  length(core_asvs_P01_P08_I),
  "\n"
)

threshold_80_P01_P08_P11_I <- total_samples_P01_P08_P11_I * 0.7
core_asvs_P01_P08_P11_I <- names(asv_sample_counts_P01_P08_P11_I[
  asv_sample_counts_P01_P08_P11_I >= threshold_80_P01_P08_P11_I
])
cat(
  "ASVs found in >70% P01_I + P08_I + P11_I samples:",
  length(core_asvs_P01_P08_P11_I),
  "\n"
)

threshold_80_P01_P08_P11_P14_I <- total_samples_P01_P08_P11_P14_I * 0.7
core_asvs_P01_P08_P11_P14_I <- names(asv_sample_counts_P01_P08_P11_P14_I[
  asv_sample_counts_P01_P08_P11_P14_I >= threshold_80_P01_P08_P11_P14_I
])
cat(
  "ASVs found in >70% P01_I + P08_I + P11_I + P14_I samples:",
  length(core_asvs_P01_P08_P11_P14_I),
  "\n"
)

threshold_80_corn_I_P01_R <- total_samples_corn_I_P01_R * 0.7
core_asvs_corn_I_P01_R <- names(asv_sample_counts_corn_I_P01_R[
  asv_sample_counts_corn_I_P01_R >= threshold_80_corn_I_P01_R
])
cat("ASVs found in >70% corn_I+P01_R", length(core_asvs_corn_I_P01_R), "\n")

threshold_80_corn_I_P01_P08_R <- total_samples_corn_I_P01_P08_R * 0.7
core_asvs_corn_I_P01_P08_R <- names(asv_sample_counts_corn_I_P01_P08_R[
  asv_sample_counts_corn_I_P01_P08_R >= threshold_80_corn_I_P01_P08_R
])
cat(
  "ASVs found in >70% corn_Interrow + P01_R + P08_R samples:",
  length(core_asvs_corn_I_P01_P08_R),
  "\n"
)

threshold_80_corn_I_P01_P08_P11_R <- total_samples_corn_I_P01_P08_P11_R * 0.7
core_asvs_corn_I_P01_P08_P11_R <- names(asv_sample_counts_corn_I_P01_P08_P11_R[
  asv_sample_counts_corn_I_P01_P08_P11_R >= threshold_80_corn_I_P01_P08_P11_R
])
cat(
  "ASVs found in >70% corn_I+P01_P08_P11_R samples:",
  length(core_asvs_corn_I_P01_P08_P11_R),
  "\n"
)

threshold_80_corn_I_P01_P08_P11_P14_R <- total_samples_corn_I_P01_P08_P11_P14_R *
  0.7
core_asvs_corn_I_P01_P08_P11_P14_R <- names(asv_sample_counts_corn_I_P01_P08_P11_P14_R[
  asv_sample_counts_corn_I_P01_P08_P11_P14_R >=
    threshold_80_corn_I_P01_P08_P11_P14_R
])
cat(
  "ASVs found in >70% corn_I+R:",
  length(core_asvs_corn_I_P01_P08_P11_P14_R),
  "\n"
)
