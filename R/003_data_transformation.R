#####################################################################
# Data Transformation and Export
#
# This script transforms the raw count data to relative abundance,
# examines basic properties of the dataset, and exports the transformed
# data for further analysis.
#
# Author: Jaejin Lee
# Modified by: Bolívar Aponte Rolón
# Date: 2025-05-05
#####################################################################

# Load required libraries
source("R/utils/000_setup.R")
physeq <- readRDS(file = "data/output/processed/sabr_physeq_object.rds")
# Ensure the phyloseq object (physeq) is loaded before running this script

#--------------------------------------------------------
# Examine raw data properties
#--------------------------------------------------------

# Display summary statistics of the raw count data
summary(as.vector(otu_table(physeq)))

# Calculate and display the sum of counts for each sample
sample_sums(physeq)

#--------------------------------------------------------
# Transform to relative abundance
#--------------------------------------------------------

# Convert counts to relative abundance (proportion of each sample)
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# Verify transformation by checking sample sums (should all be 1)
sample_sums(physeq_rel)

#--------------------------------------------------------
# Extract and export relative abundance data
#--------------------------------------------------------

# Convert OTU table to a data frame
asv_table_rel_df <- as.data.frame(otu_table(physeq_rel))

# Save the relative abundance table as a CSV file
# Update the file path to match your directory structure
write.csv(
  asv_table_rel_df,
  file = "analysis/asv_table_physeq_rel.csv",
  row.names = TRUE
)

# Display dimensions of the OTU table
dim(otu_table(physeq_rel))
