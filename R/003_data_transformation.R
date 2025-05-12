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

load(file = "data/output/processed/sabr_2023_physeq_object.rda")
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


#--------------------------------------------------------
# Rarefaction
#--------------------------------------------------------

physeq_rrfy <- multi_rarefy(
  physeq,
  depth_level = 5000,
  num_iter = 50,
  .summarize = FALSE
)

save(physeq_rrfy, file = "data/output/processed/sabr_2023_asv_table_rrafy.rda")

#Let's make a rarefied data frame with it's corresponding metadata

#--------------------------------------------------------
# Master Data Frame
# (ASVs and metadata, no taxonomical info)
#--------------------------------------------------------

taxa <- read.csv(file.path("data/output/processed/sabr_2023_taxonomy.csv")) %>%
  rename(., sequence = X) %>%
  rename_with(str_to_lower, .cols = everything()) # Clean up needed after importing from .csv

rownames(taxa) <- paste0("ASV_", 1:nrow(taxa))

load(file = "data/output/processed/sabr_2023_metadata_clean.rda")
load(file = file = "data/output/processed/sabr_2023_asv_table_rrafy.rda")

mtr_physeq <- physeq_rrfy %>%
  mutate(row_id = paste0(SampleID, sep = "_", row_number())) %>%
  relocate(., row_id, .before = SampleID) %>%
  column_to_rownames(., var = "row_id") %>%
  dplyr::left_join(
    .,
    metadata %>% rownames_to_column(., var = "SampleID"),
    by = "SampleID"
  ) %>%
  relocate(., c(16881:16890), .after = "SampleID")

#--------------------------------------------------------
# Calculate Diversity indices
#--------------------------------------------------------
mtr_rfy_physeq <- mtr_physeq %>%
  mutate(
    observed = rowSums(select(., -c(1:11)) > 0),
    shannon = vegan::diversity(select(., -c(1:11)), index = "shannon"),
    simpson = vegan::diversity(select(., -c(1:11)), index = "simpson"),
    invsimpson = vegan::diversity(select(., -c(1:11)), index = "invsimpson")
  ) %>%
  relocate(
    any_of(c("observed", "shannon", "simpson", "invsimpson")),
    .before = ASV_1
  )

save(mtr_rfy_physeq, file = "data/output/processed/sabr_2023_master_rfy_df.rda")
load(file = "data/output/processed/sabr_2023_master_rfy_df.rda")
#--------------------------------------------------------
# SCRAPS for later
#--------------------------------------------------------

#dplyr::left_join(., metadata %>% rownames_to_column(., var = "SampleID"))

# t() %>%
#   as.data.frame() %>%
#   janitor::row_to_names(row_number = 1) %>%
#   mutate(across(where(is.character), as.numeric))

# new_asv <- physeq_rrfy %>%
#   mutate(row_id = paste0(SampleID, sep = "_", row_number())) %>%
#   relocate(., row_id, .before = SampleID) %>%
#   column_to_rownames(., var = "row_id") %>%
#   select(!SampleID) %>%
#   as.matrix() %>%
#   t() %>%
#   as.data.frame() %>%
#   rownames_to_column(., var = "SampleID")

# # Phyloseq object
# physeq_rrfy_obj <- phyloseq(
#   otu_table(t(new_asv), taxa_are_rows = TRUE),
#   tax_table(as.matrix(taxa)),
#   sample_data(metadata)
# )
