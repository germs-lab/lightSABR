#####################################################################
# Data Transformation and Export
#
# This script transforms the raw count data to relative abundance,
# performs rarefaction of sequences,examines basic properties of
# the dataset, and exports the transformed data for further analysis.
#
# Author: Jaejin Lee
# Modified by: Bolívar Aponte Rolón
# Date: 2025-05-05
#####################################################################

# Load required libraries
source("R/utils/000_setup.R")
# Ensure the phyloseq object (sabr_2023_physeq) is loaded before running this script

#--------------------------------------------------------
# Examine raw data properties and Normalize
#--------------------------------------------------------

#------------------------------
#  Relative Abundance
#------------------------------

# Display summary statistics of the raw count data
summary(as.vector(otu_table(sabr_2023_physeq)))

# Calculate and display the sum of counts for each sample
sample_sums(sabr_2023_physeq)

# Convert counts to relative abundance
sabr_2023_physeq_relab <- transform_sample_counts(
  sabr_2023_physeq,
  function(x) x / sum(x)
)
dim(otu_table(sabr_2023_physeq_relab))

# Verify (should all be 1)
sample_sums(sabr_2023_physeq_relab)

asv_table_relab_df <- as.data.frame(otu_table(sabr_2023_physeq_relab))

# Save the relative abundance table as a CSV file
write.csv(
  asv_table_rel_df,
  file = "data/output/processed/sabr_2023_asv_table_relab.csv",
  row.names = TRUE
)

save(
  sabr_2023_physeq_relab,
  file = "data/output/processed/sabr_2023_physeq_relab.rda"
)
#------------------------------
# Rarefaction
#------------------------------

# Rarefaction depth is informed by 002_eda_analyses.R, see `reads_sum`

asv_table_rrfy <- multi_rarefy(
  sabr_2023_physeq,
  depth_level = 5000,
  num_iter = 50,
  .summarize = FALSE,
  set_seed = 345
)

save(
  asv_table_rrfy,
  file = "data/output/processed/sabr_2023_asv_table_rrfy.rda"
)

# Let's make a rarefied data frame with it's corresponding metadata

#--------------------------------------------------------
# Rarefied Master Data Frame
# (ASVs and metadata, no taxonomical info)
#--------------------------------------------------------

taxa <- read.csv(file.path("data/output/processed/sabr_2023_taxonomy.csv")) %>%
  rename(., sequence = X) %>%
  rename_with(str_to_lower, .cols = everything()) # Clean up needed after importing from .csv
rownames(taxa) <- paste0("ASV_", 1:nrow(taxa))

# These should be loaded already by 000_setup.R
# load(file = "data/output/processed/sabr_2023_metadata_clean.rda")
# load(file = file = "data/output/processed/sabr_2023_asv_table_rrfy.rda")

# Master DF to match metadata to ASV iterations
mtr_rrfy_df <- asv_table_rrfy %>%
  rownames_to_column(., var = "iter_id") %>%
  dplyr::left_join(
    .,
    sabr_2023_metadata_clean %>% rownames_to_column(., var = "SampleID"),
    by = "SampleID"
  ) %>%
  column_to_rownames(., var = "iter_id") %>%
  relocate(., c(16881:16890), .after = "SampleID")

## Master metadata
mtr_metadata <- mtr_rrfy_df %>% # Metadata matched to all the samples in each iteration
  rownames_to_column(., var = "iter_id") %>%
  select(c(iter_id:nitrogen_conc)) %>%
  column_to_rownames(., var = "iter_id")

## Master ASV table, rarefied
mtr_asv <- mtr_rrfy_df %>%
  select(starts_with("ASV_")) %>%
  t()

# New, rarefied phyloseq object
mtr_physeq <- phyloseq(
  otu_table(as.matrix(mtr_asv), taxa_are_rows = TRUE),
  tax_table(as.matrix(taxa)),
  sample_data(mtr_metadata)
)


save(mtr_physeq, file = "data/output/processed/sabr_2023_mtr_physeq.rda")

#--------------------------------------------------------
# Calculate Diversity indices
#--------------------------------------------------------
mtr_rrfy_df <- mtr_rrfy_df %>%
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

save(mtr_rrfy_df, file = "data/output/processed/sabr_2023_master_rrfy_df.rda")
