#####################################################################
# SECTION 2: Data Transformation and Export
#
# This section of the script transforms the raw count data to
# relative abundance, performs rarefaction of sequences,examines basic
# properties of the dataset, and exports the transformed data for
# further analysis.
#
# Author: Bolívar Aponte Rolón
# Last modified: 2025-08-21
#####################################################################

# Setup
source("R/utils/000_setup.R")

# Ensure the phyloseq object (sabr_2023_physeq) is loaded.

#--------------------------------------------------------
# Main Data Frame
# (ASVs and metadata, no taxonomical info)
#--------------------------------------------------------

# Creating main_df: no data transformation,
# just manipulation for visualization and testing.

main_raw_df <- sabr_2023_physeq %>%
  otu_table() %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(., var = "iter_id") %>%
  dplyr::left_join(
    .,
    {
      sabr_2023_metadata_clean %>% rownames_to_column(., var = "sample_id")
    },
    by = c("iter_id" = "sample_id")
  ) %>%
  mutate(sample_id = iter_id) %>%
  relocate(sample_id, .before = ASV_1) %>%
  column_to_rownames(., var = "iter_id") %>%
  relocate(., c(samples:gwc_g_g), .after = sample_id)

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
  asv_table_relab_df,
  file = "data/output/processed/asv_tables/derived/sabr_2023_asv_table_relab.csv",
  row.names = TRUE
)

save(
  sabr_2023_physeq_relab,
  file = "data/output/processed/rdata/phyloseq/sabr_2023_physeq_relab.rda"
)
#------------------------------
# Rarefaction
#------------------------------

# Rarefaction depth is informed by 002_eda_analyses.R, see `reads_sum`

asv_table_rarefied <- multi_rarefy(
  sabr_2023_physeq,
  depth_level = 5000,
  num_iter = 50,
  .summarize = FALSE,
  set_seed = 345
) %>%
  rename(., sample_id = SampleID)

save(
  asv_table_rarefied,
  file = "data/output/processed/asv_tables/derived/sabr_2023_asv_table_rarefied.rda"
)

#--------------------------------------------------------
# Rarefied Main Data Frame
# (ASVs and metadata, no taxonomical info)
#--------------------------------------------------------

taxa <- read.csv(file.path(
  "data/output/processed/metadata/sabr_2023_taxonomy.csv"
)) %>%
  rename(., sequence = X) %>%
  rename_with(str_to_lower, .cols = everything()) # Clean up needed after importing from .csv

rownames(taxa) <- paste0("ASV_", 1:nrow(taxa))

# These should be loaded already by 000_setup.R
# load(file = "data/output/processed/sabr_2023_metadata_clean.rda")
# load(file = file = "data/output/processed/sabr_2023_asv_table_rrfy.rda")

# Main rarefied DF to match metadata to ASV iterations
main_rarefied_df <- asv_table_rarefied %>%
  rownames_to_column(., var = "iter_id") %>%
  dplyr::left_join(
    .,
    {
      sabr_2023_metadata_clean %>% rownames_to_column(., var = "sample_id")
    },
    by = "sample_id"
  ) %>%
  column_to_rownames(., var = "iter_id") %>%
  relocate(., c(16881:16900), .after = "sample_id")

## Main metadata
main_rarefied_metadata <- main_rarefied_df %>% # Metadata matched to all the samples in each iteration
  rownames_to_column(., var = "iter_id") %>%
  select(c(iter_id:gwc_g_g)) %>%
  column_to_rownames(., var = "iter_id")

## Main ASV table, rarefied
main_rarefied_asv <- main_rarefied_df %>%
  select(starts_with("ASV_")) %>%
  t()

# New, rarefied phyloseq object
main_rarefied_physeq <- phyloseq(
  otu_table(as.matrix(main_rarefied_asv), taxa_are_rows = TRUE),
  tax_table(as.matrix(taxa)),
  sample_data(main_rarefied_metadata)
)

save(
  main_rarefied_physeq,
  file = "data/output/processed/rdata/phyloseq/sabr_2023_main_rarefied_physeq.rda"
)

#--------------------------------------------------------
# Calculate Diversity indices
#--------------------------------------------------------
datasets <- list(
  raw = main_raw_df,
  rarefied = main_rarefied_df
)

unwanted_cols <- 1:21 # columns to keep as metadata; ASVs start after this

compute_diversity <- function(
  dataset,
  drop = unwanted_cols,
  first_asv_col = NULL
) {
  # Abundance block (ASVs) as numeric matrix
  abund <- dataset %>% select(-all_of(drop))
  abund_mat <- as.matrix(abund)

  out <- dataset %>%
    mutate(
      observed = rowSums(abund_mat > 0, na.rm = TRUE),
      shannon = vegan::diversity(abund_mat, index = "shannon"),
      simpson = vegan::diversity(abund_mat, index = "simpson"),
      inv_simpson = vegan::diversity(abund_mat, index = "invsimpson")
    )

  # Figure out where to relocate the new columns
  if (is.null(first_asv_col)) {
    # put before the first ASV column after the dropped block
    first_asv_col <- colnames(dataset)[min(setdiff(
      seq_len(ncol(dataset)),
      drop
    ))]
  }

  out %>%
    relocate(
      any_of(c("observed", "shannon", "simpson", "inv_simpson")),
      .before = all_of(first_asv_col)
    )
}

main_dataframes <- imap(
  datasets,
  ~ compute_diversity(.x, drop = unwanted_cols, first_asv_col = "ASV_1")
)


main_raw_df <- main_dataframes$raw
main_rarefied_df <- main_dataframes$rarefied

save(
  main_raw_df,
  file = "data/output/processed/rdata/dataframes/sabr_2023_main_raw_df.rda"
)
save(
  main_rarefied_df,
  file = "data/output/processed/rdata/dataframes/sabr_2023_main_rarefied_df.rda"
)
