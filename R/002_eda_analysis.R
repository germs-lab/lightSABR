#####################################################################
# Exploratory Data Analysis of Microbial Communities
#
# This script performs basi exploratory data analysis for micrbiome data set.
# Calculates summary statistic, read counts, Good's coverage, prevalence and core
# taxa analysis (albeit brief) and  Non-Metric Multidimensional Scaling (NMDS)
# analysis on the phyloseq object to visualize the microbial community
# structure across different samples.
#
# Author: Jaejin Lee
# Modified by: Bolívar Aponte Rolón
# Date: 2025-05-05
#####################################################################

# Setup
source("R/utils/000_setup.R")

physeq <- load(file = "data/output/processed/sabr_physeq_object.rda")

#--------------------------------------------------------
# Exploration of data set
#--------------------------------------------------------
## Basic metadata exploration
colnames(tax_table(physeq))
colnames(sample_data(physeq))

#--------------------------------------------------------
# Summaries and Read Counts
#--------------------------------------------------------
metagMisc::phyloseq_summary(physeq, more_stats = F, long = F)
microbiome::summarize_phyloseq(physeq)
# All samples have at least 5000 reads

# Taxonomic distribution
percent_phyla_clean <- phyloseq_ntaxa_by_tax(
  physeq,
  TaxRank = "phylum",
  relative = F,
  add_meta_data = F
) |>
  as.data.frame() |>
  mutate(sum = sum(N.OTU)) |>
  group_by(phylum) |>
  summarise(occurance_in_samples = n())

#--------------------------------------------------------
# Read Count Analysis
#--------------------------------------------------------
names_list <- colnames(physeq@otu_table)

reads <- readcount(physeq) |>
  as.data.frame() |>
  rownames_to_column(var = "sample_id") |>
  rename(n_seqs = "readcount(physeq)") |>
  dplyr::right_join(
    rownames_to_column(metadata, var = "sample_id"),
    by = "sample_id"
  ) |>
  filter(sample_id %in% names_list)

reads_sum <- reads |>
  group_by(sample_id) |>
  filter(sample_id %in% names_list) |>
  summarise(n_seqs = sum(n_seqs)) |>
  arrange(n_seqs)

#--------------------------------------------------------
# Read Count Visualizations
#--------------------------------------------------------
ggplot(reads, aes(x = n_seqs)) +
  geom_density(position = "identity", stat = "density", fill = "red") +
  #geom_histogram(binwidth = 250, fill = "grey", color = "black") +
  coord_cartesian(xlim = c(0, 40000))

ggplot(reads, aes(x = 1, y = n_seqs)) +
  geom_jitter() +
  scale_y_log10()

ggplot(
  reads |>
    arrange(n_seqs) |>
    filter(sample_id %in% names_list),
  aes(x = 1:nrow(reads), y = n_seqs)
) +
  geom_line() +
  coord_cartesian(ylim = c(0, 100000))

reads |>
  arrange(n_seqs) |>
  filter(sample_id %in% names_list)

#--------------------------------------------------------
# Good's Coverage Analysis
#--------------------------------------------------------
## Coverage estimates
cover_chao <- phyloseq_coverage(physeq, correct_singletons = T, add_attr = T)

ps_melt <- physeq |>
  psmelt() |>
  janitor::clean_names() |>
  rename(sample_id = sample, asv = otu) |>
  select(!nitrogen_conc)

cover_goods <- ps_melt |>
  select(sample_id, plant, asv, abundance) |>
  group_by(sample_id, plant) |>
  summarise(
    n_seqs = sum(abundance),
    n_singletons = sum(abundance == 1),
    goods = 1 - (n_singletons / n_seqs)
  )

cover_goods |>
  #filter(n_seqs > 750) |>
  ggplot(aes(x = n_seqs, y = goods, color = plant)) +
  geom_point()

cover_goods |>
  filter(n_seqs > 750) |>
  arrange(goods)


#--------------------------------------------------------
# Prevalence & Core Microbiome Analysis
#--------------------------------------------------------

## Prevalence visualization
metagMisc::phyloseq_prevalence_plot(
  physeq,
  prev.trh = 0.5,
  taxcolor = "phylum",
  facet = TRUE,
  point_alpha = 0.7,
  showplot = T
)

## Host plant-specific averages
ps_average <- metagMisc::phyloseq_average(
  physeq,
  avg_type = "arithmetic",
  acomp_zero_impute = NULL,
  aldex_samples = 213,
  aldex_denom = "all",
  group = "plant",
  drop_group_zero = TRUE,
  verbose = TRUE,
  progress = "text",
)

## Core microbiome analysis
microbiome::core_abundance(
  physeq@otu_table,
  detection = 0,
  prevalence = 50 / 100,
  include.lowest = FALSE
)

microbiome::rare_abundance(
  physeq@otu_table,
  detection = 0,
  prevalence = 50 / 100,
  include.lowest = FALSE
)

core_tax <- microbiome::core_members(
  physeq,
  detection = 0,
  prevalence = 90 / 100,
  include.lowest = FALSE
)

rare_tax <- microbiome::rare_members(
  physeq,
  detection = 0,
  prevalence = 50 / 100,
  include.lowest = FALSE
)
