#####################################################################
# NMDS Analysis of Microbial Communities
#
# This script performs Non-Metric Multidimensional Scaling (NMDS)
# analysis on the phyloseq object to visualize the microbial community
# structure across different samples.
#
# Author: Jaejin Lee
# Modified by: Bolívar Aponte Rolón
# Date: 2025-05-05
#####################################################################

# Setup
source("R/utils/000_setup.R")

physeq <- readRDS(file = "data/output/processed/sabr_physeq_object.rds")
# Ensure the phyloseq object (physeq) is loaded before running this script

#--------------------------------------------------------
# Exploration of data set
#--------------------------------------------------------
## Checking sampling effort, read counts, etc.
colnames(tax_table(physeq)) # Display column names for reference
colnames(sample_data(physeq))

## Summaries
#--------------------------------------------------------
metagMisc::phyloseq_summary(physeq, more_stats = F, long = F) # metagMisc overall summary
microbiome::summarize_phyloseq(physeq) #Microbiome package summary

# How many samples with ASVs per phyla?
percent_phyla_clean <- phyloseq_ntaxa_by_tax(
  physeq,
  TaxRank = "phylum",
  relative = F,
  add_meta_data = F
) |>
  as.data.frame() |>
  mutate(sum = sum(N.OTU)) |>
  group_by(phylum) |>
  summarise(occurance_in_samples = n()) # Count the number of ASV per phylum


## Coverage
#--------------------------------------------------------
cover_chao <- phyloseq_coverage(physeq, correct_singletons = T, add_attr = T) #Coverage: Good-Turing frequency estimation (Chiu, Chao 2016)

## Prevalence plot
phyloseq_prevalence_plot(
  physeq,
  prev.trh = 0.5,
  taxcolor = "phylum",
  facet = TRUE,
  point_alpha = 0.7,
  showplot = T
)

## Average relative ASVs per host plant
ps_average <- phyloseq_average(
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

## Core ASV abundance
core_abundance(
  physeq@otu_table,
  detection = 0.1 / 100,
  prevalence = 50 / 100,
  include.lowest = FALSE
)

## Rare ASVs
rare_abundance(
  physeq@otu_table,
  detection = 0.1 / 100,
  prevalence = 50 / 100,
  include.lowest = FALSE
)

## Core taxa
core_tax <- core_members(
  physeq,
  detection = 0.1 / 100,
  prevalence = 50 / 100,
  include.lowest = FALSE
)

## Rare taxa
rare_tax <- rare_members(
  physeq,
  detection = 0.1 / 100,
  prevalence = 50 / 100,
  include.lowest = FALSE
)

# Reads per sample
names_list <- colnames(physeq@otu_table) #List of samples

reads <- readcount(physeq) |>
  as.data.frame() |>
  rownames_to_column(var = "sample_id") |>
  rename(n_seqs = "readcount(physeq)") |>
  dplyr::right_join(
    rownames_to_column(metadata, var = "sample_id"),
    by = "sample_id"
  ) |>
  filter(sample_id %in% names_list)

# Sum of reads per site per species
reads_sum <- reads |>
  group_by(plot, plant) |>
  filter(sample_id %in% names_list) |>
  summarise(n_seqs = sum(n_seqs))

# Visualizing the distribution of reads
ggplot(reads, aes(x = n_seqs)) +
  geom_histogram(binwidth = 10000, fill = "grey", color = "black") +
  coord_cartesian(xlim = c(0, 100000)) # Cut-off

ggplot(reads, aes(x = 1, y = n_seqs)) +
  geom_jitter() +
  scale_y_log10()

test <- reads |>
  arrange(n_seqs) |>
  filter(sample_id %in% names_list)

ggplot(test, aes(x = 1:nrow(test), y = n_seqs)) +
  geom_line() +
  coord_cartesian(ylim = c(0, 100000)) # Cut-off


# Arrange the samples to see where the big breaks in counts are.
reads |>
  arrange(n_seqs) |>
  filter(sample_id %in% names_list)

#Good's coverage
ps_melt <- physeq |>
  psmelt() |>
  janitor::clean_names() |>
  rename(sample_id = sample, asv = otu) |>
  select(!nitrogen_conc)


cover_goods <- ps_melt |>
  select(sample_id, plant, ASV, Abundance) |>
  group_by(sample_id) |>
  summarise(
    n_seqs = sum(Abundance),
    n_sing = sum(Abundance == 1),
    goods = 1 - (n_sing / n_seqs)
  )
cover_goods |>
  filter(n_seqs > 750) |> #Filtering samples with more than f000 reads
  ggplot(aes(x = n_seqs, y = goods)) +
  geom_point() # Similar to `phyloseq_coverage` results above but here we know how it is calculated.

cover_goods |>
  filter(n_seqs > 750) |>
  arrange(goods)


#------------------------------------------------------

plot_bar(physeq, x = "samples", fill = "phylum")


## Richness

plot_richness(
  physeq,
  x = "samples",
  measures = c("Shannon", "Simpson"),
  color = "sampling_location",
  scales = "free"
)


#--------------------------------------------------------
# Calculate Bray-Curtis distance and perform NMDS analysis
#--------------------------------------------------------

# Calculate Bray-Curtis distance matrix
bc_dist <- phyloseq::distance(physeq, method = "bray")

# Perform NMDS analysis (k=2 reduces to 2 dimensions)
# trymax=100 ensures convergence by trying up to 100 different random starts
nmds <- ordinate(physeq, method = "NMDS", distance = "bray", trymax = 100)

#--------------------------------------------------------
# Prepare data for visualization
#--------------------------------------------------------

# Convert "Date" to a factor for proper panel separation in visualization
#sample_data(physeq)$Date <- as.factor(sample_data(physeq)$Date)

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
