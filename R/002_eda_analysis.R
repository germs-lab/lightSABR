##################################################################################
# Exploratory Data Analysis of Microbial Communities
#
# This script is the first of three that performs basic exploratory data analysis
# for the SABR data set.
#
# Section 1. Exploration of dataset: Calculates summary statistic, read counts &
# Good's coverage,
#
# Section 2. Data Transformation & Export: Transforms the raw count data to
# relative abundance and performs rarefaction of sequences,
#
# Section 3. ASV Prevalence Analysis: Explores prevalence and core taxa analysis
#
#
# Author: Jaejin Lee & Bolívar Aponte Rolón
# Last modified: 2025-08-21
##################################################################################

# Setup
source("R/utils/000_setup.R")

#--------------------------------------------------------
# SECTION 1: Exploration of data set
#--------------------------------------------------------
## Basic metadata exploration
colnames(tax_table(sabr_2023_physeq))
colnames(sample_data(sabr_2023_physeq))

#--------------------------------------------------------
# Summaries and Read Counts
#--------------------------------------------------------
metagMisc::phyloseq_summary(sabr_2023_physeq, more_stats = F, long = F)
microbiome::summarize_phyloseq(sabr_2023_physeq)
# All samples have at least 5000 reads

# Taxonomic distribution
percent_phyla_clean <- phyloseq_ntaxa_by_tax(
  sabr_2023_physeq,
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
reads <- readcount(sabr_2023_physeq) |>
  as.data.frame() |>
  rownames_to_column(var = "sample_id") |>
  rename(n_seqs = "readcount(sabr_2023_physeq)") |>
  dplyr::left_join(
    rownames_to_column(sabr_2023_metadata_clean, var = "sample_id"),
    by = "sample_id"
  )

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
cover_chao <- phyloseq_coverage(
  sabr_2023_physeq,
  correct_singletons = T,
  add_attr = T
)
cover_chao <- phyloseq_coverage(
  sabr_2023_physeq,
  correct_singletons = T,
  add_attr = T
)
ps_melt <- sabr_2023_physeq |>
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

#####################################################################
# SECTION 2: Data Transformation and Export
#
# This section of the script transforms the raw count data to
# relative abundance, performs rarefaction of sequences,examines basic
# properties of the dataset, and exports the transformed data for
# further analysis.
#
#####################################################################

# Ensure the phyloseq object (sabr_2023_physeq) is loaded.

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
  file = "data/output/processed/sabr_2023_asv_table_rarefied.rda"
)

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
  file = "data/output/processed/sabr_2023_main_rarefied_physeq.rda"
)

#--------------------------------------------------------
# Calculate Diversity indices
#--------------------------------------------------------
main_rarefied_df <- main_rarefied_df %>%
  mutate(
    observed = rowSums(select(., -c(1:21)) > 0),
    shannon = vegan::diversity(select(., -c(1:21)), index = "shannon"),
    simpson = vegan::diversity(select(., -c(1:21)), index = "simpson"),
    inv_simpson = vegan::diversity(select(., -c(1:21)), index = "invsimpson")
  ) %>%
  relocate(
    any_of(c("observed", "shannon", "simpson", "inv_simpson")),
    .before = ASV_1
  )

save(
  main_rarefied_df,
  file = "data/output/processed/sabr_2023_main_rarefied_df.rda"
)


#####################################################################
# SECTION 3: ASV Prevalence Analysis
#
# This section analyzes the prevalence of ASVs across
# samples, identifying and examining ASVs prevalent in a large
# proportion of samples at different thresholds.
#
# We shift from using an S4 `phyloseq` object to an
# S4 `TreeSummarizedExperiement` object.
#####################################################################

#--------------------------------------------------------------------
# `phyloseq` PACKAGE WORKFLOW
#--------------------------------------------------------------------
#--------------------------------------------------------
# Prevalence & Core Microbiome Analysis
#--------------------------------------------------------

## Prevalence visualization
metagMisc::phyloseq_prevalence_plot(
  sabr_2023_physeq,
  sabr_2023_physeq,
  prev.trh = 0.5,
  taxcolor = "phylum",
  facet = TRUE,
  point_alpha = 0.7,
  showplot = T
)

## Host plant-specific averages
ps_average <- metagMisc::phyloseq_average(
  sabr_2023_physeq,
  sabr_2023_physeq,
  avg_type = "arithmetic",
  acomp_zero_impute = NULL,
  aldex_samples = 213,
  aldex_denom = "all",
  group = "plant",
  drop_group_zero = TRUE,
  verbose = TRUE,
  progress = "text",
)


#--------------------------------------------------------
# Visualize ASV prevalence
#--------------------------------------------------------

# Create a bar plot showing the number of samples each ASV is found in
ggplot(
  sabr_2023_physeq_relab %>%
    phyloseq::otu_table(.) %>%
    t(.) %>%
    as.data.frame(.) %>%
    {
      \(df) colSums(df > 0)
    }() %>%
    {
      \(counts) data.frame(OTU = names(counts), Sample_Counts = counts)
    }(),
  aes(x = reorder(OTU, -Sample_Counts), y = Sample_Counts)
) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Number of ASVs per sample",
    x = "ASV",
    y = "Number of Samples"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

#--------------------------------------------------------
# Identify high-prevalence ASVs (>90% samples)
#--------------------------------------------------------
# Back to using the raw count objects
## Core & rare microbiome

# The core taxa are defined as those that exceed the given population
# prevalence threshold at the given detection level.
microbiome::core_abundance(
  sabr_2023_physeq@otu_table,
  detection = 0,
  prevalence = 50 / 100,
  include.lowest = FALSE
)

# The rarity function provides the abundance of the least abundant taxa
# within each sample, regardless of the population prevalence.

microbiome::rare_abundance(
  sabr_2023_physeq@otu_table,
  detection = 0,
  prevalence = 50 / 100,
  include.lowest = FALSE
)

core_tax <- microbiome::core_members(
  sabr_2023_physeq,
  sabr_2023_physeq,
  detection = 0,
  prevalence = 90 / 100,
  include.lowest = FALSE
)

core_tax

core_tax

rare_tax <- microbiome::rare_members(
  sabr_2023_physeq,
  sabr_2023_physeq,
  detection = 0,
  prevalence = 50 / 100,
  include.lowest = FALSE
)

rare_tax


#----------------------------------------------------------------------
# Extract taxonomy for prevalent and rare ASVs at different thresholds
# and ranks
#----------------------------------------------------------------------

ps_list_og <- get_prevalent_rare(
  sabr_2023_physeq,
  thresholds = c(90, 80, 70, 60, 30),
  detection = 0 / 100,
  include.lowest = FALSE
)

ps_list_relab <- get_prevalent_rare(
  sabr_2023_physeq_relab,
  thresholds = c(90, 80, 70, 60, 30),
  detection = 0 / 100,
  include.lowest = FALSE
)

ps_list_rrfy <- get_prevalent_rare(
  mtr_physeq,
  thresholds = c(90, 80, 70, 60, 30),
  detection = 0 / 100,
  include.lowest = FALSE
)

mtr_physeq_list <- list(
  original_phyloseq_lists = ps_list_og,
  relab_phyloseq_lists = ps_list_relab,
  rarefied_phyloseq_lists = ps_list_rrfy
)

save(
  mtr_physeq_list,
  file = "data/output/processed/sabr_2023_mtr_physeq_list.rda"
)

#----------------------------------------------------------------------
# SECTION 4: mean, SD, and autocorrelations (ACF) for each feature
#----------------------------------------------------------------------

# Inputs
features <- c("gnha", "nitrate_ppm", "ammonia_ppm", "n_available", "gwc_g_g")
parameters <- c("plot", "plant", "sampling_location")
time_col <- "sampling_date"


feature_summary <- map_dfr(features, function(f) {
  set.seed(123)
  x <- suppressWarnings(as.numeric(main_rarefied_df[[f]]))
  x <- x[is.finite(x) & !is.na(x)]
  n_non_na <- length(x)

  # Safe initialization
  sh_W <- NA_real_
  sh_p <- NA_real_
  sh_n <- NA_integer_
  sh_note <- NA_character_

  if (n_non_na >= 3) {
    # Since Missing values are allowed, but the number of non-missing values must be between 3 and 5000.

    x_use <- if (n_non_na > 5000) {
      sh_note <- "Random sample of 5K"
      sample(x, 5000)
    }
    if (n_non_na < 5000) {
      x
    }
    # Extracting useful stats
    sh <- stats::shapiro.test(x_use)
    sh_W <- unname(sh$statistic)
    sh_p <- sh$p.value
    sh_n <- length(x_use)
  }
  if (n_non_na <= 3) {
    sh_note <- "Insufficient non-NA values for Shapiro-Wilk (n < 3)"
  }

  # QQ plots per feature
  qq_plot <-
    if (n_non_na > 3) {
      df_plot <- tibble(value = x)
      ggplot(df_plot, aes(sample = value)) +
        stat_qq(alpha = 0.6) +
        stat_qq_line(
          color = "red",
          linewidth = 0.6,
          distribution = stats::qnorm,
        ) +
        labs(
          title = paste("QQ plot:", str_to_title(gsub("_", " ", f))),
          x = "Theoretical Quantiles",
          y = "Sample Quantiles"
        ) +
        theme_minimal()
    }
  if (n_non_na < 3) {
    ggplot() +
      annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = paste0(f, ": insufficient data (n < 3)")
      ) +
      theme_void()
  }

  tibble(
    feature = f,
    n = n_non_na,
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    min = suppressWarnings(min(x, na.rm = TRUE)),
    max = suppressWarnings(max(x, na.rm = TRUE)),
    shapiro_W = sh_W,
    shapiro_p = sh_p,
    shapiro_n = sh_n,
    shapiro_note = sh_note,
    qqplot = list(qq_plot)
  )
})


feature_summary
feature_summary$qqplot
# Looks like there all is non-normally distributed.

# 2) Global ACF (one series per feature, aggregated to 1 value per date)
#-----------------
# Helpers
#-----------------
acf_preprocess <- function(
  df,
  feature,
  time_col = "sampling_date",
  agg_fun = mean,
  na_rm = TRUE
) {
  series <- df %>%
    filter(!is.na(.data[[time_col]])) %>%
    group_by(.data[[time_col]]) %>%
    summarise(
      value = agg_fun(
        as.numeric(.data[[feature]]),
        na.rm = na_rm
      ),
      .groups = "drop"
    ) %>%
    arrange(.data[[time_col]]) %>%
    filter(is.finite(value))

  return(series)
}


acf_plots <- function(
  df,
  feature,
  time_col = "sampling_date",
  title = NULL,
  lag_max = NULL,
  agg_fun = mean
) {
  series <- acf_preprocess(
    df,
    feature = feature,
    time_col = time_col,
    agg_fun = agg_fun
  )
  ac <- stats::acf(series$value, plot = FALSE, na.action = stats::na.pass)

  df_ac <- data.frame(
    lag = as.integer(ac$lag[, 1, 1]),
    acf = as.numeric(ac$acf[, 1, 1])
  )

  conf <- 1.96 / sqrt(ac$n.used)
  if (is.null(title)) {
    title <- "ACF"
  }

  ggplot(df_ac, aes(lag, acf)) +
    geom_hline(yintercept = 0, color = "grey50") +
    geom_hline(
      yintercept = c(-conf, conf),
      linetype = "dashed",
      color = "red"
    ) +
    geom_col(width = 0.8, fill = "steelblue") +
    labs(title = title, x = "Lag", y = "ACF") +
    theme_minimal()
}
#----------------

# Compute ACF stats
acf_global <- map_dfr(features, function(feat) {
  series <- acf_preprocess(
    main_rarefied_df,
    feature = feat,
    time_col = "sampling_date"
  )
  n <- nrow(series)

  if (n < 3) {
    return(tibble(
      feature = feat,
      lag = NA_integer_,
      acf = NA_real_,
      conf = NA_real_,
      n = n
    ))
  }

  ac <- stats::acf(
    series$value,
    type = "correlation",
    plot = FALSE,
    na.action = na.pass
  )

  tibble(
    feature = feat,
    lag = as.integer(ac$lag[, 1, 1]),
    acf = as.numeric(ac$acf[, 1, 1]),
    conf = 1.96 / sqrt(ac$n.used),
    n = ac$n.used
  )
})


acf_global <- acf_global %>%
  mutate(
    plots = map(
      feature,
      ~ acf_plots(
        df = main_rarefied_df,
        feature = .x,
        time_col = "sampling_date",
        title = paste("ACF:", str_to_title(gsub("_", " ", .x)))
      )
    )
  )
acf_global$plots[[1]]
acf_global$plots[[20]]

# # 3) Grouped ACF (by each parameter), same aggregation to 1 value per date
# acf_by_group <- map_dfr(parameters, function(par) {
#   if (!par %in% names(df0)) {
#     return(tibble())
#   }

#   df0 %>%
#     group_by(.data[[par]]) %>%
#     group_modify(
#       ~ {
#         grp <- .x %>% arrange(date_tmp)
#         map_dfr(features, function(f) {
#           series <- grp %>%
#             group_by(date_tmp) %>%
#             summarise(
#               value = mean(.data[[f]], na.rm = TRUE),
#               .groups = "drop"
#             ) %>%
#             arrange(date_tmp) %>%
#             filter(!is.na(value))

#           n <- nrow(series)
#           if (n < 3) {
#             return(tibble(
#               feature = f,
#               lag = NA_integer_,
#               acf = NA_real_,
#               conf = NA_real_,
#               n = n
#             ))
#           }
#           ac <- stats::acf(series$value, plot = FALSE, na.action = na.pass)
#           tibble(
#             feature = f,
#             lag = as.integer(ac$lag[, 1, 1]),
#             acf = as.numeric(ac$acf[, 1, 1]),
#             conf = 1.96 / sqrt(n),
#             n = n
#           )
#         })
#       }
#     ) %>%
#     ungroup() %>%
#     rename(group = !!par) %>%
#     mutate(parameter = par)
# })

# Outputs:
# - feature_summary: one row per feature with n/mean/sd
# - acf_global: ACF per feature (lags, acf, and ±conf bounds)
# - acf_by_group: ACF per feature within each group and parameter
