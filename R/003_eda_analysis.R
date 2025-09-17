##################################################################################
# Exploratory Data Analysis of Microbial Communities
#
# This script is the first of three that performs basic exploratory data analysis
# for the SABR data set.
#
# Section 1. Exploration of dataset: Calculates summary statistic, read counts &
# Good's coverage,
#
# Section 2. ASV Prevalence Analysis: Explores prevalence and core taxa analysis
#
# Section 3: Summary statistics, autocorrelations (ACF) and correlation coefficients
#  for each feature

# Author: Jaejin Lee & Bolívar Aponte Rolón
# Last modified: 2025-08-21
##################################################################################

# Setup
source("R/utils/000_setup.R")

#--------------------------------------------------------
# SECTION 1: Exploration of data set
#--------------------------------------------------------
## Basic metadata exploration
base::colnames(tax_table(sabr_2023_physeq))
base::colnames(sample_data(sabr_2023_physeq))
names_list <- base::rownames(sample_data(sabr_2023_physeq))

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

percent_phyla_clean

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

# reads |>
#   arrange(n_seqs) |>
#   filter(sample_id %in% names_list)

#--------------------------------------------------------
# Good's Coverage Analysis
#--------------------------------------------------------
## Coverage estimates
cover_chao <- phyloseq_coverage(
  sabr_2023_physeq,
  correct_singletons = T,
  add_attr = T
)

cover_chao

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
# SECTION 2: ASV Prevalence Analysis
#
# This section analyzes the prevalence of ASVs across
# samples, identifying and examining ASVs prevalent in a large
# proportion of samples at different thresholds.
#
# We shift from using an S4 `phyloseq` object to an
# S4 `TreeSummarizedExperiement` object.
#####################################################################

#--------------------------------------------------------
# Prevalence & Core Microbiome Analysis
#--------------------------------------------------------

## Prevalence visualization
metagMisc::phyloseq_prevalence_plot(
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
  avg_type = "arithmetic",
  acomp_zero_impute = NULL,
  aldex_samples = 213,
  aldex_denom = "all",
  group = "plant",
  drop_group_zero = TRUE,
  verbose = TRUE,
  progress = "text",
)
ps_average


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
core_abun <- microbiome::core_abundance(
  sabr_2023_physeq@otu_table,
  detection = 0,
  prevalence = 50 / 100,
  include.lowest = FALSE
)

# The rarity function provides the abundance of the least abundant taxa
# within each sample, regardless of the population prevalence.

rare_abun <- microbiome::rare_abundance(
  sabr_2023_physeq@otu_table,
  detection = 0,
  prevalence = 50 / 100,
  include.lowest = FALSE
)

core_tax <- microbiome::core_members(
  sabr_2023_physeq,
  detection = 0,
  prevalence = 90 / 100,
  include.lowest = FALSE
)

core_tax

rare_tax <- microbiome::rare_members(
  sabr_2023_physeq,
  detection = 0,
  prevalence = 50 / 100,
  include.lowest = FALSE
)

#rare_tax

#----------------------------------------------------------------------
# Extract taxonomy for prevalent and rare ASVs at different thresholds
# and ranks
#----------------------------------------------------------------------

ps_list_raw <- get_prevalent_rare(
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

ps_list_rarefied <- get_prevalent_rare(
  main_rarefied_physeq,
  thresholds = c(90, 80, 70, 60, 30),
  detection = 0 / 100,
  include.lowest = FALSE
)

main_physeq_list <- list(
  raw_phyloseq_lists = ps_list_raw,
  relab_phyloseq_lists = ps_list_relab,
  rarefied_phyloseq_lists = ps_list_rarefied
)
ps_list_raw

# save(
#   main_physeq_list,
#   file = "data/output/processed/rdata/phyloseq/sabr_2023_main_physeq_list.rda"
# )

#----------------------------------------------------------------------
# SECTION 3: mean, SD, and autocorrelations (ACF) for each feature
#----------------------------------------------------------------------

# Inputs
features <- c(
  "gnha",
  "nitrate_ppm",
  "ammonia_ppm",
  "amo_a_337",
  "cnor_b",
  #  "n_available",
  "gwc_g_g"
)
parameters <- c("plot", "plant", "sampling_location")
time_col <- "sampling_date"
use_log_y <- TRUE
#-------------------

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
# Use log10 scale, or other normalization technique, when plotting or testing in models

# LEt's look at some boxplots

feature_boxplots <- features %>%
  set_names() %>%
  map(function(feat) {
    parameters %>%
      set_names() %>%
      map(function(par) {
        df <- sabr_2023_metadata_clean %>%
          {
            .[(rownames(.) %in% main_rarefied_df$sample_id), , drop = FALSE]
          } %>%
          mutate(
            !!feat := as.numeric(.data[[feat]]),
            !!par := as.factor(.data[[par]])
          )

        df <- if (use_log_y) {
          df %>% filter(!is.na(.data[[feat]]), .data[[feat]] > 0)
        } else {
          df %>% filter(!is.na(.data[[feat]]))
        }

        p <- ggplot(
          df,
          aes(x = .data[[par]], y = .data[[feat]], color = .data[[par]])
        ) +
          geom_boxplot(
            #,
            outlier.shape = NA,
            alpha = 0.35
          ) +
          geom_jitter(
            aes(color = .data[[par]]),
            width = 0.2,
            height = 0,
            alpha = 0.25
          ) +
          labs(
            title = paste(
              str_to_title(gsub("_", " ", feat)),
              "by",
              str_to_title(gsub("_", " ", par))
            ),
            x = str_to_title(gsub("_", " ", par)),
            y = paste("log10(", str_to_title(gsub("_", " ", feat)), ")"),
            color = str_to_title(gsub("_", " ", par))
          ) +
          scale_color_manual(values = pals::glasbey(32)) +
          #scale_fill_manual(values = pals::glasbey(32))
          theme_minimal() +
          theme(legend.position = "right")

        if (use_log_y) {
          p <- p + scale_y_log10(labels = scales::label_number())
        }

        p
      })
  })


purrr::iwalk(
  feature_boxplots,
  ~ {
    feature_name <- .y
    cat(rep("-", 40), "\n")
    cat("FEATURE:", toupper(feature_name), "\n")
    cat(rep("-", 40), "\n")
    print(.x)
    cat("\n\n")
  }
)


# 2) Global ACF (one series per feature, aggregated to 1 value per date)
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
  distinct(feature) %>% # one row per feature
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
acf_global$plots


# 3) Grouped ACF (by each parameter), same aggregation to 1 value per date
acf_by_group <- map_dfr(parameters, function(par) {
  if (!par %in% names(main_rarefied_df)) {
    return(tibble())
  }

  main_rarefied_df %>%
    group_by(.data[[par]]) %>%
    group_modify(
      ~ {
        map_dfr(features, function(f) {
          series <- acf_preprocess(
            main_rarefied_df,
            feature = f,
            time_col = "sampling_date"
          )

          n <- nrow(series)
          if (n < 3) {
            return(tibble(
              feature = f,
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
            feature = f,
            lag = as.integer(ac$lag[, 1, 1]),
            acf = as.numeric(ac$acf[, 1, 1]),
            conf = 1.96 / sqrt(n),
            n = n
          )
        })
      }
    ) %>%
    ungroup() %>%
    rename(group = !!par) %>%
    mutate(parameter = par)
})


# Correlation coefficients

metadata_matrix <- sabr_2023_metadata_clean |>
  select(c(gnha, amo_a_337:gwc_g_g))

cor(
  metadata_matrix,
  use = "complete.obs",
  method = 'pearson'
)


ggcorrplot::ggcorrplot(
  cor(
    metadata_matrix,
    use = "complete.obs",
    method = 'pearson'
  ),
  method = "square",
  hc.order = TRUE,
  lab = TRUE,
  digits = 3,
  type = "lower",
  p.mat = ggcorrplot::cor_pmat(metadata_matrix),
  insig = "blank",
  colors = c("#6D9EC1", "white", "#E46726")
)

# To be expected, nitrate_ppm, ammonia_ppm and N_available are highly correlated.
# Specially nitrate_ppm and n_available. For downstream analysses we will focus
# ammonia_ppm and nitrate_ppm.
