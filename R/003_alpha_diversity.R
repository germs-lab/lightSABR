#####################################################################
# Community Diversity Analyses
#
# This script explores general community diversity patterns in SABR data set.
# We explore, richness, alpha and beta diversity using NMDS, PCoA, PERMANOVAS , etc.
#
# We briefly explore coomunity diversity with the phyloseq object but move on with
# the TreeSummarized Experiment object

# Author: Jaejin Lee & Bolívar Aponte Rolón
# Last modified: 2025-08-21
#####################################################################

# Load required libraries
source("R/utils/000_setup.R")


#------------------------------------------------------
# Community Diversity: Alpha Diversity
#--------------------------------------------------------
#----------------------
# Features to test
#----------------------
feat <- c(
  "year",
  "plot",
  "sampling_location",
  "replicate",
  "sampling_date",
  "fertilization"
)
#----------------------

# Temporary/visualization datasets
temp_alpha <- estimate_richness(
  sabr_2023_physeq,
  measures = c("Observed", "Shannon", "Simpson", "InvSimpson")
) %>%
  janitor::clean_names() %>%
  rownames_to_column(var = "sample_id") %>%
  mutate(
    # estimate_richness() messes up sample names, now we need to clean them
    sample_id = gsub(
      pattern = "\\.",
      replacement = "-",
      x = sample_id,
    )
  ) %>%
  dplyr::left_join(
    rownames_to_column(sabr_2023_metadata_clean, var = "sample_id"),
    by = "sample_id"
  ) %>%
  mutate(sampling_date = as.character(sampling_date)) %>%
  tidyr::pivot_longer(
    cols = c(observed, shannon, simpson, inv_simpson),
    names_to = "Diversity_Index",
    values_to = "Value"
  ) %>%
  select(-c(nitrogen_conc, dna_conc_ng_ul, x16s_copies_g, amo_a_337, cnor_b))


temp_mtr_long <- mtr_rrfy_df %>%
  mutate(sampling_date = as.character(sampling_date)) %>%
  tidyr::pivot_longer(
    cols = c(observed, shannon, simpson, inv_simpson),
    names_to = "Diversity_Index",
    values_to = "Value"
  )


# Master list with Alpha diversity plots for raw and rarefied data
# Alpha diversity measures for relative abundance dataset are the same as raw dataset.
# Define the datasets and their corresponding names
datasets <- list(
  raw = temp_alpha,
  rarefied = temp_mtr_long
)

alpha_values <- list(
  raw = 0.8,
  rarefied = 0.03
)

# Create the master list using nested purrr::map()
alpha_plots_master <- purrr::map(
  names(datasets),
  ~ {
    dataset_name <- .x
    dataset <- datasets[[.x]]
    alpha_val <- alpha_values[[.x]]

    feature_plots <- purrr::map(
      feat,
      ~ {
        alpha_plots(
          dataset,
          feat = !!rlang::sym(.x),
          .x = !!rlang::sym(.x),
          .y = Value,
          .color = !!rlang::sym(.x),
          .facet = Diversity_Index,
          title = paste("Alpha Diversity by", .x, "-", toupper(dataset_name)),
          subtitle = "Faceted by Diversity Index",
          x_lab = .x,
          y_lab = "Alpha Diversity Measure",
          legend_by = .x,
          alpha = alpha_val
        )
      }
    )

    names(feature_plots) <- feat
    return(feature_plots)
  }
)

# Name the top-level list
names(alpha_plots_master) <- names(datasets)


# Fixing `sampling_date`
alpha_plots_master$raw$sampling_date <- alpha_plots_master$raw$sampling_date +
  theme(axis.text.x = element_blank())
alpha_plots_master$rarefied$sampling_date <- alpha_plots_master$rarefied$sampling_date +
  theme(axis.text.x = element_blank())

# Print
purrr::iwalk(
  alpha_plots_master,
  ~ {
    dataset_name <- .y
    cat("\n", rep("=", 60), "\n")
    cat("DATASET:", toupper(dataset_name), "\n")
    cat(rep("=", 60), "\n\n")

    purrr::iwalk(
      .x,
      ~ {
        feature_name <- .y
        cat(rep("-", 40), "\n")
        cat("FEATURE:", toupper(feature_name), "\n")
        cat(rep("-", 40), "\n")
        print(.x)
        cat("\n\n")
      }
    )
  }
)
