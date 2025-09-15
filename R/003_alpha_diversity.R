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


temp_main_long <- main_rarefied_df %>%
  mutate(sampling_date = as.character(sampling_date)) %>%
  tidyr::pivot_longer(
    cols = c(observed, shannon, simpson, inv_simpson),
    names_to = "Diversity_Index",
    values_to = "Value"
  )


# Main list with Alpha diversity plots for raw and rarefied data
# Alpha diversity measures for relative abundance dataset are the same as raw dataset.
# Define the datasets and their corresponding names
#----------------------
# Features to test
#----------------------
cat_feat <- c(
  "year",
  "plot",
  "sampling_location",
  "replicate",
  "sampling_date",
  "fertilization"
)
cont_feat <- c(
  "gnha",
  "nitrate_ppm",
  "ammonia_ppm",
  "n_available",
  "gwc_g_g"
)
#----------------------

datasets <- list(
  # raw = temp_alpha,
  rarefied = temp_main_long
)

alpha_values <- list(
  # raw = 0.8,
  rarefied = 0.03
)

# Create the main list using nested purrr::map()
alpha_cat_plots_main <- purrr::map(
  names(datasets),
  ~ {
    dataset_name <- .x
    dataset <- datasets[[.x]]
    alpha_val <- alpha_values[[.x]]

    feature_plots <- purrr::map(
      cat_feat,
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

    names(feature_plots) <- cat_feat
    return(feature_plots)
  }
)

# Name the top-level list
names(alpha_cat_plots_main) <- names(datasets)


# Fixing `sampling_date`
# alpha_cat_plots_main$raw$sampling_date <- alpha_cat_plots_main$raw$sampling_date +
#   theme(axis.text.x = element_blank())
alpha_cat_plots_main$rarefied$sampling_date <- alpha_cat_plots_main$rarefied$sampling_date +
  theme(axis.text.x = element_blank())

# Print
purrr::iwalk(
  alpha_cat_plots_main,
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

# Auto correlation between feastures/covariates Correlations
# How does the data look like for features?

main_rarefied_df |>
  filter(!is.na(nitrate_ppm), !is.na(gwc_g_g)) |>
  ggplot(aes(x = nitrate_ppm, y = observed)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.03)
scale_x_log10() +
  scale_y_log10()


parameters <- c("plot", "plant", "sampling_location")

cont_feat <- c(
  "gnha",
  "nitrate_ppm",
  "ammonia_ppm",
  "n_available",
  "gwc_g_g"
)

use_log_y <- TRUE

alpha_param_plots <- cont_feat %>%
  set_names() %>%
  map(function(feat) {
    parameters %>%
      set_names() %>%
      map(function(par) {
        df <- main_rarefied_df %>%
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
          geom_boxplot(outlier.shape = NA, alpha = 0.35) +
          geom_jitter(width = 0.2, height = 0, alpha = 0.25) +
          labs(
            title = paste(
              str_to_title(gsub("_", " ", feat)),
              "by",
              str_to_title(gsub("_", " ", par))
            ),
            x = str_to_title(gsub("_", " ", par)),
            y = str_to_title(gsub("_", " ", feat)),
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

alpha_param_plots$gnha$plot
alpha_param_plots$gnha$sampling_location

# or
# alpha_nitrate_plots[["shannon"]]
# alpha_nitrate_plots[["inv_simpson"]]
