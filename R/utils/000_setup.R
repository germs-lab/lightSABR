#####################################################################
# Common Functions and Utilities for Microbial Community Analysis
#
# This file contains shared functions, constants, and utilities that
# are used across multiple analysis scripts.
#
# Author: Bolívar Aponte Rolón
# Date: 2025-05-05
#####################################################################

# Package and Environment setup

invisible(
  lapply(
    c(
      "conflicted",
      "phyloseq",
      "vegan",
      "tidyverse",
      "janitor",
      "mia",
      "microbiome",
      "metagMisc",
      "BRCore"
    ),
    library,
    character.only = TRUE
  )
)

# List files and source each
list.files(here::here("R/functions"), pattern = "\\.R$", full.names = TRUE) %>%
  purrr::map(source)

# Objects
list.files(
  here::here("data/output/processed"),
  full.names = TRUE,
  recursive = TRUE,
  pattern = "\\.rda$"
) %>%
  purrr::walk(~ load(.x, envir = .GlobalEnv))


# Solve known conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("intersect", "base")
conflict_prefer("survival", "cluster")
conflict_prefer("setdiff", "Biostrings")
