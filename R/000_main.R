#####################################################################
# Main Script for Microbial Community Analysis
#
# This script coordinates the execution of all analysis modules
# for the microbial community analysis project.
#
# Author: Bolívar Aponte Rolón
# Date: 2025-05-05
#####################################################################

# Set working directory if needed
# setwd("/path/to/your/project")

# Load required libraries
source("R/utils/000_setup.R")
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(UpSetR)
library(reshape2)

# Create output directories if they don't exist
dir.create("outputs", showWarnings = FALSE)
dir.create("analysis", showWarnings = FALSE)

# Load your phyloseq object
# This should be prepared in a separate script not included in this analysis
# For example:
# source("scripts/00_data_preparation.R")
# You should have a phyloseq object named 'ps' loaded in your environment

# If you need to load from a saved RDS file
# ps <- readRDS("data/processed/phyloseq_object.rds")

# Check that the phyloseq object is properly loaded
if (!exists("ps")) {
  stop("Phyloseq object 'ps' not found. Please load your data first.")
}

# Run the individual analysis scripts
cat("Starting NMDS analysis...\n")
source("01_nmds_analysis.R")

cat("Transforming data and exporting...\n")
source("02_data_transformation.R")

cat("Analyzing ASV prevalence...\n")
source("03_asv_prevalence.R")

cat("Analyzing core microbiome...\n")
source("04_core_microbiome.R")

cat("Analyzing prevalence vs. abundance...\n")
source("05_prevalence_abundance.R")

cat("Performing comparative analysis...\n")
source("06_comparative_analysis.R")

cat(
  "Analysis complete! Results have been saved to the 'outputs' and 'analysis' directories.\n"
)

# Optional: Save the workspace for future reference
# save.image("analysis/workspace.RData")
