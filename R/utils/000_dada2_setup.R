# Setup

# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of CRAN packages to install
cran_pkgs <- c("usethis", "devtools", "Rcpp")

# List of Bioconductor packages to install
bioc_pkgs <- c("Rsamtools", "dada2", "ShortRead", "Biostrings")

# Install missing CRAN packages
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install missing Bioconductor packages
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}
# Load packages
library("usethis")
library("devtools")
library("Rcpp")
library("dada2")
library("ShortRead")
library("Biostrings") # This will install other packages and dependencies.
