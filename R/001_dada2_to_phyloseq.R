#####################################################################
# DADA2 to Phyloseq Workflow for Amplicon Sequencing Data
#
# This script processes raw amplicon sequencing data (16S rRNA) through
# the DADA2 pipeline and creates a phyloseq object for downstream analysis.
# The workflow includes quality filtering, error learning, ASV inference,
# chimera removal, and taxonomy assignment.
#
# Author: Jaejin Lee
# Modified by: Bolívar Aponte Rolón
# Date: 2025-05-05
#####################################################################

#--------------------------------------------------------
# Setup and Library Loading
#--------------------------------------------------------

# Load setup scripts
# source("R/utils/000_dada2_setup.R") # Only run the first time to install packages
source("R/utils/000_setup.R")

# Load required libraries
library("usethis")
library("devtools")
library("Rcpp")
library("dada2")
library("ShortRead")
library("Biostrings")

#--------------------------------------------------------
# Step 1: Set up file paths and identify sequence files
#--------------------------------------------------------

# Set working directory for input files
in_path <- "data/input/raw_sequences"
out_dir <- "data/output"

# List all fastq files
all_files <- list.files(in_path, pattern = "\\.fastq\\.gz$", full.names = TRUE)

# Separate forward (R1) and reverse (R2) read files
fnFs <- sort(grep("_R1", all_files, value = TRUE))
fnRs <- sort(grep("_R2", all_files, value = TRUE))
print(fnFs)
print(fnRs)

# Define output paths for filtered files
filtFs <- file.path(out_dir, "filtered", basename(fnFs))
filtRs <- file.path(out_dir, "filtered", basename(fnRs))

#--------------------------------------------------------
# Step 2: Quality assessment and filtering
#--------------------------------------------------------

# Check quality profiles of the first sample
plotQualityProfile(fnFs[1]) # Quality check for first R1 file
plotQualityProfile(fnRs[1]) # Quality check for first R2 file

# Filter and trim reads based on quality parameters
filter_stats <- filterAndTrim(
  fnFs,
  filtFs,
  fnRs,
  filtRs,
  truncLen = c(140, 140), # Truncate reads at position 140
  maxN = 0, # Maximum number of N bases allowed
  maxEE = c(4, 4), # Maximum expected errors allowed
  truncQ = 2, # Truncate reads at first quality score <= 2
  rm.phix = TRUE, # Remove PhiX spike-in sequences
  compress = TRUE, # Compress output files
  multithread = TRUE # Use multiple processors
)
print(filter_stats)

# Verify output files were created
length(filtFs)
length(filtRs)

#--------------------------------------------------------
# Step 3: Learn error rates from the data
#--------------------------------------------------------

# Learn error rates from filtered reads
errF <- learnErrors(filtFs, multithread = TRUE) # Forward reads
errR <- learnErrors(filtRs, multithread = TRUE) # Reverse reads

# Visualize error models to ensure they fit well
# Forward read error model
plotErrors(errF, nominalQ = TRUE) +
  scale_y_continuous(
    trans = "log10",
    limits = c(1e-6, 1),
    breaks = c(1e-1, 1e-3, 1e-5),
    labels = c("-1", "-3", "-5")
  ) +
  labs(y = "Occurrence frequency (log10)")

# Reverse read error model
plotErrors(errR, nominalQ = TRUE) +
  scale_y_continuous(
    trans = "log10",
    limits = c(1e-6, 1),
    breaks = c(1e-1, 1e-3, 1e-5),
    labels = c("-1", "-3", "-5")
  ) +
  labs(y = "Occurrence frequency (log10)")

#--------------------------------------------------------
# Step 4: Dereplicate reads
#--------------------------------------------------------

# Dereplicate forward reads (combine identical sequences)
derepFs <- derepFastq(filtFs, verbose = TRUE)

# Dereplicate reverse reads
derepRs <- derepFastq(filtRs, verbose = TRUE)

# Extract sample names from filenames
sample.names <- gsub("^(.*?)_S\\d+.*$", "\\1", basename(filtFs))
print(sample.names)

# Assign sample names to dereplicated sequences
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#--------------------------------------------------------
# Step 5: Infer ASVs (Amplicon Sequence Variants)
#--------------------------------------------------------

# Denoise forward reads to infer ASVs
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)

# Denoise reverse reads to infer ASVs
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#--------------------------------------------------------
# Step 6: Merge paired-end reads
#--------------------------------------------------------

# Merge forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Summarize merging results
merger_summary <- sapply(mergers, function(x) sum(x$abundance))
print(merger_summary)

#--------------------------------------------------------
# Step 7: Evaluate merging efficiency
#--------------------------------------------------------

# Number of reads after filtering
reads_in <- sapply(dadaFs, function(x) sum(x$denoised))

# Number of merged reads
reads_merged <- sapply(mergers, function(x) sum(x$abundance))

# Calculate merging rate (percentage of reads that were successfully merged)
merge_rate <- reads_merged / reads_in * 100
merge_rate

#--------------------------------------------------------
# Step 8: Filter low-coverage samples
#--------------------------------------------------------

# Identify samples below different read count thresholds
low_samples_5000 <- names(merger_summary[merger_summary < 5000])
low_samples_10000 <- names(merger_summary[merger_summary < 10000])

# Identify samples below 50% of mean sequencing depth
mean_seq <- mean(merger_summary)
low_samples_mean <- names(merger_summary[merger_summary < (mean_seq * 0.5)])

# Print threshold results
cat("Sample names <5000 seqs:\n", low_samples_5000, "\n\n")
cat("#Samples <5000 seqs:", length(low_samples_5000), "\n")

cat("Sample names <10000 seqs:\n", low_samples_10000, "\n\n")
cat("#Samples <10000 seqs:", length(low_samples_10000), "\n")

cat(
  "Samples names - below 50% of the mean sequencing depth:\n",
  low_samples_mean,
  "\n\n"
)
cat(
  "#Samples - below 50% of the mean sequencing depth:",
  length(low_samples_mean),
  "\n"
)

# Visualize sequencing depth distribution with threshold lines
depth_data <- data.frame(sequencing_depth = merger_summary)

ggplot(depth_data, aes(x = sequencing_depth)) +
  geom_histogram(
    aes(y = ..density..),
    binwidth = 1000,
    fill = "steelblue",
    alpha = 0.6
  ) +
  geom_density(color = "darkred", size = 1) +
  geom_vline(
    xintercept = c(5000, 10000, mean(merger_summary) * 0.5),
    color = "red",
    linetype = "dashed",
    size = 0.8
  ) +
  geom_vline(
    xintercept = mean(merger_summary),
    color = "black",
    linetype = "dashed",
    size = 0.8
  ) +
  labs(
    title = "Sequencing Depth Distribution",
    x = "Sequencing Depth",
    y = "Density"
  ) +
  theme_minimal()

# Remove samples with fewer than 5000 sequences
mergers_filtered <- mergers[!names(mergers) %in% low_samples_5000]
cat("#Samples after exclusion:", length(mergers_filtered), "\n")

#--------------------------------------------------------
# Step 9: Create ASV table and filter by length
#--------------------------------------------------------

# Create initial ASV table
seqtab <- makeSequenceTable(mergers_filtered)

# Check sequence length distribution
table(nchar(getSequences(seqtab)))

# Visualize sequence length distribution
seq_lengths <- nchar(getSequences(seqtab))
seq_length_table <- table(seq_lengths)

seq_length_df <- as.data.frame(seq_length_table)
colnames(seq_length_df) <- c("Length", "Count")

ggplot(seq_length_df, aes(x = factor(Length), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = 20, color = "red", linetype = "dotted", size = 1) +
  scale_y_continuous(
    trans = "log10",
    breaks = c(1, 10, 100, 1000, 10000),
    labels = c("0", "1", "2", "3", "4")
  ) +
  labs(
    title = "Sequence Length Distribution",
    x = "Sequence Length (bp)",
    y = "Log10(Frequency)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Filter sequences to keep only those within the expected length range
expected_range <- seq_lengths >= 251 & seq_lengths <= 257
seqtab_filtered <- seqtab[, expected_range]

# Compare sequence counts before and after length filtering
cat("Number of sequences before filtering:", ncol(seqtab), "\n")
cat("Number of sequences after filtering:", ncol(seqtab_filtered), "\n")

#--------------------------------------------------------
# Step 10: Remove chimeric sequences
#--------------------------------------------------------

# Identify and remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(
  seqtab_filtered,
  method = "consensus", # Use consensus method to identify chimeras
  multithread = TRUE, # Use multiple processors
  verbose = TRUE # Print progress information
)

# Check sequence length distribution after chimera removal
table(nchar(getSequences(seqtab.nochim)))

# Compare sequence counts before and after chimera removal
cat("Before removing chimeras:", sum(seqtab_filtered), "\n")
cat("After removing chimeras:", sum(seqtab.nochim), "\n")

# Save file
save(
  seqtab.nochim,
  file = file.path(out_dir, "processed/sabr_2023_asv_table.rda")
)
write.csv(
  seqtab.nochim,
  file.path(out_dir, "processed/sabr_2023_asv_table.csv")
) # Long file name but it indicates this file has gone through all the steps in the pipeline.

load(
  file.path(out_dir, "processed/sabr_2023_asv_table.rda")
)

#--------------------------------------------------------
# Step 11: Assign taxonomy using SILVA database
#--------------------------------------------------------

# Assign taxonomy using the SILVA reference database
taxa <- assignTaxonomy(
  seqtab.nochim,
  "data/input/SILVA/silva_nr99_v138.2_toSpecies_trainset.fa.gz",
  multithread = TRUE
)

# Replace NA values with "Unclassified"
taxa[is.na(taxa)] <- "Unclassified"

# Examine taxonomy assignments
head(taxa)
str(taxa)
dim(taxa)
colnames(taxa)
summary(taxa)

# Save object
write.csv(
  taxa,
  file.path(out_dir, "processed/sabr_2023_taxonomy.csv")
)

taxa <- read.csv(file.path(out_dir, "processed/sabr_2023_taxonomy.csv")) %>%
  rename(., sequence = X) %>%
  rename_with(str_to_lower, .cols = everything()) # Clean up needed after importing from .csv

rownames(taxa) <- paste0("ASV_", 1:nrow(taxa))

#--------------------------------------------------------
# Step 12: Create phyloseq object
#--------------------------------------------------------

# Create ASV table
## Cleaning ASV names for FASTA file

asv_fasta <- seqtab2fasta(seqtab.nochim)

seqtab.nochim <- t(seqtab.nochim) # Retaining sequences and asigning shorthand ASV names

row.names(seqtab.nochim) <- sub(">", "", asv_fasta$asv_headers)

## Save FASTA file
write(
  asv_fasta$asv_fasta,
  file.path(out_dir, "processed/sabr_2023_asv.fa")
)

asv <- otu_table(seqtab.nochim, taxa_are_rows = TRUE)

# Create taxonomy table
tax <- tax_table(as.matrix(taxa))

# Create initial phyloseq object
physeq <- phyloseq(asv, tax)
physeq

#--------------------------------------------------------
# Step 13: Add sample metadata
#--------------------------------------------------------

# Load metadata from CSV file
metadata <- readxl::read_xlsx(
  "data/input/metadata_CABBI_SABR2023_DNA.xlsx"
) |>
  janitor::clean_names() |>
  mutate(
    id = gsub(
      pattern = "^(.*?)_R\\d+.*$",
      replacement = "\\1",
      x = id
    ),
    across(where(is.numeric), as.factor)
  ) |>
  distinct(id, .keep_all = TRUE) |>
  column_to_rownames(var = "id")


sabr_2023_metadata_clean <- metadata # Name change for saving purposes
save(
  sabr_2023_metadata_clean,
  file = "data/output/processed/sabr_2023_metadata_clean.rda"
)

# Check for sample name consistency between phyloseq and metadata
head(sample_names(physeq))
head(rownames(metadata))
base::setdiff(sample_names(physeq), rownames(metadata)) # In phyloseq but not in metadata
base::setdiff(rownames(metadata), sample_names(physeq)) # In metadata but not in phyloseq

# Convert metadata to phyloseq sample_data format
sampledata <- sample_data(metadata)
sampledata

# Add metadata to phyloseq object
physeq <- merge_phyloseq(physeq, sampledata)
physeq

# Verify phyloseq components
otu_table(physeq) # ASV table
sample_data(physeq) # Sample metadata
tax_table(physeq) # Taxonomy table

# #--------------------------------------------------------
# # Step 14: Save results
# #--------------------------------------------------------

# # Save phyloseq object as RDS file
save(
  physeq,
  file = "data/output/processed/sabr_2023_physeq_object.rda"
)
