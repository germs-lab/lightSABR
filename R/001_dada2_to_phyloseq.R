# Setup
# source("R/utils/000_dada2_setup.R") # Only run the first time to install packages
source("R/utils/000_setup.R")

library("usethis")
library("devtools")
library("Rcpp")
library("dada2")
library("ShortRead")
library("Biostrings") 



# set working directory
path <- "data/input/raw_sequences"

# list all files
all_files <- list.files(path, pattern = "\\.fastq\\.gz$", full.names = TRUE)

# separate R1 and R2 files
fnFs <- sort(grep("_R1", all_files, value = TRUE))
fnRs <- sort(grep("_R2", all_files, value = TRUE))
print(fnFs)
print(fnRs)

# set path for filtered files
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))

plotQualityProfile(fnFs[1]) # quality check for first R1 file
plotQualityProfile(fnRs[1]) # quality check for first R2 file

# filtering and trimming
filter_stats <- filterAndTrim(
  fnFs,
  filtFs,
  fnRs,
  filtRs,
  truncLen = c(140, 140),
  maxN = 0,
  maxEE = c(4, 4),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)
print(filter_stats)

length(filtFs)
length(filtRs)

# learning errors
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# plot the original error model
plotErrors(errF, nominalQ = TRUE) +
  scale_y_continuous(
    trans = "log10",
    limits = c(1e-6, 1),
    breaks = c(1e-1, 1e-3, 1e-5),
    labels = c("-1", "-3", "-5")
  ) +
  labs(y = "Occurrence frequency (log10)")
plotErrors(errR, nominalQ = TRUE) +
  scale_y_continuous(
    trans = "log10",
    limits = c(1e-6, 1),
    breaks = c(1e-1, 1e-3, 1e-5),
    labels = c("-1", "-3", "-5")
  ) +
  labs(y = "Occurrence frequency (log10)")

# forward reads dereplication
derepFs <- derepFastq(filtFs, verbose = TRUE)

# reverse reads dereplication
derepRs <- derepFastq(filtRs, verbose = TRUE)

# remove R1/R2 and file extension
sample.names <- gsub("_R[12]\\.fastq\\.gz", "", basename(filtFs))
print(sample.names)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Denoising
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# merging forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# summarizing merging result
merger_summary <- sapply(mergers, function(x) sum(x$abundance))
print(merger_summary)

---
  # number of reads after filtering
  reads_in
<- sapply(dadaFs, function(x) sum(x$denoised))

# number of merged reads
reads_merged <- sapply(mergers, function(x) sum(x$abundance))

# calculating merging rate
merge_rate <- reads_merged / reads_in * 100
merge_rate

---
  # applying an absolute threshold
  low_samples_5000
<- names(merger_summary[merger_summary < 5000])
low_samples_10000 <- names(merger_summary[merger_summary < 10000])

# identifying samples below 50% of the mean sequencing depth
mean_seq <- mean(merger_summary)
low_samples_mean <- names(merger_summary[merger_summary < (mean_seq * 0.5)])

# print threshold results
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

# visualizing merger_summary with potent thresholds
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

# excluding samples with <5000 seqs
mergers_filtered <- mergers[!names(mergers) %in% low_samples_5000]
cat("#Samples after exclusion:", length(mergers_filtered), "\n")

---
  # creating initial ASV table
  seqtab
<- makeSequenceTable(mergers_filtered)

# check the distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# get sequence lengths
seq_lengths <- nchar(getSequences(seqtab))
seq_length_table <- table(seq_lengths)

seq_length_df <- as.data.frame(seq_length_table)
colnames(seq_length_df) <- c("Length", "Count")

ggplot(seq_length_df, aes(x = factor(Length), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = 20, color = "red", linetype = "dotted", size = 1) +
  scale_y_continuous(
    trans = "log10",
    breaks = c(1, 10, 100, 1000, 10000), # 실제 로그 값 설정
    labels = c("0", "1", "2", "3", "4")
  ) + # 로그 값 라벨링
  labs(
    title = "Sequence Length Distribution",
    x = "Sequence Length (bp)",
    y = "Log10(Frequency)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# define the expected length range
expected_range <- seq_lengths >= 251 & seq_lengths <= 257

# Filter sequences within the expected length range
seqtab_filtered <- seqtab[, expected_range]

# Compare the number of sequences before and after filtering
cat("Number of sequences before filtering:", ncol(seqtab), "\n")
cat("Number of sequences after filtering:", ncol(seqtab_filtered), "\n")


# removing chimeras
seqtab.nochim <- removeBimeraDenovo(
  seqtab_filtered,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

# check the distribution of sequence lengths after chimera removal
table(nchar(getSequences(seqtab.nochim)))

# number of sequences before and after removing chimeras
cat("Before removing chimeras:", sum(seqtab_filtered), "\n")
cat("After removing chimeras:", sum(seqtab.nochim), "\n")


# Taxonomy assignment using the SILVA database
taxa <- assignTaxonomy(
  seqtab.nochim,
  "/Users/jaejinlee/Files/Data/databases/silva_nr99_v138.2_toSpecies_trainset.fa.gz",
  multithread = TRUE
)
taxa[is.na(taxa)] <- "Unclassified"

head(taxa)

str(taxa)
dim(taxa)
colnames(taxa)
summary(taxa)

# creating final ASV table
asv <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)

# creating taxonomy table
tax <- tax_table(as.matrix(taxa))

# creating Phyloseq object
ps <- phyloseq(asv, tax)
ps

# loading metadata
metadata <- read.csv(
  "/Users/jaejinlee/Files/Data/2023SABR_amplicon/SABR2023_metadata.csv",
  row.names = 1
)
head(metadata)

head(sample_names(ps)) # check sample names in ps
head(rownames(metadata)) # check sample names in metadata
setdiff(sample_names(ps), rownames(metadata)) # present in in, absent in metadata
setdiff(rownames(metadata), sample_names(ps)) # present in metadata, absent in ps

# converting metadata to sample_data
sampledata <- sample_data(metadata)
sampledata

# adding metadata to ps
ps <- merge_phyloseq(ps, sampledata)
ps

sample_data(ps)

otu_table(ps) # asv table
sample_data(ps) # metadata
tax_table(ps) # taxonomy

saveRDS(
  ps,
  file = "/Users/jaejinlee/Files/Data/2023SABR_amplicon/analysis/ps_object.rds"
)
# when calling the ps, use: ps <- readRDS(file = "/Users/jaejinlee/Files/Data/2023SABR_amplicon/analysis/ps_object.rds")

# converting asv table to data frame
asv_data <- as.data.frame(otu_table(ps))

# save as a csv file
write.csv(
  asv_data,
  file = "/Users/jaejinlee/Files/Data/2023SABR_amplicon/analysis/asv_table.csv",
  row.names = TRUE
)

