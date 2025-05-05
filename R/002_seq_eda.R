----------
  ### NMDS
  # calculating Bray-Curtis distance
  bc_dist
<- phyloseq::distance(ps, method = "bray")

# NMDS analysis (k=2 reduces to 2 dimensions)
nmds <- ordinate(ps, method = "NMDS", distance = "bray", trymax = 100)

# Visualizing NMDS results
# Use the "Date" column in the metadata to create four panels by sampling date
# Convert "Date" to a factor for proper panel separation
sample_data(ps)$Date <- as.factor(sample_data(ps)$Date)
colnames(sample_data(ps))

# Generate NMDS plot
nmds_plot <- plot_ordination(ps, nmds, color = "Plant", shape = "Location") +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~Date, ncol = 2) + # Separate panels by Date (4 panels)
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

# Print the plot
print(nmds_plot)


summary(as.vector(otu_table(ps)))
sample_sums(ps)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
sample_sums(ps_rel)

otu_table(ps_rel)
# converting otu_table(ps_rel) to a data frame
asv_table_rel_df <- as.data.frame(otu_table(ps_rel))

# save as a csv file
write.csv(
  asv_table_rel_df,
  file = "/Users/jaejinlee/Files/Data/2023SABR_amplicon/analysis/asv_table_ps_rel.csv",
  row.names = TRUE
)

dim(otu_table(ps_rel))

# Calculating number of samples containing each ASV
asv_table_data <- as.data.frame(otu_table(ps_rel))
asv_sample_counts <- colSums(asv_table_data > 0)
asv_count_df <- data.frame(
  OTU = names(asv_sample_counts),
  Sample_Counts = asv_sample_counts
)

# visialization
ggplot(asv_count_df, aes(x = reorder(OTU, -Sample_Counts), y = Sample_Counts)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Number of Samples Each OTU is Found In",
    x = "ASV",
    y = "Number of Samples"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), # OTU가 많아 x축 레이블 제거
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# total number of samples
total_samples <- nrow(otu_table_data)
total_samples

# ASV found in >90% of samples
threshold <- total_samples * 0.9
high_prevalence_asvs <- names(asv_sample_counts[asv_sample_counts >= threshold])

cat("ASVs found in >90% samples:", length(high_prevalence_asvs), "\n")
cat("list of ASVs found in >90% samples:\n", high_prevalence_asvs, "\n")

# extracting taxonomy table
tax_table_data <- tax_table(ps_rel)

# taxonomy for ASVs found in >90% samples
high_prevalence_taxa <- tax_table_data[high_prevalence_asvs, ]
print(high_prevalence_taxa)

# a new table for ASVs found in >90% samples
asv_high_prev <- otu_table_data[high_prevalence_asvs, ]
asv_high_prev <- otu_table(asv_high_prev, taxa_are_rows = TRUE)

ps_high_prev <- merge_phyloseq(ps_rel, asv_high_prev)

# NMDS
bc_dist_high_prev <- phyloseq::distance(ps_high_prev, method = "bray")
nmds_high_prev <- ordinate(
  ps_high_prev,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

nmds_plot_high_prev <- plot_ordination(
  ps_high_prev,
  nmds_high_prev,
  color = "Plant",
  shape = "Location"
) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~Date, ncol = 2) +
  labs(
    title = "NMDS of Highly Prevalent OTUs",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_high_prev)


### ASVs found in <70% samples

total_samples <- nrow(otu_table_data)
total_samples

# calculating ASVs found in <70% samples
threshold_u70 <- total_samples * 0.7
low_prevalence_asvs_u70 <- names(asv_sample_counts[
  asv_sample_counts < threshold_u70
])

low_prevalence_taxa_u70 <- tax_table_data[low_prevalence_asvs_u70, ]

cat("ASVs found in <70% samples:", length(low_prevalence_asvs_u70), "\n")
print(low_prevalence_taxa_u70)

# a new table for ASVs found in <70% samples
asv_low_prev_u70 <- otu_table_data[, low_prevalence_asvs_u70]
asv_low_prev_u70 <- otu_table(asv_low_prev_u70, taxa_are_rows = FALSE)

ps_low_prev_u70 <- merge_phyloseq(ps_rel, asv_low_prev_u70)

# NMDS
bc_dist_low_prev_u70 <- phyloseq::distance(ps_low_prev_u70, method = "bray")
nmds_low_prev_u70 <- ordinate(
  ps_low_prev_u70,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

# visualization
nmds_plot_low_prev_u70 <- plot_ordination(
  ps_low_prev_u70,
  nmds_low_prev_u70,
  color = "Plant",
  shape = "Location"
) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~Date, ncol = 2) +
  labs(
    title = "NMDS of Low Prevalence OTUs (Present in <70% of Samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_low_prev_u70)


### ASVs found in <30% samples

threshold_u30 <- total_samples * 0.3
low_prevalence_asvs_u30 <- names(asv_sample_counts[
  asv_sample_counts < threshold_u30
])

low_prevalence_taxa_u30 <- tax_table_data[low_prevalence_asvs_u30, ]

cat("ASVs found in <30% samples:", length(low_prevalence_asvs_u30), "\n")
print(low_prevalence_taxa_u30)

# a new table for ASVs found in <30% samples
asv_low_prev_u30 <- otu_table_data[, low_prevalence_asvs_u30]
asv_low_prev_u30 <- otu_table(asv_low_prev_u30, taxa_are_rows = FALSE)

ps_low_prev_u30 <- merge_phyloseq(ps_rel, asv_low_prev_u30)

# NMDS
bc_dist_low_prev_u30 <- phyloseq::distance(ps_low_prev_u30, method = "bray")
nmds_low_prev_u30 <- ordinate(
  ps_low_prev_u30,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

# visualization
nmds_plot_low_prev_u30 <- plot_ordination(
  ps_low_prev_u30,
  nmds_low_prev_u30,
  color = "Plant",
  shape = "Location"
) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~Date, ncol = 2) +
  labs(
    title = "NMDS of Low Prevalence OTUs (Present in <30% of Samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_low_prev_u30)


### Core Microbiome
## Corn, plot by plot

ps_corn <- subset_samples(ps_rel, Plant == "Corn")
ps_corn
sample_data(ps_corn)
otu_table(ps_corn)
tax_table(ps_corn)

# Phyloseq objects for each Corn plot
ps_corn_P01 <- subset_samples(ps_corn, Plot == "P01")
ps_corn_P08 <- subset_samples(ps_corn, Plot == "P08")
ps_corn_P11 <- subset_samples(ps_corn, Plot == "P11")
ps_corn_P14 <- subset_samples(ps_corn, Plot == "P14")

# Corn-Plot_Location-Interrow
ps_corn_P01_I <- subset_samples(ps_corn_P01, Location == "Interrow")
ps_corn_P08_I <- subset_samples(ps_corn_P08, Location == "Interrow")
ps_corn_P11_I <- subset_samples(ps_corn_P11, Location == "Interrow")
ps_corn_P14_I <- subset_samples(ps_corn_P14, Location == "Interrow")

# Corn-Plot_Location-Row
ps_corn_P01_R <- subset_samples(ps_corn_P01, Location == "Row")
ps_corn_P08_R <- subset_samples(ps_corn_P08, Location == "Row")
ps_corn_P11_R <- subset_samples(ps_corn_P11, Location == "Row")
ps_corn_P14_R <- subset_samples(ps_corn_P14, Location == "Row")

# 생성된 객체 확인
ps_corn_P01_I
ps_corn_P08_I
ps_corn_P11_I
ps_corn_P14_I
ps_corn_P01_R
ps_corn_P08_R
ps_corn_P11_R
ps_corn_P14_R

## P01_I
# ASV table
asv_table_data_P01_I <- as.data.frame(otu_table(ps_corn_P01_I))
dim(asv_table_data_P01_I)

# num of samples for P01_I
total_samples_P01_I <- nrow(asv_table_data_P01_I)

# counting number of samples containing each OTU
asv_sample_counts_P01_I <- colSums(asv_table_data_P01_I > 0)

# ASVs found in >80% samples
threshold_80_P01_I <- total_samples_P01_I * 0.8
core_asvs_P01_I <- names(asv_sample_counts_P01_I[
  asv_sample_counts_P01_I >= threshold_80_P01_I
])

cat("ASVs found in >80% P01_I samples:", length(core_asvs_P01_I), "\n")
head(asv_table_data_P01_I[, core_asvs_P01_I])

## P01_I, P08_I
ps_corn_P01_P08_I <- merge_phyloseq(ps_corn_P01_I, ps_corn_P08_I)

asv_table_data_P01_P08_I <- as.data.frame(otu_table(ps_corn_P01_P08_I))

total_samples_P01_P08_I <- nrow(asv_table_data_P01_P08_I)

asv_sample_counts_P01_P08_I <- colSums(asv_table_data_P01_P08_I > 0)

threshold_80_P01_P08_I <- total_samples_P01_P08_I * 0.8
core_asvs_P01_P08_I <- names(asv_sample_counts_P01_P08_I[
  asv_sample_counts_P01_P08_I >= threshold_80_P01_P08_I
])

cat(
  "ASVs found in >80% P01_I + P08_I samples:",
  length(core_asvs_P01_P08_I),
  "\n"
)

# NMDS for P01-I and P08-I
nmds_P01_P08_I <- ordinate(
  ps_corn_P01_P08_I,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

nmds_plot_P01_P08_I <- plot_ordination(
  ps_corn_P01_P08_I,
  nmds_P01_P08_I,
  color = "Plot",
  shape = "Date"
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "NMDS (P01-IR + P08-IR, 13 samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_P01_P08_I)


## P01_I, P08_I, P11_I
ps_corn_P01_P08_P11_I <- merge_phyloseq(
  ps_corn_P01_I,
  ps_corn_P08_I,
  ps_corn_P11_I
)

asv_table_data_P01_P08_P11_I <- as.data.frame(otu_table(ps_corn_P01_P08_P11_I))

total_samples_P01_P08_P11_I <- nrow(asv_table_data_P01_P08_P11_I)

asv_sample_counts_P01_P08_P11_I <- colSums(asv_table_data_P01_P08_P11_I > 0)

threshold_80_P01_P08_P11_I <- total_samples_P01_P08_P11_I * 0.8
core_asvs_P01_P08_P11_I <- names(asv_sample_counts_P01_P08_P11_I[
  asv_sample_counts_P01_P08_P11_I >= threshold_80_P01_P08_P11_I
])

cat(
  "ASVs found in >80% P01_I + P08_I + P11_I samples:",
  length(core_asvs_P01_P08_P11_I),
  "\n"
)

# NMDS for P01_I, P08_I, P11_I
nmds_P01_P08_P11_I <- ordinate(
  ps_corn_P01_P08_P11_I,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

nmds_plot_P01_P08_P11_I <- plot_ordination(
  ps_corn_P01_P08_P11_I,
  nmds_P01_P08_P11_I,
  color = "Plot",
  shape = "Date"
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "NMDS (P01-IR + P08-IR + P11-IR, 20 samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_P01_P08_P11_I)


## P01_I, P08_I, P11_I, P14_I
ps_corn_P01_P08_P11_P14_I <- merge_phyloseq(
  ps_corn_P01_I,
  ps_corn_P08_I,
  ps_corn_P11_I,
  ps_corn_P14_I
)

asv_table_data_P01_P08_P11_P14_I <- as.data.frame(otu_table(
  ps_corn_P01_P08_P11_P14_I
))
total_samples_P01_P08_P11_P14_I <- nrow(asv_table_data_P01_P08_P11_P14_I)

asv_sample_counts_P01_P08_P11_P14_I <- colSums(
  asv_table_data_P01_P08_P11_P14_I > 0
)

threshold_80_P01_P08_P11_P14_I <- total_samples_P01_P08_P11_P14_I * 0.8
core_asvs_P01_P08_P11_P14_I <- names(asv_sample_counts_P01_P08_P11_P14_I[
  asv_sample_counts_P01_P08_P11_P14_I >= threshold_80_P01_P08_P11_P14_I
])

cat(
  "ASVs found in >80% P01_I + P08_I + P11_I + P14_I samples:",
  length(core_asvs_P01_P08_P11_P14_I),
  "\n"
)

# NMDS
nmds_P01_P08_P11_P14_I <- ordinate(
  ps_corn_P01_P08_P11_P14_I,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

nmds_plot_P01_P08_P11_P14_I <- plot_ordination(
  ps_corn_P01_P08_P11_P14_I,
  nmds_P01_P08_P11_P14_I,
  color = "Plot",
  shape = "Date"
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "NMDS (P01-IR + P08-IR + P11-IR + P14-IR, 27 samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_P01_P08_P11_P14_I)


### Corn-Row 추가
ps_corn_P01_R
ps_corn_P08_R
ps_corn_P11_R
ps_corn_P14_R

## Corn_Interrow + P01_R
ps_corn_I_P01_R <- merge_phyloseq(
  ps_corn_P01_I,
  ps_corn_P08_I,
  ps_corn_P11_I,
  ps_corn_P14_I,
  ps_corn_P01_R
)

asv_table_data_corn_I_P01_R <- as.data.frame(otu_table(ps_corn_I_P01_R))

total_samples_corn_I_P01_R <- nrow(asv_table_data_corn_I_P01_R)

asv_sample_counts_corn_I_P01_R <- colSums(asv_table_data_corn_I_P01_R > 0)

threshold_80_corn_I_P01_R <- total_samples_corn_I_P01_R * 0.8
core_asvs_corn_I_P01_R <- names(asv_sample_counts_corn_I_P01_R[
  asv_sample_counts_corn_I_P01_R >= threshold_80_corn_I_P01_R
])

cat("ASVs found in >80% corn_I+P01_RI", length(core_asvs_corn_I_P01_R), "\n")

# NMDS
nmds_corn_I_P01_R <- ordinate(
  ps_corn_I_P01_R,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

nmds_plot_corn_I_P01_R <- plot_ordination(
  ps_corn_I_P01_R,
  nmds_corn_I_P01_R,
  color = "Plot",
  shape = "Date"
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "NMDS (Corn Interrow + P01-R, 33 samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_corn_I_P01_R)


## Corn_Interrow + P01_R + P08_R
ps_corn_I_P01_P08_R <- merge_phyloseq(
  ps_corn_P01_I,
  ps_corn_P08_I,
  ps_corn_P11_I,
  ps_corn_P14_I,
  ps_corn_P01_R,
  ps_corn_P08_R
)

asv_table_data_corn_I_P01_P08_R <- as.data.frame(otu_table(ps_corn_I_P01_P08_R))

total_samples_corn_I_P01_P08_R <- nrow(asv_table_data_corn_I_P01_P08_R)

asv_sample_counts_corn_I_P01_P08_R <- colSums(
  asv_table_data_corn_I_P01_P08_R > 0
)

threshold_80_corn_I_P01_P08_R <- total_samples_corn_I_P01_P08_R * 0.8
core_asvs_corn_I_P01_P08_R <- names(asv_sample_counts_corn_I_P01_P08_R[
  asv_sample_counts_corn_I_P01_P08_R >= threshold_80_corn_I_P01_P08_R
])

cat(
  "ASVs found in >80% corn_Interrow + P01_R + P08_R samples:",
  length(core_asvs_corn_I_P01_P08_R),
  "\n"
)

# NMDS
nmds_corn_I_P01_R_P08_R <- ordinate(
  ps_corn_I_P01_P08_R,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

nmds_plot_corn_I_P01_R_P08_R <- plot_ordination(
  ps_corn_I_P01_P08_R,
  nmds_corn_I_P01_R_P08_R,
  color = "Plot",
  shape = "Date"
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "NMDS (Corn Interrow + P01-R + P08-R, 40 samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_corn_I_P01_R_P08_R)

## Corn_Interrow + P01_R + P08_R + P11_R
ps_corn_I_P01_P08_P11_R <- merge_phyloseq(
  ps_corn_P01_I,
  ps_corn_P08_I,
  ps_corn_P11_I,
  ps_corn_P14_I,
  ps_corn_P01_R,
  ps_corn_P08_R,
  ps_corn_P11_R
)

asv_table_data_corn_I_P01_P08_P11_R <- as.data.frame(otu_table(
  ps_corn_I_P01_P08_P11_R
))

total_samples_corn_I_P01_P08_P11_R <- nrow(asv_table_data_corn_I_P01_P08_P11_R)

asv_sample_counts_corn_I_P01_P08_P11_R <- colSums(
  asv_table_data_corn_I_P01_P08_P11_R > 0
)

threshold_80_corn_I_P01_P08_P11_R <- total_samples_corn_I_P01_P08_P11_R * 0.8
core_asvs_corn_I_P01_P08_P11_R <- names(asv_sample_counts_corn_I_P01_P08_P11_R[
  asv_sample_counts_corn_I_P01_P08_P11_R >= threshold_80_corn_I_P01_P08_P11_R
])

cat(
  "ASVs found in >80% corn_I+P01_P08_P11_R samples:",
  length(core_asvs_corn_I_P01_P08_P11_R),
  "\n"
)

# NMDS
nmds_corn_I_P01_R_P08_R_P11_R <- ordinate(
  ps_corn_I_P01_P08_P11_R,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

nmds_plot_corn_I_P01_R_P08_R_P11_R <- plot_ordination(
  ps_corn_I_P01_P08_P11_R,
  nmds_corn_I_P01_R_P08_R_P11_R,
  color = "Plot",
  shape = "Date"
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "NMDS (Corn Interrow + P01-R + P08-R + P11-R, 47 samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_corn_I_P01_R_P08_R_P11_R)


## Corn_Interrow + P01_R + P08_R + P11_R + P14_R
ps_corn_I_P01_P08_P11_P14_R <- merge_phyloseq(
  ps_corn_P01_I,
  ps_corn_P08_I,
  ps_corn_P11_I,
  ps_corn_P14_I,
  ps_corn_P01_R,
  ps_corn_P08_R,
  ps_corn_P11_R,
  ps_corn_P14_R
)

asv_table_data_corn_I_P01_P08_P11_P14_R <- as.data.frame(otu_table(
  ps_corn_I_P01_P08_P11_P14_R
))

total_samples_corn_I_P01_P08_P11_P14_R <- nrow(
  asv_table_data_corn_I_P01_P08_P11_P14_R
)

asv_sample_counts_corn_I_P01_P08_P11_P14_R <- colSums(
  asv_table_data_corn_I_P01_P08_P11_P14_R > 0
)

threshold_80_corn_I_P01_P08_P11_P14_R <- total_samples_corn_I_P01_P08_P11_P14_R *
  0.8
core_asvs_corn_I_P01_P08_P11_P14_R <- names(asv_sample_counts_corn_I_P01_P08_P11_P14_R[
  asv_sample_counts_corn_I_P01_P08_P11_P14_R >=
    threshold_80_corn_I_P01_P08_P11_P14_R
])

cat(
  "ASVs found in >80% corn_I+R:",
  length(core_asvs_corn_I_P01_P08_P11_P14_R),
  "\n"
)

# NMDS
nmds_corn_I_R <- ordinate(
  ps_corn_I_P01_P08_P11_P14_R,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

nmds_plot_corn_I_R <- plot_ordination(
  ps_corn_I_P01_P08_P11_P14_R,
  nmds_corn_I_R,
  color = "Plot",
  shape = "Date"
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "NMDS (Corn Interrow + Row, 53 samples)",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(nmds_plot_corn_I_R)

#### number count - diffrernt thesholds
threshold_80_P01_I <- total_samples_P01_I * 0.7
core_asvs_P01_I <- names(asv_sample_counts_P01_I[
  asv_sample_counts_P01_I >= threshold_80_P01_I
])
cat("ASVs found in >80% P01_I samples:", length(core_asvs_P01_I), "\n")

threshold_80_P01_P08_I <- total_samples_P01_P08_I * 0.7
core_asvs_P01_P08_I <- names(asv_sample_counts_P01_P08_I[
  asv_sample_counts_P01_P08_I >= threshold_80_P01_P08_I
])
cat(
  "ASVs found in >80% P01_I + P08_I samples:",
  length(core_asvs_P01_P08_I),
  "\n"
)

threshold_80_P01_P08_P11_I <- total_samples_P01_P08_P11_I * 0.7
core_asvs_P01_P08_P11_I <- names(asv_sample_counts_P01_P08_P11_I[
  asv_sample_counts_P01_P08_P11_I >= threshold_80_P01_P08_P11_I
])
cat(
  "ASVs found in >80% P01_I + P08_I + P11_I samples:",
  length(core_asvs_P01_P08_P11_I),
  "\n"
)

threshold_80_P01_P08_P11_P14_I <- total_samples_P01_P08_P11_P14_I * 0.7
core_asvs_P01_P08_P11_P14_I <- names(asv_sample_counts_P01_P08_P11_P14_I[
  asv_sample_counts_P01_P08_P11_P14_I >= threshold_80_P01_P08_P11_P14_I
])
cat(
  "ASVs found in >80% P01_I + P08_I + P11_I + P14_I samples:",
  length(core_asvs_P01_P08_P11_P14_I),
  "\n"
)

threshold_80_corn_I_P01_R <- total_samples_corn_I_P01_R * 0.7
core_asvs_corn_I_P01_R <- names(asv_sample_counts_corn_I_P01_R[
  asv_sample_counts_corn_I_P01_R >= threshold_80_corn_I_P01_R
])
cat("ASVs found in >80% corn_I+P01_RI", length(core_asvs_corn_I_P01_R), "\n")

threshold_80_corn_I_P01_P08_R <- total_samples_corn_I_P01_P08_R * 0.7
core_asvs_corn_I_P01_P08_R <- names(asv_sample_counts_corn_I_P01_P08_R[
  asv_sample_counts_corn_I_P01_P08_R >= threshold_80_corn_I_P01_P08_R
])
cat(
  "ASVs found in >80% corn_Interrow + P01_R + P08_R samples:",
  length(core_asvs_corn_I_P01_P08_R),
  "\n"
)

threshold_80_corn_I_P01_P08_P11_R <- total_samples_corn_I_P01_P08_P11_R * 0.7
core_asvs_corn_I_P01_P08_P11_R <- names(asv_sample_counts_corn_I_P01_P08_P11_R[
  asv_sample_counts_corn_I_P01_P08_P11_R >= threshold_80_corn_I_P01_P08_P11_R
])
cat(
  "ASVs found in >80% corn_I+P01_P08_P11_R samples:",
  length(core_asvs_corn_I_P01_P08_P11_R),
  "\n"
)

threshold_80_corn_I_P01_P08_P11_P14_R <- total_samples_corn_I_P01_P08_P11_P14_R *
  0.7
core_asvs_corn_I_P01_P08_P11_P14_R <- names(asv_sample_counts_corn_I_P01_P08_P11_P14_R[
  asv_sample_counts_corn_I_P01_P08_P11_P14_R >=
    threshold_80_corn_I_P01_P08_P11_P14_R
])
cat(
  "ASVs found in >80% corn_I+R:",
  length(core_asvs_corn_I_P01_P08_P11_P14_R),
  "\n"
)


####

# prevalance vs abundance

# all ASVs
asv_table_data <- as.data.frame(otu_table(ps_rel))
asv_table_data

write.csv(
  asv_table_data,
  file = "/Users/jaejinlee/Files/Data/2023SABR_amplicon/analysis/asv_table_ps_rel_all.csv",
  row.names = TRUE
)

total_samples <- nrow(asv_table_data)

# calculating prevalence
asv_prevalence <- colSums(asv_table_data > 0) / total_samples * 100

# mean adundance for each ASV
asv_abundance <- colMeans(asv_table_data)

# creating a data frame
prevalence_abundance_df <- data.frame(
  OTU = names(asv_prevalence),
  Prevalence = asv_prevalence,
  Abundance = asv_abundance
)

library(scales)
# ggplot(adding 90% line)
ggplot(prevalence_abundance_df, aes(x = Prevalence, y = Abundance)) +
  geom_point(alpha = 0.6, color = "gray60") + # 전체 OTU (회색)
  geom_point(
    data = prevalence_abundance_df[prevalence_abundance_df$Prevalence >= 90, ],
    aes(x = Prevalence, y = Abundance),
    color = "blue",
    alpha = 0.8
  ) + # dots in blue
  scale_y_log10(labels = function(y) round(log10(y), 0)) + # log value
  geom_vline(xintercept = 90, linetype = "dashed", color = "red", size = 1) + # adding 90% prevalence line
  labs(
    title = "Prevalence vs. Abundance (All ASVs, 90% Core Highlighted)",
    x = "Prevalence (%)",
    y = "Mean Relative Abundance (Log scale)"
  ) +
  theme_minimal()


# list of blue dots (>90% prevalence)
blue_asvs <- prevalence_abundance_df$OTU[
  prevalence_abundance_df$Prevalence >= 90
]

cat("ASVs >90% corn samples:", length(blue_asvs), "\n")

# taxnoomy table
tax_table_data <- as.data.frame(tax_table(ps_rel))

# blue_asv Taxonomy
blue_asvs_taxa <- tax_table_data[blue_asvs, ]
print(blue_asvs_taxa)

# converting to a data framce
tax_table_data <- as.data.frame(tax_table(ps_rel))

blue_asvs_taxa <- tax_table_data[blue_asvs, , drop = FALSE]
asv_table_data <- as.data.frame(otu_table(ps_rel))

total_samples <- nrow(asv_table_data)

asv_prevalence <- colSums(asv_table_data > 0) / total_samples * 100 # 백분율(%)
asv_abundance <- colMeans(asv_table_data)

blue_asvs_df <- blue_asvs_taxa %>%
  mutate(
    OTU_Sequence = rownames(blue_asvs_taxa), # OTU 시퀀스
    Prevalence = asv_prevalence[blue_asvs], # Prevalence (%)
    Mean_Relative_Abundance = asv_abundance[blue_asvs] # Mean Relative Abundance
  ) %>%
  select(OTU_Sequence, everything()) # OTU 시퀀스를 첫 컬럼으로 정렬

print(blue_asvs_df)

# save as a CSV file
write.csv(
  blue_asvs_df,
  file = "/Users/jaejinlee/Files/Data/2023SABR_amplicon/analysis/blue_asvs_table.csv",
  row.names = FALSE
)


### Considering both prevalance (>80% samples) and abundance (top 1%)
# Prevalence calculation
asv_prevalence_1 <- colSums(asv_table_data > 0) / total_samples * 100

# Abundance calculation
asv_abundance_1 <- colMeans(asv_table_data)

# Calculate the threshold for top 20% abundance
abundance_threshold_1 <- quantile(asv_abundance_1, 0.99)

# Create a data frame with Prevalence and Abundance
prevalence_abundance_df_1 <- data.frame(
  OTU = names(asv_prevalence_1),
  Prevalence = asv_prevalence_1,
  Abundance = asv_abundance_1
)

# Plot Prevalence vs Abundance with both criteria highlighted
ggplot(prevalence_abundance_df_1, aes(x = Prevalence, y = Abundance)) +
  geom_point(alpha = 0.6, color = "gray60") +
  geom_point(
    data = prevalence_abundance_df_1[
      prevalence_abundance_df_1$Prevalence >= 80 &
        prevalence_abundance_df_1$Abundance >= abundance_threshold_1,
    ],
    aes(x = Prevalence, y = Abundance),
    color = "blue",
    alpha = 0.8
  ) +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(.x))
  ) + # log10 value
  geom_vline(xintercept = 80, linetype = "dashed", color = "red", size = 1) +
  geom_hline(
    yintercept = abundance_threshold_1,
    linetype = "dashed",
    color = "green",
    size = 1
  ) +
  labs(
    title = "Prevalence vs. Abundance (80% Prevalence & Top 1% Abundance Highlighted)",
    x = "Prevalence (%)",
    y = "Mean Relative Abundance (Log scale)"
  ) +
  theme_minimal()


# Extract OTUs that meet both criteria (Prevalence ≥ 80% and Top 20% Abundance)
selected_otus_1 <- prevalence_abundance_df_1$OTU[
  prevalence_abundance_df_1$Prevalence >= 80 &
    prevalence_abundance_df_1$Abundance >= abundance_threshold_1
]

# Count the number of selected OTUs
cat(
  "OTUs meeting both criteria (Prevalence ≥ 80% & Top 20% Abundance):",
  length(selected_otus_1),
  "\n"
)

# Retrieve taxonomy information for selected OTUs
selected_otus_taxa_1 <- tax_table_data[selected_otus_1, , drop = FALSE]

# Create a final data frame with taxonomy, prevalence, and abundance
selected_otus_df_1 <- selected_otus_taxa_1 %>%
  mutate(
    OTU_Sequence = rownames(selected_otus_taxa_1), # OTU sequence
    Prevalence = asv_prevalence_1[selected_otus_1], # Prevalence (%)
    Mean_Relative_Abundance = asv_abundance_1[selected_otus_1] # Mean relative abundance
  ) %>%
  select(OTU_Sequence, everything()) # Move OTU sequence to the first column

# Print the final data frame
print(selected_otus_df_1)

# Save the final table as a CSV file
write.csv(
  selected_otus_df_1,
  file = "/Users/jaejinlee/Files/Data/2023SABR_amplicon/selected_otus_table_1.csv",
  row.names = FALSE
)


#### heatmap type analysis
### compareing asv lists
#  각 단계의 ASV 리스트 추출
asv_P01_I <- core_asvs_P01_I
asv_P01_P08_I <- core_asvs_P01_P08_I
asv_P01_P08_P11_I <- core_asvs_P01_P08_P11_I
asv_P01_P08_P11_P14_I <- core_asvs_P01_P08_P11_P14_I

# 각 단계에서 제거된 ASV 찾기
removed_P08 <- setdiff(asv_P01_I, asv_P01_P08_I)
removed_P11 <- setdiff(asv_P01_P08_I, asv_P01_P08_P11_I)
removed_P14 <- setdiff(asv_P01_P08_P11_I, asv_P01_P08_P11_P14_I)

# 혹시 새로 추가된 ASV 있는지 확인
added_P08 <- setdiff(asv_P01_P08_I, asv_P01_I)
added_P11 <- setdiff(asv_P01_P08_P11_I, asv_P01_P08_I)
added_P14 <- setdiff(asv_P01_P08_P11_P14_I, asv_P01_P08_P11_I)

# 결과 출력
cat("P01 → P01+P08 단계에서 제거된 ASV 개수:", length(removed_P08), "\n")
cat(
  "P01+P08 → P01+P08+P11 단계에서 제거된 ASV 개수:",
  length(removed_P11),
  "\n"
)
cat(
  "P01+P08+P11 → P01+P08+P11+P14 단계에서 제거된 ASV 개수:",
  length(removed_P14),
  "\n"
)

cat("\n각 단계에서 제거된 ASV 목록:\n")
cat("Removed from P08:", removed_P08, "\n\n")
cat("Removed from P11:", removed_P11, "\n\n")
cat("Removed from P14:", removed_P14, "\n\n")

cat("\n혹시 새롭게 추가된 ASV가 있는지 확인:\n")
cat("Added in P08:", added_P08, "\n\n")
cat("Added in P11:", added_P11, "\n\n")
cat("Added in P14:", added_P14, "\n\n")

cat("P01 → P01+P08 단계에서 추가된 ASV 개수:", length(added_P08), "\n")
cat("P01+P08 → P01+P08+P11 단계에서 추가된 ASV 개수:", length(added_P11), "\n")
cat(
  "P01+P08+P11 → P01+P08+P11+P14 단계에서 추가된 ASV 개수:",
  length(added_P14),
  "\n"
)


## plot별 core asv 비교
library(UpSetR)
library(reshape2)

core_asvs_list <- list(
  P01 = core_asvs_P01_I,
  P08 = core_asvs_P01_P08_I,
  P11 = core_asvs_P01_P08_P11_I,
  P14 = core_asvs_P01_P08_P11_P14_I
)
core_asvs_list
# 각 ASV가 어떤 Plot에서 존재하는지 0/1 매트릭스로 변환
asv_matrix <- data.frame(
  ASV = unique(unlist(core_asvs_list)), # 모든 unique ASV 리스트
  P01 = as.integer(unique(unlist(core_asvs_list)) %in% core_asvs_P01_I),
  P08 = as.integer(unique(unlist(core_asvs_list)) %in% core_asvs_P01_P08_I),
  P11 = as.integer(unique(unlist(core_asvs_list)) %in% core_asvs_P01_P08_P11_I),
  P14 = as.integer(
    unique(unlist(core_asvs_list)) %in% core_asvs_P01_P08_P11_P14_I
  )
)
asv_matrix$P14


# ASV를 rowname으로 설정
rownames(asv_matrix) <- asv_matrix$ASV
asv_matrix <- asv_matrix[, -1] # ASV 컬럼 제거 (UpSetR에 필요 없음)

# 필요한 패키지 로드
library(tidyverse)

# 데이터 프레임을 long format으로 변환하면서 ASV 이름 변경
asv_matrix_long <- asv_matrix %>%
  rownames_to_column(var = "ASV") %>% # 기존 ASV 이름 유지
  mutate(ASV_new = sprintf("C_ASV%03d", row_number())) %>% # 새로운 ASV 이름 생성
  pivot_longer(
    cols = -c(ASV, ASV_new),
    names_to = "Sample",
    values_to = "Presence"
  ) %>%
  mutate(ASV_new = factor(ASV_new, levels = unique(.$ASV_new))) # 원래 순서 유지

# 히트맵 생성
ggplot(asv_matrix_long, aes(x = Sample, y = ASV_new, fill = factor(Presence))) +
  geom_tile(color = "gray70", size = 0.3) + # ASV별 구분선 추가
  scale_fill_manual(values = c("0" = "white", "1" = "skyblue")) +
  scale_x_discrete(
    labels = c(
      "P01" = "P01 (7)",
      "P08" = "P01+P08 (13)",
      "P11" = "P01+P08+P11 (19)",
      "P14" = "P01+P08+P11+P14 (26)"
    )
  ) + # X축 샘플 이름 변경
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), # 전체 제목 크기 조정 (20)
    axis.text.y = element_text(
      size = 8,
      hjust = 1,
      vjust = 0.5,
      margin = margin(r = -20)
    ), # Y축 레이블과 그래프 사이 간격 줄이기
    axis.text.x = element_text(size = 14, angle = 0, vjust = 0.5, hjust = 0.5), # X축 샘플 이름 크기 조정
    axis.title.x = element_text(size = 16, face = "bold"), # X축 제목 크기 조정
    axis.title.y = element_text(size = 16, face = "bold"), # Y축 제목 크기 조정
    axis.ticks.y = element_blank(), # Y축 눈금 제거
    panel.grid = element_blank(), # 배경 격자 제거
    legend.position = "none" # 범례 제거
  ) +
  labs(
    title = "Core ASV (existing >80% samples)",
    x = "Corn-Interrow Plots",
    y = "ASV"
  )
