#####################################################################
# NMDS Analysis of Microbial Communities
#
# This script performs Non-Metric Multidimensional Scaling (NMDS)
# analysis on the phyloseq object to visualize the microbial community
# structure across different samples.
#
# Author: Bolívar Aponte Rolón
# Date: 2025-05-05
#####################################################################

# Load required libraries
library(phyloseq)
library(ggplot2)

# Ensure the phyloseq object (ps) is loaded before running this script

#--------------------------------------------------------
# Calculate Bray-Curtis distance and perform NMDS analysis
#--------------------------------------------------------

# Calculate Bray-Curtis distance matrix
bc_dist <- phyloseq::distance(ps, method = "bray")

# Perform NMDS analysis (k=2 reduces to 2 dimensions)
# trymax=100 ensures convergence by trying up to 100 different random starts
nmds <- ordinate(ps, method = "NMDS", distance = "bray", trymax = 100)

#--------------------------------------------------------
# Prepare data for visualization
#--------------------------------------------------------

# Convert "Date" to a factor for proper panel separation in visualization
sample_data(ps)$Date <- as.factor(sample_data(ps)$Date)
colnames(sample_data(ps)) # Display column names for reference

#--------------------------------------------------------
# Generate and display NMDS plot
#--------------------------------------------------------

# Create NMDS plot with samples colored by Plant, shaped by Location, and faceted by Date
nmds_plot <- plot_ordination(ps, nmds, color = "Plant", shape = "Location") +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~Date, ncol = 2) + # Separate panels by Date
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

# Display the plot
print(nmds_plot)

# Optionally, save the plot
# ggsave("outputs/nmds_plot.pdf", nmds_plot, width = 10, height = 8)
