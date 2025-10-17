# Safely convert a phyloseq object to long relative-abundance form
relab_prevalence_long <- function(ps, prevalence_label) {
  if (is.null(ps)) {
    return(NULL)
  }
  # Transform to relative abundance
  ps_rel <- transform_sample_counts(ps, function(x) {
    if (sum(x) == 0) x else x / sum(x)
  })

  df <- psmelt(ps_rel) %>%
    mutate(
      prevalence_level = prevalence_label
    ) %>%
    # Keep only needed columns; rename for clarity
    rename(
      sample_id = Sample,
      asv = OTU,
      rel_abundance = Abundance
    )
}
