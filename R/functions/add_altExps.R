#' Add prevalent and rare taxa subsets as alternative experiments to a TreeSummarizedExperiment
#'
#' Iterates over specified prevalence thresholds and adds both prevalent and rare taxa
#' subsets as alternative experiments to the provided TreeSummarizedExperiment object.
#'
#' @param ps A TreeSummarizedExperiment (from mia) object.
#' @param prev_name Optional string prefix for prevalent altExp names. If NULL, uses "prevalent".
#' @param rare_name Optional string prefix for rare altExp names. If NULL, uses "rare".
#' @param thresholds Numeric vector of prevalence thresholds (as percentages, e.g. c(90, 80, 70, 60)).
#' @param rank Taxonomic rank for prevalence calculation (default "phylum").
#' @param rank_name Optional string to append to altExp names, e.g. "Phylum" (default NULL).
#' @param detection Detection threshold for mia functions (default 1/100).
#' @param assay.type Assay type for mia functions (default "counts").
#'
#' @return TreeSummarizedExperiment with new altExps added for each threshold.
#' @examples
#' ps <- add_prevalent_rare_altExps(ps)
add_prevalent_rare_altExps <- function(
  ps,
  prev_name = NULL,
  rare_name = NULL,
  thresholds = c(90, 80, 70, 60),
  rank = "phylum",
  rank_name = NULL,
  detection = 1 / 100,
  assay.type = "counts"
) {
  for (thr in thresholds) {
    # Build altExp names
    prev_suffix <- if (!is.null(rank_name)) {
      paste0("_", rank_name, "_", thr)
    } else {
      paste0("_", thr)
    }
    prev_exp_name <- if (!is.null(prev_name)) {
      paste0(prev_name, prev_suffix)
    } else {
      paste0("prevalent", prev_suffix)
    }
    rare_exp_name <- if (!is.null(rare_name)) {
      paste0(rare_name, prev_suffix)
    } else {
      paste0("rare", prev_suffix)
    }

    # Add prevalent taxa subset
    altExp(ps, prev_exp_name) <- subsetByPrevalent(
      ps,
      rank = rank,
      assay.type = assay.type,
      detection = detection,
      prevalence = thr / 100
    )

    # Add rare taxa subset (inverse)
    altExp(ps, rare_exp_name) <- subsetByRare(
      ps,
      rank = rank,
      assay.type = assay.type,
      detection = detection,
      prevalence = thr / 100
    )
  }
  return(ps)
}
