#' Create FASTA format from ASV sequences in a DADA2 sequence table
#'
#' This function converts the ASV sequences from a DADA2 sequence table
#' into a FASTA format vector with headers in the form ">ASV_1", ">ASV_2", etc.
#'
#' @param seqtab.nochim A DADA2 sequence table (the output from removeBimeraDenovo)
#' @return A character vector representing FASTA format with alternating headers and sequences
#' @examples
#' # Assuming you have a seqtab.nochim object from DADA2:
#' fasta_output <- seqtab2fasta(seqtab.nochim)
#' # Write to a file:
#' # writeLines(fasta_output, "asv_sequences.fasta")

seqtab2fasta <- function(seqtab.nochim) {
  # Extract ASV sequences from column names
  asv_seqs <- colnames(seqtab.nochim)

  # Create ASV headers
  asv_headers <- vector(dim(seqtab.nochim)[2], mode = "character")
  for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep = "_")
  }

  # Interleave headers and sequences to create FASTA format
  asv_fasta <- c(rbind(asv_headers, asv_seqs))

  return(list(
    asv_fasta = asv_fasta,
    asv_headers = asv_headers
  ))
}
