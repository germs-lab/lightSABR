# Function to clean sample names in rownames
clean_rownames <- function(df) {
  # Get current rownames
  old_rownames <- rownames(df)

  # Create new cleaned rownames
  new_rownames <- gsub("^(.*?)_S\\d+.*$", "\\1", old_rownames)

  # Assign new rownames
  rownames(df) <- new_rownames

  return(df)
}
