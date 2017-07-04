
## write a fasta file from a sequence of files
write_fasta <- function(sequences, names, filepath) {

  maxwidth <- 20000
  if (any(nchar(sequences) >= maxwidth)) {
    warning("sequence will be written out in multiple lines. Max seq. length = ", max(nchar(sequences)),
            " > ", maxwidth)
  }
  
  seq <- sequences
  names(seq) <- names
  obj <- Biostrings::DNAStringSet(seq)

  Biostrings::writeXStringSet(obj, filepath = filepath, width = maxwidth)
}


##' Convert data.table containing one sequence per row into a DNAStringSet
##' @param dt data.table containing columns id_col and seq_col
##' @param id_col column name for the identifier (has to be unique)
dt2DNAStringSet <- function(dt, id_col = "ID", seq_col = "seq") {
  set <- dt[, seq_col, with = F][[1]] %>% as("DNAStringSet")
  rownames <- dt[, id_col, with = F][[1]]
  stopifnot(!any(duplicated(rownames)))
  names(set) <- rownames
  attr(set, "id_col") <- id_col
  attr(set, "seq_col") <- seq_col
  return(set)
}


##' Import narrowPeak .bed files
##'
##' From: https://charlesjb.github.io/How_to_import_narrowPeak/
##' 
##' @param file bed file path
##' @return GRanges
import_narrowPeak <- function(file) {
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                            qValue = "numeric", peak = "integer")
  rtracklayer::import(file, format = "BED", extraCols = extraCols_narrowPeak)
}

##' Import narrowPeak .bed files
##'
##' From: https://charlesjb.github.io/How_to_import_narrowPeak/
##' 
##' @param file bed file path
##' @return GRanges
import_broadPeak <- function(file) {
  extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric",
                           qValue = "numeric")
  rtracklayer::import(file, format = "BED", extraCols = extraCols_broadPeak)
}
