##' Generate counts for all the k-mers in the sequence
##' @param seq character vector of letters A, C, G, T
##' @param k Number of bases used to generate 
get_kmers_X <- function(seq, k) {
  library(Biostrings)
  dseq <- DNAStringSet(seq)
  X <- oligonucleotideFrequency(dseq, k)
  return(X)
}

##' str_pad
pad_start <- function(seq, max_n = max(nchar(seq))) {
  stopifnot(max_n >= max(nchar(seq)))
  seq <- stringr::str_pad(seq, max_n, side = "left", pad = "N")
  ## seq <- sapply(seq, function(word) paste0(word, NN, collapse = ''))
  
  return(seq)
}
