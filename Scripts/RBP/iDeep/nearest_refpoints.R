#'---
#' Title: Generate positions.csv
#' Author: Ziga Avsec
#'---
## Input:
## - fasta with ranges embeded in the names
## - version=hg19
##
## Output:
## - csv file with the same number of rows as the bed file

version <- "hg19"
source("../../../.Rprofile", chdir=TRUE)
source("functions.R")


gr <- get_human_anno(version=version, use_mito=TRUE)
gr_positions <- get_position_landmarks(gr)

generate_position_csv <- function(in_fa, out_csv) {
  gr_peak <- get_peak_centers(in_fa)

  dt_dist <- sapply(gr_positions, function(gr_anno) {
    ids <- GenomicRanges::nearest(gr_peak, gr_anno)
    pos_nearest <- start(gr_anno)[ids]
    rel_distance <- ifelse(strand(gr_peak) == "+", 1, -1) * (pos_nearest - start(gr_peak))
    return(rel_distance)  
  }, simplify=FALSE) %>% as.data.table

  dt_dist[, class := gr_peak$class]
  write_csv(dt_dist, out_csv)  
}

generate_position_csv(in_fa = snakemake@input[["fa"]], out_csv = snakemake@output[["csv"]])

## pblapply(dir_list, generate_position_csv)

## --------------------------------------------
