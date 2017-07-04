#'---
#' Title: script
#'---
#!/usr/bin/env R

## Input:
## - bed
## - version=hg19
##
## Output:
## - csv file with the same number of rows as the bed file

version <- "hg19"
source("Scripts/RBP/iDeep/functions.R")

gr <- get_human_anno(version=version, use_mito=TRUE)

gr_positions <- get_position_landmarks(gr)


## list all the directories
dir_list <- list.dirs("~/github/iDeep/datasets/clip/") %>% grep("sample_0$", ., value=TRUE)

generate_position_csv <- function(dir) {
  ## TODO - what happens if there is no matching positions?
  gr_peak <- get_peak_centers(dir)
  dt_dist <- sapply(gr_positions, function(gr_anno) {
    ids <- GenomicRanges::nearest(gr_peak, gr_anno)
    pos_nearest <- start(gr_anno)[ids]
    rel_distance <- ifelse(strand(gr_peak) == "+", 1, -1) * (pos_nearest - start(gr_peak))
    return(rel_distance)  
  }, simplify=FALSE) %>% as.data.table

  dt_dist[, class := gr_peak$class]
  write_csv(dt_dist, file.path(dir, "positions.csv"))  
}

## Run for all
pblapply(dir_list, generate_position_csv)

## --------------------------------------------
