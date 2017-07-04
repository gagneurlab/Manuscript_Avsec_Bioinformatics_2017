#!/usr/bin/env Rscript
#'---
#' title: Generate fasta sequences
#'---
##' usage script.R rbp_name

## generate fasta sequences for all proteins
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) stop("One argument required")
rbp_name <- args[1]

eclip_load_modelling_data(rbp_name)     #get all required data

types <- c("train", "test", "valid")
ispos <- c(TRUE, FALSE)

vars <- expand.grid(types = types, ispos = ispos, stringsAsFactors = FALSE)

get_DNAStringSet <- function(dt) {
  dss <- dt[[info$sequence]] %>% DNAStringSet
  names(dss) <- dt[, paste0(gene_id,"|", TSS_distance,"|", polya_distance)]
  return(dss)
}

to_fasta <- function(sset, train_type, is_positive, info) {
  posneg <- ifelse(is_positive, "pos", "neg")

  basepath <- file.path("data/encode/eclip/processed/design_matrix/",
                        train_type, "fasta")
  stopifnot(dirname(basepath) == dirname(info[[paste0("path_", train_type)]]))
  dir.create(basepath, showWarnings = FALSE)

  endpath <- file.path(basepath, paste0(info$rbp,"_", posneg,".fa"))
  writeXStringSet(x = sset,
                  filepath = endpath,
                  format = "fasta",
		  append = FALSE,
                  width = 202)
  flog.info(paste0("Successfully wrote fasta to: ", endpath))
}


for (i in 1:nrow(vars)) {
  train_type <- vars[i,"types"]
  ispos <- vars[i,"ispos"]
  get(paste0("dt_", train_type))[binding_site == ispos] %>% 
    get_DNAStringSet %>%
    to_fasta(train_type, ispos, info)
}
