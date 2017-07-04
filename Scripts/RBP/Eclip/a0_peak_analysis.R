#'---
#' title: a0_peak_analysis.R
#'---

flog.info("reading info")
grpeaks <- readRDS("data/encode/eclip/processed/peak_center-gene_mapping_with_sequence.rds")

dtpeaks <- grpeaks %>% as.data.frame %>% as.data.table

## make directory if it doesn't exist

##' 
##' ### Overall summary statistics for true binding sites
##' 
##' in median, we have 2500 peaks per rbp:
dtpeaks_t <- dtpeaks[binding_site == TRUE]
dtpeaks_t[,.N , by = rbp][,N] %>% summary

crbp <- "UPF1"
##'
##' ### Restrict to just one rbp: `r crbp`
dtpeaks_t<- dtpeaks_t[rbp %in% crbp]

##' number of peak positions
qplot(as.factor(N), data = dtpeaks_t[, .N, by = gene_name][order(N)])

##' median distance between positions ~ 200 
dtpeaks_t[, mad(polya_distance)[.N > 1], by = gene_name][, V1[V1 < 2000]] %>% hist(breaks = 100)


## ##'
## ##' Restrict to just a single gene
## cgene <- "VBP1"

## ## number of peaks per gene:
## dtpeaks[gene_name %in% cgene]

##--------------------------------------------

## they just negatively sampled the peak positions
## - my negative peaks are those who don't overlap with any experimentally found peaks?
