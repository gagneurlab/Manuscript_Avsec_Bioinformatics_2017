#'---
#' title: Meta-gene plots for encode
#' author: Ziga Avsec
#' wb:
#'   input: ["data/encode/eclip/processed/peak_center-gene_mapping.rds",
#'           "data/encode/eclip/processed/protein_peak_overlaps.rds"]
#'---

## #' output: 
## #'   html_document:
## #'     pandoc_args: ["+RTS", "-K64m","-RTS"]
## #'     toc: yes
## #'     css: ../baderd_html.css
## #'---


#' 
##' ## Goal
##'
##' - See if there is any positional preference for the peaks wrt start/stop positions
##' 
##'
##' ## Conclusions
##'
##' - some RBP's exhibit quite strong positional preference
##' b
##' 
##' 
##' ## TODO's
##' 
##' * DONE : Generate negative peak samples
##' * CANCELLED : check with torkler plot your negative peaks
##' * TODO : Build a model for predicting parclip peaks
##'     - with positional preference
##'     - without positional preference
##' * TODO : Create the sequence-based torkler plot and identify proteins looking more to position
##'     - also check where the gain is the largest
##' * DONE : Implement CONCISE to work with such data and run it
##' 
##'
##' ---------------------------------------------------------------

## ##+ include_default_knitr_values, child='../mertes_html_default_include.Rmd'
##' 
##+ set_workdir, echo=F
library(knitr)
library(rmarkdown)
opts_knit$set(root.dir = getwd())
opts_chunk$set(echo=FALSE, cache=F)
options(width=140)
#'
#' ## Explore
##' 
##' Read in the Encode eclip data. Obtained from [this](https://www.encodeproject.org/matrix/?type=Experiment&status=released&assay_slims=RNA+binding&assay_title=eCLIP&target.investigated_as=RNA+binding+protein&assembly=GRCh38&files.file_type=bed+narrowPeak) link.

## grlist <- read_encode_eclip_bed()
## mdata <- read_encode_eclip_meta()
## anno <- get_human_anno(version = "hg38")
## fa <- get_human_fasta(version = "hg38")
## genes <- anno[anno$type == "gene"]

dt_genes <- readRDS("data/encode/eclip/processed/peak_center-gene_mapping.rds")
dt_genes_t <- dt_genes[binding_site == TRUE]
dt_genes_f <- dt_genes[binding_site == FALSE]

##' 
##' ### Site locations
##'
##'  Lower panel = true binding sites,
##'
##' 
##' #### Natural scale
##' 
#+ fig.width = 13, fig.height = 5
some_genes <- dt_genes_t[, rbp] %>% unique %>% .[1:10]
qplot(polya_distance, data = dt_genes[rbp %in% some_genes & polya_distance > -2000],
      bins = 30) + tilt_xlab + 
  facet_grid(facets = binding_site~rbp)
qplot(TSS_distance, data = dt_genes[rbp %in% some_genes & TSS_distance < 2000],
      bins = 30) + tilt_xlab +
  facet_grid(facets = binding_site~rbp)
#'
#' #### Log-scale
#+ fig.width = 13, fig.height = 5
qplot(- polya_distance + 1, data = dt_genes[rbp %in% some_genes],
      bins = 30, log="x") + tilt_xlab + 
  facet_grid(facets = binding_site~rbp, scales="free")
qplot(TSS_distance + 1, data = dt_genes[rbp %in% some_genes],
      bins = 30, log ="x") + tilt_xlab + 
  facet_grid(facets = binding_site~rbp, scales="free")


###################################################################################################
## 
## grlist <- read_encode_eclip_bed()
## mdata <- read_encode_eclip_meta()
## fa <- get_human_fasta(version = "hg38")
anno <- get_human_anno(version = "hg38")
genes_gr<- anno[anno$type == "gene"]

##' 
##' ## Torkler plots
protein_peak_overlaps <- readRDS("data/encode/eclip/processed/protein_peak_overlaps.rds")

#+
signal_gr <- protein_peak_overlaps[protein_peak_overlaps$protein == "GEMIN5"]

torkler_plot(genes_gr, protein_peak_overlaps[protein_peak_overlaps$protein == "GEMIN5"],
             bp_bucket_size = 500, gene_bucket_size = 100, saturate_at = 10,
             gene_length_quantile = .8) +
  ggtitle("GEMIN5")

torkler_plot(genes_gr, protein_peak_overlaps[protein_peak_overlaps$protein == "SUB1"],
             bp_bucket_size = 500, gene_bucket_size = 100, saturate_at = 10,
             gene_length_quantile = .8) +
  ggtitle("SUB1")


torkler_plot(genes_gr, protein_peak_overlaps[protein_peak_overlaps$protein == "UPF1"],
             bp_bucket_size = 500, gene_bucket_size = 100, saturate_at = 10,
             gene_length_quantile = .8) +
  ggtitle("UPF1")

## we can safely assume that RBP binding must happen on the mRNA itself

##' Is this specific positional placement due to motif distribution or due to other positional effects we don't know?
##'
##' If later, can we model it and learn it?
