##--------------------------------------------
## required packages

## data + convenience libraries
library(futile.logger)
library(data.table)
## library(reshape2) ## functons like melt etc by Hadley Wickham

## hadleyverse
library(tidyr)
library(forcats) # factor manipulation
library(readr)
library(stringr)


library(Matrix)
library(magrittr)
library(broom)

## plotting libraries
library(ggplot2)
library(ggsignif)
library(ggrepel)
## library(gplots)                 #heatmap.2 function
## library(scales)                 #plot(..., col = alpha(color, 0.5)) function for base r plotting

## Parallel processing
library(parallel)
library(BiocParallel)
library(pbapply)
library(pbmcapply)

## ML libraries
library(caret)
library(glmnet)
library(mgcv)
library(pROC)
library(PRROC)

## Genomic libraries
library(Biostrings)
## library(Rsamtools)		# provides an interface to BAM files
## library(GenomicAlignments)	# Representation and manipulation of short genomic alignment
library(GenomicRanges)		# The GenomicRanges package serves as the foundation for representing genomic
library(rtracklayer)            #read in all the formats
library(GenomicFeatures)
# library(VariantAnnotation)  	# handling vcf/bcf files and some other variant annotations

library(cowplot)
library(ggsci)
library(dtplyr)
##--------------------------------------------
## Global variables

STD_CHROMOSOMES <- paste0("chr", c(1:22, "X","Y"))

## PHhome_ase <- "~/ase-grant"
PHhome <- "."
PHdata="./data"
PHdatao="./data-offline"
##--------------------------------------------
## User defined functions

## IMPORTANT: local functions have to be loaded last, as I could override some function definitions
## get the function_get_helmholtz_file

source_all <- function(dir) {
  sapply(list.files(path = dir, pattern=".*\\.R$", full.names=TRUE), source, .GlobalEnv)
}

## general functions
source_all("Scripts/functions/avsec_utils")

## project specific functions
source_all("./Scripts/functions")
##--------------------------------------------
## Knitr config
library(knitr)
library(rmarkdown)
opts_knit$set(root.dir = getwd())
opts_chunk$set(echo=FALSE, cache=F)
