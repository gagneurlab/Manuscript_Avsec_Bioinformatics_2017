##--------------------------------------------
## required packages

## data + convenience libraries
library(futile.logger)
library(tictoc)
library(data.table)

## hadleyverse
library(tidyr)
library(forcats)
library(readr)
library(stringr)


library(Matrix)
library(magrittr)
library(broom)

## plotting libraries
library(ggplot2)
library(ggsignif)
library(ggrepel)
library(scales)

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
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(Rsamtools)

library(cowplot)
library(ggsci)
library(dtplyr)
##--------------------------------------------
## Global variables
STD_CHROMOSOMES <- paste0("chr", c(1:22, "X","Y"))

PHhome <- "."
PHdata="./data"
PHdatao="./data-offline"
##--------------------------------------------
## User defined functions
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
