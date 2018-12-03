## annotation relevant functions
##
## get_<organism>_<type>
##
## <type>:
## _anno = read in the annotation GTF file
## _fasta = read in the fasta file
## _conservation = read in the conservation score
##
## author: avsec
##
## TODO: add version argument for being back compatible
## TODO: add attribute "version" to the returned object

library(magrittr)
library(Biostrings)
library(rtracklayer)
library(magrittr)


## human
get_human_anno <- function(only_regular_chr = TRUE, use_mito=FALSE, only_protein_coding = TRUE, version = "hg19") {
  if (version == "hg38") {
    ## use GRCh38.p7 - corresponds to gencode v25
    rds_gff_file <- "fasta/gencode.v25.annotation.gtf.rds"
  } else {
    stop("Invalid version")
  }

  anno <- readRDS(rds_gff_file)

  if (isTRUE(use_mito)) {
    mito_chr <- c("chrM")
  } else {
    mito_chr <- c()
  }
  if (only_regular_chr) {
    anno <- anno[seqnames(anno) %in% c(paste0("chr", 1:22), "chrX", "chrY", mito_chr)]
    seqlevels(anno) <- seqlevels(anno)[seqlevels(anno) %in% unique(seqnames(anno))]
  }
  
  if (only_protein_coding) {
    anno <- anno[anno$gene_type == "protein_coding"]
  }
  
  return(anno)
}

get_human_fasta <- function(version = "hg19") {
  if (version == "hg38"){
    ## use GRCh38.p7 - corresponds to gencode v25
    fafile <- "fasta/GRCh38.p7.genome.fa"
  } else {
    stop("Invalid version")
  }

  set_my <- getSeq(open(FaFile(fafile)))
  ## discard second entry if it exists 
  names(set_my) <- tstrsplit(names(set_my), " ")[[1]]
  return(set_my)
}

regular_human_chromosome <- function() {
  return(c(paste0("chr", 1:22), "chrX", "chrY"))
}
