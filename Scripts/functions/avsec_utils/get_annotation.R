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

## cerevisiae
get_cer_anno <- function(only_regular_chr = TRUE, only_protein_coding = TRUE,
                         gtf_file = "/data/ouga01/home/ag_gagneur/share/genomes/sacCer/sacCer3/annotation/Saccharomyces_cerevisiae.R64-1-1.31.gtf") {
  anno <- import(gtf_file)

  sl <- seqlevels(anno)
  sl[sl != "Mito"] <- sl[sl != "Mito"] %>% as.roman %>% as.integer %>% paste0("chr", .)
  sl[sl == "Mito"] <- "chrM"
  seqlevels(anno) <- sl
    
  if (only_regular_chr) {
    anno <- anno[seqnames(anno) != c("chrM")]
  }
  
  if (only_protein_coding) {
    anno <- anno[anno$gene_biotype == "protein_coding"]
  }
  
  return(anno)
}

get_cer_fasta <- function() {
  fasta_dir <- "/data/ouga01/home/ag_gagneur/share/genomes/sacCer/sacCer3/"
  fa_file <- open(FaFile(file.path(fasta_dir, "sacCer3.fa")))
  set_saccer3 <- getSeq(fa_file)

  ## clean up the chromosome names
  sl <- names(set_saccer3)
  sl <- gsub("^chr", "", sl)
  sl[sl != "M"] <- sl[sl != "M"] %>% as.roman %>% as.integer
  sl <- paste0("chr", sl)
  names(set_saccer3) <- sl
  
  return(set_saccer3)
}

get_cer_conservation <- function(bw_path = "/s/genomes/sacCer/sacCer3/conservation/phastCons/sacCer3.phastCons7way.bw") {
  
  wig <- import(bw_path)
  sl <- seqlevels(wig)
  sl <- gsub("^chr", "", sl)
  sl[sl != "M"] <- sl[sl != "M"] %>% as.roman %>% as.integer
  sl <- paste0("chr", sl)
  seqlevels(wig) <- sl
  
  return(wig)
}

## TODO - refactor
normalize_cer_chr <- function(sl) {
  sl <- gsub("^chr", "", sl)
  sl[sl != "M"] <- sl[sl != "M"] %>% as.roman %>% as.integer
  sl <- paste0("chr", sl)  
  return(sl)
}

## pombe
get_pombe_anno <- function(only_regular_chr = TRUE, only_protein_coding = TRUE) {
  anno_file <- "/s/genomes/sacPombe/ASM294v2.31/Schizosaccharomyces_pombe.ASM294v2.31.gtf"
  anno <- import(anno_file)
  sl <- seqlevels(anno)
  suppressWarnings(sl_NA <- as.roman(sl))
  sl[!is.na(sl_NA)] <- as.numeric(sl_NA[!is.na(sl_NA)])
  sl <- paste0("chr", sl)
  seqlevels(anno) <- sl
  
  if (only_regular_chr) {
    anno <- anno[seqnames(anno) %in% paste0("chr", 1:3)]
  }
  if (only_protein_coding) {
    anno <- anno[anno$gene_biotype == "protein_coding"]
  }

  return(anno)
}

get_pombe_fasta <- function() {
  file <- open(FaFile("/s/genomes/sacPombe/ASM294v2.31/Schizosaccharomyces_pombe.ASM294v2.31.dna.genome.fa"))
  set <- getSeq(file)

  sl <- names(set)
  sl <- sl %>% strsplit(" ") %>% lapply(function(x) x[[1]][1]) %>% unlist
  suppressWarnings(sl_NA <- as.roman(sl))
  sl[!is.na(sl_NA)] <- as.numeric(sl_NA[!is.na(sl_NA)])
  sl <- paste0("chr", sl)
  names(set) <- sl

  return(set)
}

get_pombe_conservation <- function() {
  file <-  "/s/genomes/sacPombe/ASM294v2.31/conservation/phyloP_like/Sp.cons.wig"
  wig <- import(file)

  seqlevels(wig) <- seqlevels(wig) %>% as.roman %>% as.numeric %>% paste0("chr", .)
  return(wig)
}

## human
get_human_anno <- function(only_regular_chr = TRUE, use_mito=FALSE, only_protein_coding = TRUE, version = "hg19") {

  ## gff_file <- "/data/ouga01/home/ag_gagneur/share/genomes/human/hg19/gencode/gencode.v24lift37.annotation.gff3.gz"
  ## anno <- import(gff_file)
  if (version == "hg19") {
    rds_gff_file <- "/s/genomes/human/hg19/gencode/gencode.v24lift37.annotation.rds"
  } else if (version == "hg38") {
    ## use GRCh38.p7 - corresponds to gencode v25
    rds_gff_file <- "/s/genomes/human/hg38/GRCh38.p7/gencode.v25.annotation.gtf.rds"
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
  if (version == "hg19") {
    fafile <- "/s/genomes/human/hg19/fasta/hg19.fa"
  } else if (version == "hg38"){
    ## use GRCh38.p7 - corresponds to gencode v25
    fafile <- "/s/genomes/human/hg38/GRCh38.p7/GRCh38.p7.genome.fa"
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
