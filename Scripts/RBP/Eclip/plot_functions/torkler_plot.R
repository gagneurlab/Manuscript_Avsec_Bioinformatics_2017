#'---
#' title: Helper functions for torkler plots
#' author: Avsec
#'---
## load the data
get_human_anno <- function(only_regular_chr = TRUE, only_protein_coding = TRUE, version = "hg19") {
  ## TODO updated file

  ## anno <- import(gff_file)
  if (version == "hg19") {
    rds_gff_file <- "/s/genomes/human/hg19/gencode/gencode.v24lift37.annotation.rds"
  } else if (version == "hg38") {
    ## use GRCh38.p7 - corresponds to gencode v25
    rds_gff_file <- "data/annotation/gencode.v25.annotation.gtf.rds"
  } else {
    stop("Invalid version")
  }

  anno <- readRDS(rds_gff_file)

  if (only_regular_chr) {
    anno <- anno[seqnames(anno) %in% c(paste0("chr", 1:22), "chrX", "chrY")]
    seqlevels(anno) <- seqlevels(anno)[seqlevels(anno) %in% unique(seqnames(anno))]
  }
  
  if (only_protein_coding) {
    anno <- anno[anno$gene_type == "protein_coding"]
  }
  
  return(anno)
}

torkler_plot_anno <- get_human_anno(version = "hg38")
torkler_plot_genes_gr <- torkler_plot_anno[torkler_plot_anno$type == "gene"]
torkler_plot_protein_peak_overlaps <- readRDS("data/encode/eclip/processed/protein_peak_overlaps.rds")

top_genes_tork_plot <- function(proteins_of_choice, 
                                bp_bucket_size = 5* 500, gene_bucket_size = 2* 100, saturate_at = 50,
                                gene_length_quantile = .8, line_alpha=0.4,
                                cylab="Genes sorted by length",
                                cxlab="Distance to gene start (nt)") {
  ##' ## Torkler plots
  protein_peak_overlaps <- readRDS("data/encode/eclip/processed/protein_peak_overlaps.rds")
  #+
  library(scales)
  
  dt_all <- sapply(proteins_of_choice, function(rbp) {
    torkler_plot_table(torkler_plot_genes_gr, torkler_plot_protein_peak_overlaps[
      torkler_plot_protein_peak_overlaps$protein == rbp],
      bp_bucket_size = bp_bucket_size, gene_bucket_size = gene_bucket_size,
      saturate_at = saturate_at,
      gene_length_quantile = gene_length_quantile)}, simplify=FALSE) %>% rbindlist(idcol="rbp")
  dt_all[, rbp := fct_relevel(rbp, proteins_of_choice)]

  ggplot(dt_all, aes(position_bucket, gene_bucket)) +
    geom_raster(aes(fill = signal), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    xlab("Position wrt gene start [bp]") +
    ylab("Genes sorted by length") +
    facet_grid(~rbp) + 
    theme_cowplot() +
    geom_line(aes(gene_end_bucket, gene_bucket), alpha=line_alpha)+
    geom_line(aes(gene_start_bucket, gene_bucket), alpha=line_alpha) +
    xlab(cxlab) + ylab(cylab) +
    legend_off + 
    scale_y_continuous(labels = scales::scientific_format(digits=2)) + 
    scale_x_continuous(breaks = c(0, 50000), labels = c(0, 50000))
}


## --------------------------------------------
