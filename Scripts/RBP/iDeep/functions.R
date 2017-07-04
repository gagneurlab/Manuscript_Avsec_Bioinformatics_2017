#'---
#' Title: functions
#'---
## get the bed file

## bed_file_gz <- file.path(dir, "positions.bedGraph.gz")
## bed_file <- file.path(dir, "positions.bedGraph")
## bed_file_final <- file.path(dir, "positions2.bedGraph")
## if (file.exists(bed_file_gz)) {
##   system(paste0("gunzip ", bed_file_gz))
## }
## if (file.exists(bed_file_final)) {
##   system(paste0("tail -n +2 ", bed_file, " > ", bed_file_final))
## }


## bf <- rtracklayer::import(bed_file_final)
## stopifnot(all(anno[, center ] == start(bf)))

get_peak_centers <- function(dir) {
  library("Biostrings")
  fa <- readDNAStringSet(file.path(dir, "sequences.fa.gz"))
  anno <- names(fa) %>% tstrsplit(",|;") %>% as.data.table
  setnames(anno, c("chr", "strand", "start", "end", "class"))
  anno[, start := as.numeric(start)]
  anno[, end := as.numeric(end)]
  anno[, center := (end + start)/2 + 1]

  anno[, chr := trimws(chr)]
  gr_peaks <- anno %>% as("GRanges")
  gr_peaks <- resize(gr_peaks, width = 1, fix="center") %>% shift(1)
  stopifnot(all(start(gr_peaks) == gr_peaks$center))
  return(gr_peaks)
}

get_position_landmarks <- function(gr) {
  gr_TSS <- gr[gr$type == "transcript"] %>% resize(width=1, fix="start")
  gr_polya <- gr[gr$type == "transcript"] %>% resize(width=1, fix="end")

  gr_exon_intron <- gr[gr$type == "exon"] %>% resize(width=1, fix="end")
  gr_intron_exon <- gr[gr$type == "exon"] %>% resize(width=1, fix="start")

  gr_start_codon <- gr[gr$type == "start_codon"] %>% resize(width=1, fix="start")
  gr_stop_codon <- gr[gr$type == "stop_codon"] %>% resize(width=1, fix="end")

  gr_gene_start <- gr[gr$type == "gene"] %>% resize(width=1, fix="start")
  gr_gene_end <- gr[gr$type == "gene"] %>% resize(width=1, fix="end")
  
  ## remove exons that are also on the transcript ends
  gr_exon_intron <- gr_exon_intron[!gr_exon_intron %in% gr_polya]
  gr_intron_exon <- gr_intron_exon[!gr_intron_exon %in% gr_TSS]
  
  return(list(TSS = gr_TSS,
       polya = gr_polya,
       exon_intron = gr_exon_intron,
       intron_exon = gr_intron_exon,
       start_codon = gr_start_codon,
       stop_codon = gr_stop_codon,
       gene_start = gr_gene_start,
       gene_end = gr_gene_end
       ))
}

dir2rbp <- function(dir) {
  ## TODO - get the file from Amin
  dt_hash <- fread("data/encode/eclip/processed/predictions/iClip_dir_rbp_mapping.csv")
  setkey(dt_hash, "dir")
  return(dt_hash[dir][, rbp])
}

tork_plot_iclip <- function(dirs, gr_genes,
                            bp_bucket_size = 5* 500, gene_bucket_size = 2* 100, saturate_at = 50,
                            gene_length_quantile = .8, line_alpha=0.4,
                            cylab="Genes sorted by length",
                            cxlab="Distance to gene start (nt)") {
  rbps <- dir2rbp(tstrsplit(dirs, "/")[[12]])
  ##' ## Torkler plots
  #+
  library(scales)
  
  dt_all <- sapply(dirs, function(dir) {
    gr_peak <- get_peak_centers(dir) %>% as("GRanges")
    torkler_plot_table(gr_genes, gr_peak,
                       bp_bucket_size = bp_bucket_size, gene_bucket_size = gene_bucket_size,
                       saturate_at = saturate_at,
                       gene_length_quantile = gene_length_quantile)},
    simplify=FALSE) %>% rbindlist(idcol="dir")
  print(dt_all[, tstrsplit(dir, "/")])
  dt_all[, rbp := dir2rbp(tstrsplit(dir, "/")[[12]])]
  ## dt_all[, rbp := fct_relevel(rbp, rbps)]
  dt_all[, dir := fct_relevel(dir, dirs)]

  ggplot(dt_all, aes(position_bucket, gene_bucket)) +
    geom_raster(aes(fill = signal), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    xlab("Position wrt gene start [bp]") +
    ylab("Genes sorted by length") +
    ## facet_grid(~dir) + 
    theme_cowplot() +
    geom_line(aes(gene_end_bucket, gene_bucket), alpha=line_alpha)+
    geom_line(aes(gene_start_bucket, gene_bucket), alpha=line_alpha) +
    xlab(cxlab) + ylab(cylab) +
    legend_off + 
    scale_y_continuous(labels = scales::scientific_format(digits=2)) 
    ## scale_x_continuous(breaks = c(0, 50000), labels = c(0, 50000))
}
