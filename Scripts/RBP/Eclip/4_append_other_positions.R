#'---
#' title: Append other positional features
#' author: Avsec
#'---

rbp <- "UPF1"

source("Scripts/RBP/iDeep/functions.R")

anno <- get_human_anno(version = "hg38")

gr_positions <- get_position_landmarks(anno)
gene_start <- gr_positions$gene_start
gene_end <- gr_positions$gene_end

strand2num <- function(x) {
  ifelse(x=="+", +1, -1)
}

get_peaks <- function(dt, gr_positions) {
  gene_start <- gr_positions$gene_start
  gene_end <- gr_positions$gene_end

  dt_gene_start <- gene_start %>% as.data.frame %>% as.data.table %>% .[, .(gene_start=start, seqnames, strand, gene_id)]
  dt_gene_end <- gene_end %>% as.data.frame %>% as.data.table %>% .[, .(gene_end=start, seqnames, strand, gene_id)]
  dt <- merge(dt, dt_gene_start, by =c("seqnames", "strand", "gene_id"), all.x=TRUE, sort=FALSE)
  dt <- merge(dt, dt_gene_end, by =c("seqnames", "strand", "gene_id"), all.x=TRUE, sort=FALSE)
  dt[, pos := gene_start + strand2num(strand) * TSS_distance]
  dt[, pos_alt := gene_end + strand2num(strand) * polya_distance]
  stopifnot(dt[, all(pos == pos_alt)])
  dt[, pos_alt := NULL]
  dt[, start := pos]
  dt[, end := pos]
  return(dt %>% as("GRanges"))
}

splits <- c("train", "valid", "test")

dir_list <- list.files(paste0("data/eclip/processed/design_matrix/", splits), full.names = T) %>%
  grep("\\.csv$", ., value=TRUE)
dir_list <- dir_list[!grepl("extended", dir_list)]

## find the nearest feature point

file <- dir_list[1]
generate_additional_position_csv <- function(file) {
  print(file)
  dt_init <- fread(file)

  dt_init[, id := 1:.N]
  gr_peak <- get_peaks(dt_init, gr_positions)
  for (feat in names(gr_positions)) {
    gr_anno <- gr_positions[[feat]]
    dt_anno <- gr_anno %>% as.data.frame %>% as.data.table %>% subset(select=c("seqnames", "start", "end", "strand", "gene_id"))
    dt_peak <- gr_peak %>% as.data.frame %>% as.data.table
    dt_anno <- dt_anno %>% setkeyv(c("seqnames", "strand", "gene_id", "start"))
    dt_peak <- dt_peak %>% setkeyv(c("seqnames", "strand", "gene_id", "start"))
    dt_peak_closest <- dt_anno[dt_peak, roll="nearest"] %>% setkey(NULL) %>% unique
    dt_peak_closest[, pos_nearest := end]
    dt_peak_closest <- dt_peak_closest[order(id)]
    
    rel_distance <- dt_peak_closest[, strand2num(strand) * (pos_nearest - pos)]
    dt_init[, (feat) := rel_distance]
  }
  stopifnot(dt_init[, all(TSS_distance == - gene_start)], 
            dt_init[, all(polya_distance == - gene_end)])
  write_csv(dt_init, gsub("\\.csv$", "_extended.csv", file))
}

pblapply(dir_list, generate_additional_position_csv)
