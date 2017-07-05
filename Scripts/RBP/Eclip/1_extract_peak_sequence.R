#'---
#' title: Extract peak sequence
#'---
## #' wb:
## #'   input: ["data/eclip/processed/peak_center-gene_mapping.rds"]
## #'   output: ["data/eclip/processed/peak_center-gene_mapping_with_sequence.rds"]
## #'---

SEQUENCE_LENGTH <- 101

flog.info("read the rds data")
dt_genes <- readRDS("data/eclip/processed/peak_center-gene_mapping.rds")
flog.info("read fasta")
fa <- get_human_fasta(version = "hg38")

## fa <- get_human_fasta(version = "hg38")

## extract the underlying sequence
flog.info("extract the underlying sequence")
setnames(dt_genes, c("start", "end"), c("gene_start", "gene_end"))
dt_genes[, width := NULL]
grpeaks <- makeGRangesFromDataFrame(dt_genes,
                                    keep.extra.columns = TRUE, start.field = "rbp_peak_center", end.field = "rbp_peak_center")

grpeaks <- resize(grpeaks, width = SEQUENCE_LENGTH, fix = "center")
flog.info("extract sequence")
grpeaks <- easy_Views(fa, grpeaks, reverse_seq = TRUE) #append sequence
flog.info("write to disk")
saveRDS(grpeaks, "data/eclip/processed/peak_center-gene_mapping_with_sequence.rds")
