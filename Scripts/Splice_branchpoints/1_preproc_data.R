#'---
#' title: Convert long to wide format
#' author: Žiga Avsec
#' wb:
#'   input: ["data/Concise/Splice_branchpoints/processed/branchpointer/train/branchpoint_df_HCN_w_id.csv",
#'          "data/Concise/Splice_branchpoints/processed/branchpointer/test/branchpoint_df_HCN_w_id.csv"]
#'   output: ["data/Concise/Splice_branchpoints/processed/branchpointer/test/wide_data.csv",
#'            "data/Concise/Splice_branchpoints/processed/branchpointer/train/wide_data.csv"]
#'---
GENOME_FA <- "~/github/splice_branchpoints/data/inputs/Homo_sapiens.GRCh37.67.dna.toplevel.fa"
flog.info("read the genome")
genome <- readDNAStringSet(GENOME_FA)
flog.info("done")
## fix chromosome names
chrs <- genome %>% names %>% tstrsplit(" ") %>% .[[1]]
chrs <- paste0("chr", chrs)
exclude_chr <- "chrY"
keep <- chrs %in% regular_human_chromosome() & !chrs %in% exclude_chr
names(genome) <- chrs
genome <- genome[keep]

dt_train <- fread("data/Concise/Splice_branchpoints/processed/branchpointer/train/branchpoint_df_HCN_w_id.csv")
dt_test <- fread("data/Concise/Splice_branchpoints/processed/branchpointer/test/branchpoint_df_HCN_w_id.csv")

## main function
long_to_wide <- function(dt, genome) {
  dt_comp <- complete(dt, exon_id.2, dist.2) %>% as.data.table
  ## strand and chr are the same witin exon
  dt_comp[, chr := na.omit(chr)[1], by = exon_id.2]
  dt_comp[, strand := na.omit(strand)[1], by = exon_id.2]

  dt_comp[, set] %>% table(useNA = "always")
  dt_comp

  ## branchpoint distribution per 
  dt_comp[, sum(set == "HC", na.rm=TRUE), by = exon_id.2][, table(V1)] %>%
    barplot(xlab = "Number of branchpoints per intron region")

  ## merge the sequence columns into one
  seq_cols <- grep("^seq_", colnames(dt_comp), val=T)
  dt_comp[, seq := do.call(paste0, c(.SD[, seq_cols, with = F]))]
  dt_comp[, (seq_cols) := NULL]



  ## add the position and restrict start/end to the interval
  ## all have the first position
  ## sometimes the first position is missing

  dt_comp[, which(is.na(set)), by = exon_id.2][, table(V1)] %>% barplot
  dt_comp[dist.2==18, is.na(set), by = exon_id.2][, table(V1)]
  dt_comp[dist.2==44, is.na(set), by = exon_id.2][, table(V1)]


  ## fill in the positions
  dt_comp[, start := if (strand == "-") {fill_incr_na(start)} else {fill_incr_na(start, -1)}, by = .(exon_id.2, strand)]
  dt_comp[, end := if (strand == "-") {fill_incr_na(end)} else {fill_incr_na(end, -1)}, by = .(exon_id.2, strand)]
  dt_comp[, position := end]
  dt_comp[, start := min(end, na.rm = TRUE), by = exon_id.2]
  dt_comp[, end := max(end, na.rm = TRUE), by = exon_id.2]
  ## fix the coordinates
  ## expand by few bases
  dt_comp[, start := start -5]
  dt_comp[, end := end +5]


  dt_comp[nchar(seq)> 12 , seq := gsub("A", "", seq)]
  stopifnot(all(dt_comp[, nchar(seq) == 11]))

  ## which columns still have NA's
  dt_comp %>% sapply(. %>% is.na %>% any)

  ## get ranges
  gr <- dt_comp[, .(exon_id.2, chr, start = start, end=end, strand)] %>%
    unique %>% as("GRanges")
  gr_seq <- easy_Views(genome, gr, reverse_seq=TRUE)

  ## add wide ranges
  gr_wide <- dt_comp[, .(exon_id.2, chr, start = start-100, end=end+100, strand)] %>%
    unique %>% as("GRanges")
  gr_seq_wide <- easy_Views(genome, gr_wide, reverse_seq=TRUE, colname = "seq_wide")
  gr_seq$seq_wide <- gr_seq_wide$seq_wide

  ## check the sequence width:
  stopifnot(nchar(gr_seq$seq) == 37) # 5 more on each side


  ## convenient visualization tool
  str_pad_lr <- function(seq, n_left=0, n_right=1, pad = " ") {
    len <- nchar(seq)[1]
    stopifnot(uniqueN(nchar(seq)) == 1)
    
    pad1 <- str_pad(seq, n_left + nchar(seq), side = "left", pad = pad)
    return(str_pad(pad1, n_right + nchar(pad1), side = "right", pad = pad))
  }

  ## ## train dataset
  ## test_exon_neg <- "ENST00000000412.3.4"
  ## test_exon_pos <- "ENST00000577161.1.19"
  ## ## test dataset
  ## test_exon_neg <- "ENST00000040738.5.13"
  ## test_exon_pos <- "ENST00000226299.4.4"
  ## ## visually check if it is working
  ## test_exon <- test_exon_neg
  ## seq <- dt_comp[exon_id.2 == test_exon][, seq]
  ## seqs <- str_pad_lr(seq, n_right= 0:(length(seq)-1), n_left = (length(seq)-1):0)
  ## extracted <- gr_seq[gr_seq$exon_id.2 == test_exon]$seq
  ## as.data.frame(c(extracted, seqs, extracted))
  ## test_exon <- test_exon_pos
  ## seq <- dt_comp[exon_id.2 == test_exon][, seq]
  ## seqs <- str_pad_lr(seq, n_right= 0:(length(seq)-1), n_left = (length(seq)-1):0)
  ## extracted <- gr_seq[gr_seq$exon_id.2 == test_exon]$seq
  ## as.data.frame(c(extracted, seqs, extracted))

  dtseq <- gr_seq %>% as.data.frame %>% as.data.table
  setnames(dtseq, "seqnames", "chr")
  dtseq[, width := NULL]
  dt_comp <- merge(dt_comp[, -"seq", with = F],
                         dtseq, all.x = TRUE, by = c("chr", "start", "end", "strand", "exon_id.2"), sort=FALSE)

  ## convert set to 1, 0
  dt_comp[, set := as.integer(set == "HC")]
  setnames(dt_comp, "set", "is_branchpoint")
  flog.info("reverse order of dist.2")
  dt_comp <- dt_comp[order(exon_id.2, -dist.2)]
  ## long to wide
  flog.info("long to wide")
  stopifnot(dt_comp[,.N, by = .(chr, start, end, strand, exon_id.2, seq, seq_wide)][, all(N==27)])
  dt_comp[, position := NULL]
  dt_wide <- dt_comp[, lapply(.SD, function(x) paste(x, collapse=",")), by = .(chr, start, end, strand, exon_id.2, seq, seq_wide)]
  return(dt_wide)
}

flog.info("long-> wide for train")
dt_train_wide <- long_to_wide(dt_train, genome)
flog.info("long-> wide for test")
dt_test_wide <- long_to_wide(dt_test, genome)

flog.info("save the data")
write_csv(dt_train_wide, "data/Concise/Splice_branchpoints/processed/branchpointer/train/wide_data.csv")
write_csv(dt_test_wide, "data/Concise/Splice_branchpoints/processed/branchpointer/test/wide_data.csv")
flog.info("done")
## DONE - extract the whole sequences - increase the range by the window size
## DONE - compare them with the old ones (pad + align them)
## DONE - remove the sequence and concatenate everything
## DONE - write the same pipeline for the test set
## MAYBE - fill NA's with either extremely high or low values or using fill()
