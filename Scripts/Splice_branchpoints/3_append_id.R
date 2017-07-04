#'---
#' title: Append id
#' author: Žiga Avsec
#' wb:
#'  input: ["data/Concise/Splice_branchpoints/processed/branchpointer/train/branchpoint_df_HCN.csv",
#'           "data/Concise/Splice_branchpoints/processed/branchpointer/test/branchpoint_df_HCN.csv"]
#'  output: ["data/Concise/Splice_branchpoints/processed/branchpointer/train/branchpoint_df_HCN_w_id.csv",
#'          "data/Concise/Splice_branchpoints/processed/branchpointer/test/branchpoint_df_HCN_w_id.csv"]
#'---
dt_raw <- fread("data/Concise/Splice_branchpoints/processed/branchpointer/train/branchpoint_df_HCN.csv")
dt_raw_test <- fread("data/Concise/Splice_branchpoints/processed/branchpointer/test/branchpoint_df_HCN.csv")

dt_raw[order(new_ID)][, is.sorted(dist.1), by = new_ID]
dt_raw[, new_ID] %>% head
dt_raw[, .N, by = new_ID][, table(N)]
## TODO - how to we see if we are from the same batch?

is.sorted(1:10)

process_dt <- function(dt) {
  message("process dt")
  ## TODO - make a function
  dt[, V1 := NULL]
  dt <- separate(dt, "new_ID", into = c("chr", "pos", "class"), sep = "_")
  stopifnot(all(dt[, class == set]))
  dt[, class := NULL]
  dt[, strand := ifelse(grepl("\\+", pos), "+", "-")]
  stopifnot(dt[, all(grepl("\\-", pos) == (strand == "-"))])
  dt <- separate(dt, "pos", into = c("start", "end"))
  dt[, start := as.integer(start)]
  dt[, end := as.integer(end)]
  stopifnot(dt[, all(start == end -1)])  
  return(dt)
}

dt <- process_dt(dt_raw)
dt_test <- process_dt(dt_raw_test)





## TODO - check if things are sorted
## - dist.1
## - dist.2
## - dist.


## TODO - extract the sequences from the fasta file
## TODO - compare them with the refrences


## TODO - merge with 
br_data_path <- "~/github/splice_branchpoints/data/"

chrom=c("chr1","chr2","chr3", "chr4","chr5",
        "chr6", "chr7", "chr8", "chr9","chr10",
        "chr11", "chr12", "chr13", "chr14","chr15",
        "chr16", "chr17", "chr18", "chr19","chr20",
        "chr21","chr22", "chrX","chrY")

dt_all <- lapply(chrom, function(chr)
  fread(file=paste0(br_data_path, "/outputs/branchpoint_df_with_seq_r2_", chr,".csv"))) %>%
  rbindlist(use.names = TRUE)
dt_all
dt_all[, type.1] %>% table


TF=apply(dt_all, 2, function(x) x=="N")
TF_v=apply(TF,1, any)
rm_row =which(TF_v==T)
dt_all <- dt_all[-rm_row]

dt_all[, transcript_id.1 == transcript_id.2] %>% table
dt_all[, gene_id.1 == gene_id.2] %>% all

## TODO - why different
## not all are of length 27
dt_all[, .N, by = exon_id.2][, table(N)] %>% barplot
dt_all[, set] %>% table  ## No LC data
44-18
dt_metainfo  <- dt_all[, .(chr = Chromosome, start = Start,
                           end = Stop, strand = Strand, exon_id.2, ppt_start)]

## append exon_id.2
message("merge the tables")
dt <- merge(dt, dt_metainfo, by = c("chr", "start", "end", "strand", "ppt_start"),
            all.x = TRUE, sort=FALSE)

dt_test <- merge(dt_test, dt_metainfo, by = c("chr", "start", "end", "strand", "ppt_start"),
                 all.x = TRUE, sort=FALSE)
message("done!")
dt[, .N, by = exon_id.2][, table(N)] %>% barplot
## TODO - order it within exon_id.2 by dist.2
## TODO - merge with dummy data (exon_id.2 x 10:...)

## TODO - get the start location(s) for each exon_id.2 (dist.2 should be 18)
##        - for all the distances, check that they match if you elongate them
##           - synthetically create a table with everything except the set and seq_*
##           - compare that they are the same
## TODO - extract the sequence data
## TODO - compare the extracted sequences 

## Final dataset:
## chr, start, end, strand, ppt_run_length

## TODO - this ppt_run_length can be different...
dt[, .(nl = uniqueN(ppt_run_length)), by = exon_id.2][, nl] %>% table %>% barplot

## TODO - use NA's for all the other values like ppt_run_length etc...

## TODO - can we use NA's in keras

## Sort the values
dt <- dt[order(exon_id.2, dist.2)]


## TODO - save the values
write_csv(dt, "data/Concise/Splice_branchpoints/processed/branchpointer/train/branchpoint_df_HCN_w_id.csv")
write_csv(dt_test, "data/Concise/Splice_branchpoints/processed/branchpointer/test/branchpoint_df_HCN_w_id.csv")
