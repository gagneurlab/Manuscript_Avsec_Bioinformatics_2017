#'---
#' title:Script for matching metagenes with eclip peak centers
#'--- 

##' requires: eclip data with metadata, ...
##' produces "data/eclip/processed/peak_center-gene_mapping.rds"

## config:
## ----
## number of times more negative examples to generate than positive sequences
print("start")
options(warn=1)
N_TIMES_NEGATIVE <- 4
SEQUENCE_LENGTH <- 101
## ----

message("read anno")
anno <- get_human_anno(version = "hg38")
fa <- get_human_fasta(version = "hg38")
genes <- anno[anno$type == "gene"]
message("read raw data")
grlist <- read_encode_eclip_bed()
mdata <- read_encode_eclip_meta()

cell_lines <- mdata[, unique(sample)]

##' ## Merge peaks by samples OR ovelap peaks from the same protein
message("find consensus peaks")
grlist_by_protein <- sapply(mdata[, unique(protein)],
                            function(prot)
                              sapply(cell_lines, function(cl) {
                                grlist[mdata[prot == protein & sample == cl, file]] 
                              }, simplify = FALSE), simplify = FALSE)

# intersect for replicates, append for different cell-lines
protein_peak_overlaps <- pblapply(grlist_by_protein, . %>% lapply(intersect_center)
                                  %>% gr_rbindlist(idcol="sample")) %>% gr_rbindlist(idcol="protein")
message("save protein_peak_overlaps")
saveRDS(protein_peak_overlaps,
        "data/eclip/processed/protein_peak_overlaps.rds")
## --------------------------------------------

protein_peak_overlaps <- readRDS("data/eclip/processed/protein_peak_overlaps.rds")
## save protein_peak_overlaps

## for each protein - get the peak positions wrt each gene
ol <- findOverlaps(genes, protein_peak_overlaps)
genes_ol <- genes[queryHits(ol)]
peaks_ol <- protein_peak_overlaps[subjectHits(ol)]

genes_ol$rbp_peak_center <- start(peaks_ol)
genes_ol$rbp <- peaks_ol$protein
mcols(genes_ol) <- mcols(genes_ol)[,c("gene_id", "gene_name",
                                      "rbp_peak_center", "rbp")]

genes_ol$binding_site <- TRUE

## stats

length(protein_peak_overlaps)

## not mapped to any of the genes:
1 - uniqueN(subjectHits(ol)) / length(protein_peak_overlaps)
## [1] 0.09943839
## OLD - [1] 0.1013787
## Number of peaks with a single match
data.table(pid = subjectHits(ol))[, .N, by =pid][, mean(N ==  1)]
## Number of genes per peak
length(subjectHits(ol))/uniqueN(subjectHits(ol))

## append generate negative training examples
## 0. choose rbp
## 1. choose gene
## 2. uniformly sample position wrt gene
## 3. discard positions overlapping with any previous positive hits (100 bp distance)

all_rbps <- genes_ol$rbp %>% unique
crbp <- "EFTUD2"
message("generate negative examples")
genes_ol_all <- lapply(all_rbps, function(crbp) {
  cgenes_ol <- genes_ol[genes_ol$rbp == crbp]

  ## simulate random positions
  n_examples <- length(cgenes_ol) * N_TIMES_NEGATIVE
  rel_positions <- runif(n_examples, 0, 1)
  genes <- sample(cgenes_ol$gene_id, n_examples, replace = TRUE)

  ## put the positions into the context of genes
  cgenes_unique<- copy(cgenes_ol)
  cgenes_unique$rbp_peak_center <- NULL
  cgenes_unique <- cgenes_unique %>% unique %>%
    as.data.frame %>% as.data.table %>%
    setkey("gene_id")
  cgenes_negative <- cgenes_unique[genes] %>%   #gene sampling
    .[, rbp_peak_center := round((end-start) *
                                   rel_positions) + start] %>% #pos sampling
    .[, binding_site := FALSE] %>%
    setkey(NULL)

  ## Find the set of negative sequences overlapping true binding sites
  ## ---------
  ## granges of negative sequences. width = 1
  cgenes_negative_gr <- makeGRangesFromDataFrame(cgenes_negative,
                                                 keep.extra.columns = T)
  start(cgenes_negative_gr) <- cgenes_negative_gr$rbp_peak_center
  end(cgenes_negative_gr) <- cgenes_negative_gr$rbp_peak_center
  ## granges of positive sequences. width = 101
  cgenes_positive_gr <- copy(cgenes_ol)
  start(cgenes_positive_gr) <- cgenes_positive_gr$rbp_peak_center
  end(cgenes_positive_gr) <- cgenes_positive_gr$rbp_peak_center
  cgenes_positive_gr <- resize(cgenes_positive_gr, fix = "center", width = SEQUENCE_LENGTH)
  ## get examples to discard
  ol <- findOverlaps(cgenes_negative_gr, cgenes_positive_gr)
  discard <- queryHits(ol) %>% unique

  ## discard overlapping regions
  cgenes_negative <- cgenes_negative[- discard]
  cgenes_negative_gr <- makeGRangesFromDataFrame(cgenes_negative,
                                                 keep.extra.columns = T)
  ## merge the two regions
  ## -----
  ## match colnames
  colorder <- mcols(cgenes_ol) %>% colnames
  mcols(cgenes_negative_gr) <- mcols(cgenes_negative_gr)[, colorder]
  ## merge

  cgenes_final <- c(cgenes_ol,cgenes_negative_gr)
  return(cgenes_final)
}) %>% GRangesList %>% unlist
message("done")
## -----------------------------------
message("define boundaries and save the file")
dt_genes <- genes_ol_all %>% as.data.frame %>% as.data.table
stopifnot(nrow(dt_genes[strand == "*"]) == 0)
dt_genes[, TSS_distance := ifelse(strand == '+',
                                  rbp_peak_center - start,
                                  end - rbp_peak_center)]
dt_genes[, polya_distance := ifelse(strand == '+',
                                    rbp_peak_center - end,
                                    start - rbp_peak_center)]
## write_csv(dt_genes, "data/eclip/processed/peak_center-gene_mapping.csv")
saveRDS(dt_genes, "data/eclip/processed/peak_center-gene_mapping.rds")
message("done")
