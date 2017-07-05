#'---
#' title: Create design_matrix
#'---
flog.info("reading info")
grpeaks <- readRDS("data/encode/eclip/processed/peak_center-gene_mapping_with_sequence.rds")

dtpeaks <- grpeaks %>% as.data.frame %>% as.data.table

flog.info("create directory")
all_rbp <- dtpeaks[, unique(rbp)]
CORE_DIR <- "data/encode/eclip/processed/design_matrix/"
sapply(file.path(CORE_DIR, c("train", "valid", "test", "meta_info")), dir.create, showWarnings = FALSE)
## make directory if it doesn't exist

flog.info("run for all rbp's")
pblapply(all_rbp, function(crbp) {
  ## create info
  info <- list(metainfo = c("seqnames", "strand", "gene_id", "gene_name", "rbp"),
               sequence = "seq",
               response = "binding_site",
               distance_cols = c("TSS_distance", "polya_distance"),
               path_train = file.path(CORE_DIR, "train", paste0(crbp, ".csv")),
               path_valid = file.path(CORE_DIR, "valid", paste0(crbp, ".csv")),
               path_test = file.path(CORE_DIR, "test", paste0(crbp, ".csv")),
               config_file = file.path(CORE_DIR, "meta_info", paste0(crbp, ".json")),
               rbp = crbp
               )

  ## split by chromosomes using roughly a 60/20/20 split
  all_seqnames <- dtpeaks[, unique(seqnames)] %>% as.character
  info$valid_chr <- paste0("chr", c(1, 3))
  info$test_chr <- paste0("chr", c(2, 4, 6, 8, 10))
  info$train_chr <- setdiff(all_seqnames, c(info$test_chr, info$valid_chr))

  ## subset the table wrt rows and columns. note the fraction of total
  incl_features <- c(info$metainfo, info$sequence, info$response, info$distance_cols)
  dt <- dtpeaks[rbp == crbp]
  n_all <- nrow(dt)
  dt_train <- dt[seqnames %in% info$train_chr][, incl_features, with = F]
  dt_valid <- dt[seqnames %in% info$valid_chr][, incl_features, with = F]
  dt_test  <- dt[seqnames %in% info$test_chr][, incl_features, with = F]
  info$train_frac <- nrow(dt_train) / n_all
  info$valid_frac <- nrow(dt_valid) / n_all
  info$test_frac <- nrow(dt_test) / n_all

  ## save the tables
  write_csv(dt_train, info$path_train)
  write_csv(dt_valid, info$path_valid)
  write_csv(dt_test, info$path_test)
  ## save info
  write(jsonlite::toJSON(info, pretty = TRUE, auto_unbox = TRUE), info$config_file)
  return(TRUE)
})

flog.info("done")
