## 
##' Read the eclip encode narrowPeak data from encode:
##' <https://www.encodeproject.org/matrix/?type=Experiment&status=released&assay_slims=RNA+binding&assay_title=eCLIP&target.investigated_as=RNA+binding+protein&assembly=GRCh38&files.file_type=bed+narrowPeak>
##'
##' @return List of GRanges
##' @author Žiga Avsec
read_encode_eclip_bed <- function() {
  datadir <- "./data/eclip/raw"
  files <- read_encode_eclip_meta()[, File_accession]
  stopifnot(!any(duplicated(files)))
  bfiles <- file.path(datadir, paste0(files, ".bed"))

  grlist <- if (requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pblapply(bfiles, import_narrowPeak)
  } else {
    lapply(bfiles, import_narrowPeak)
  }
  names(grlist) <- files
  return(grlist)
}

##' Get encode eclip meta-data table.
##'
##' @param only_grch38 Restrict the metadata table only to GRCh38 experiments?
##' @return data.table with metadata about eclip experiments. See also \Code{read_encode_eclip_bed}.
##' @author Žiga Avsec
read_encode_eclip_meta <- function(only_grch38 = TRUE) {
  datadir <- "./data/eclip/raw"
  suppressWarnings(mdata <- fread(file.path(datadir, "metadata.tsv")))
  setnames(mdata, make_names(colnames(mdata)) %>% gsub("\\(s\\)", "", .))
  if (isTRUE(only_grch38)) {
    mdata <- mdata[Assembly == "GRCh38"]
  }
  ## add protein 
  mdata[, c("protein", "protein_Organism") := tstrsplit(Experiment_target, "-")]
  mdata[, sample := Biosample_term_name]
  mdata[, file := File_accession]
  return(mdata)
}
