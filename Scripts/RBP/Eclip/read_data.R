#'---
#' title: read the data
#'---
## 
## 
datadir <- "./data/encode/eclip"
## unzip all the gz files
bbfiles <- list.files(datafir, pattern = 'bed$', full.names = TRUE)
ids <- sapply(bbfiles, basename) %>% gsub("\\.bed", "", .)

##' import narrow-peak files
##'
##' - [Peak format link](https://genome.ucsc.edu/FAQ/FAQformat.html#format12)
##'
##' Test if import worked
gr <- import_narrowPeak(bbfiles[1])
gr

## name <-> File accession

mdata <- fread(file.path(datadir, "metadata.tsv"))
setnames(mdata, make_names(colnames(mdata)) %>% gsub("\\(s\\)", "", .))
mdata[Assembly == "GRCh38"]
mdata %>% head(2)

##' all downloaded files are in the table
mdata[, mean(File_accession %in% ids)]
mdata[, mean(ids %in% File_accession)]


##' Column description:
##'
##' - Experiment_target - used protein
##' all are human

##' RBP's used:
mdata[, Experiment_target] %>% tstrsplit("-") %>% lapply(table)
##' All were done in human.
##' 
##' In total, they used 112 different targets
mdata[, Experiment_target] %>% tstrsplit("-") %>% .[[1]] %>% table %>% length
##'
##' With 5-6 experiments per RBP per sample.
mdata[, Experiment_target] %>% tstrsplit("-") %>% .[[1]] %>% table %>% mean

##' Data of interest:
##' For each sample & each target we have 2 experiements
mdata[Assembly == "GRCh38"][, .N, by = .(Biosample_term_name, Experiment_target)]
mdata[Assembly == "GRCh38"][, .N, by = .(Biosample_term_name, Experiment_target)][, table(N)]

##'
##' Where is the difference
mdata[Assembly == "GRCh38"] %>% head

##' Columns where they differ:
diff_cols <- mdata[Assembly == "GRCh38"][1:2][, mult_value_col(.SD)]

##' these are 2 biological replicates
mdata[Assembly == "GRCh38"][1:2, diff_cols, with = F] %>% as.data.frame

##' ## Check biological replicas
##' 
##' Briefly check the 2 biological replicas
##' 
bbfiles_sub <- grep("ENCFF078UAE|ENCFF746WOI", bbfiles, value = T)
grlist <- lapply(bbfiles_sub, import_narrowPeak)
names(grlist) <- sapply(bbfiles_sub, basename) %>% gsub("\\.bed", "", .)

a <- grlist[[1]]
b <- grlist[[2]]
ol <- findOverlaps(a,b)
length(a)
ol %>% queryHits %>% uniqueN

length(b)
ol %>% subjectHits %>% uniqueN

## roughly 20% of peaks are in good agreement
