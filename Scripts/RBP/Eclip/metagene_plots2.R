#'---
#' title: Meta-gene plots for encode 2
#' author: Ziga Avsec
#' wb:
#'   input: ["data/encode/eclip/processed/peak_center-gene_mapping.rds",
#'           "data/encode/eclip/processed/protein_peak_overlaps.rds"]
#'---

## #' output: 
## #'   html_document:
## #'     pandoc_args: ["+RTS", "-K64m","-RTS"]
## #'     toc: yes
## #'     css: ../baderd_html.css
## #'---


#' 
##' ## Goal
##'
##' - Select 6 RBP's showing prefernce for:
##'    - poly-a
##'    - TSS
##'    - no-preference
##' 
##' ## Conclusions
##'
##' - The following RBP's got selected:
##'   - UPF1
##'   - PUM2
##'   - DDX3X
##'   - NKRF
##'   - TARDBP
##'   - SUGP2
##' 
##' 
##' 
##' ## TODO's
##' 
##' * DONE : Generate negative peak samples
##' * CANCELLED : check with torkler plot your negative peaks
##' * TODO : Build a model for predicting parclip peaks
##'     - with positional preference
##'     - without positional preference
##' * TODO : Create the sequence-based torkler plot and identify proteins looking more to position
##'     - also check where the gain is the largest
##' * DONE : Implement CONCISE to work with such data and run it
##' 
##'
##' ---------------------------------------------------------------

## ##+ include_default_knitr_values, child='../mertes_html_default_include.Rmd'
##' 
##+ set_workdir, echo=F
library(knitr)
library(rmarkdown)
opts_knit$set(root.dir = getwd(), messages=FALSE)
opts_chunk$set(echo=FALSE, cache=F)
options(width=140)
library(cowplot)
#'
#' ## Explore
##' 
##' Read in the Encode eclip data. Obtained from [this](https://www.encodeproject.org/matrix/?type=Experiment&status=released&assay_slims=RNA+binding&assay_title=eCLIP&target.investigated_as=RNA+binding+protein&assembly=GRCh38&files.file_type=bed+narrowPeak) link.

source_all("Scripts/RBP/Eclip/plot_functions")
## grlist <- read_encode_eclip_bed()
## mdata <- read_encode_eclip_meta()
## anno <- get_human_anno(version = "hg38")
## fa <- get_human_fasta(version = "hg38")
## genes <- anno[anno$type == "gene"]

dt_genes <- readRDS("data/encode/eclip/processed/peak_center-gene_mapping.rds")

#' 
#' ## Choosing 6 RBP's
#'
#' ### Make a statistical test comparing binding vs non-binding sites
dt_t <- dt_genes[, c(t.test(TSS_distance[binding_site == TRUE]/width[binding_site == TRUE],
                            TSS_distance[binding_site == FALSE]/width[binding_site == FALSE]
                            )[c("statistic", "p.value")], "N"= .N), by = rbp] %>%
  .[order(-abs(statistic),-N)]
#' 
#' ### Poly-A preference
#+ echo=TRUE, results='show'
dt_t[statistic > 0] %>% head(10)
#'
#' ### TSS preference
#+ echo=TRUE, results='show'
dt_t[statistic < 0] %>% head(10)
#'
#' ### Neutral elements
#+ echo=TRUE, results='show'
dt_t[p.value > 0.05]
#'
#' ### RBP choice:
#+ echo=TRUE, results='show'
proteins_of_choice <- c(dt_t[N>10000 & statistic > 0, rbp[1:2]],
                        dt_t[N>10000 & statistic < 0, rbp[1:2]],
                        dt_t[N>10000 & p.value > 0.05][order(-p.value)][, rbp[1:2]])
print(proteins_of_choice)

##' 
##' ### Site locations
##'
##'  Lower panel = true binding sites,
##'
##' 
##' #### Natural scale
##' 
#+ fig.width = 13, fig.height = 3
some_genes <- proteins_of_choice
library(forcats)
dt_genes[, rbp := fct_relevel(rbp, proteins_of_choice)]
qplot(polya_distance, data = dt_genes[rbp %in% some_genes & polya_distance > -2000],
      bins = 30) + tilt_xlab + 
  facet_grid(facets = binding_site~rbp)
qplot(TSS_distance, data = dt_genes[rbp %in% some_genes & TSS_distance < 2000],
      bins = 30) + tilt_xlab +
  facet_grid(facets = binding_site~rbp)
#'
#' #### Log-scale
#+ fig.width = 13, fig.height = 3
qplot(- polya_distance + 1, data = dt_genes[rbp %in% some_genes],
      bins = 30, log="x") + tilt_xlab + 
  facet_grid(facets = binding_site~rbp, scales="free")
qplot(TSS_distance + 1, data = dt_genes[rbp %in% some_genes],
      bins = 30, log ="x") + tilt_xlab + 
  facet_grid(facets = binding_site~rbp, scales="free")


##' 
##' ## Torkler plots for selected RBP's
#+ fig.width=12, fig.height=3
top_genes_tork_plot(proteins_of_choice)
a <- 1

## --------------------------------------------
#'
#' ## All RBP's
## fraction
#+ fig.height=15, fig.width=10
ggplot(dt_genes,aes(x = TSS_distance/ width, color = binding_site)) +
  geom_density() +
  facet_wrap(~rbp, ncol = 8) +
  xlab("Relative gene position") + 
  ## ggtitle("Poly-a distance") +
  legend_top + 
  tilt_xlab


dt_genes_wp <- merge(dt_genes, dt_t, by="rbp")
dt_genes_wp[, rbp := fct_reorder(rbp, -statistic)]
ggplot(dt_genes_wp ,aes(x = TSS_distance/ width, color = binding_site)) +
  geom_density() +
  facet_wrap(~rbp, ncol = 8) +
  xlab("Relative gene position") + 
  ## ggtitle("Poly-a distance") +
  legend_top + 
  tilt_xlab

#' RBP's are reordered by their t-test test statistic.

## Density plots
## truncated version
## ggplot(dt_genes[polya_distance > -2000],aes(x = polya_distance, color = binding_site)) +
##   geom_density() +
##   facet_wrap(~rbp, ncol = 8) +
##   ggtitle("Poly-a distance") +
##   tilt_xlab

## ## polya-log
## ggplot(dt_genes,aes(x = -polya_distance + 1, color = binding_site)) +
##   scale_x_log10() + 
##   geom_density() +
##   facet_wrap(~rbp, ncol = 8) +
##   ggtitle("Poly-a distance") +
##   tilt_xlab

## ## tss - log
## ggplot(dt_genes,aes(x = TSS_distance, color = binding_site)) +
##   scale_x_log10() + 
##   geom_density() +
##   facet_wrap(~rbp, ncol = 8) +
##   ggtitle("Poly-a distance") +
##   tilt_xlab



## TODO - choose the RBP panels
## Maybe color them:
##  - TSS = green
##  - poly-a = blue
##  - None = Red

##--------------------------------------------
#'
#' ## All new features

dir_list <- list.files("/s/project/deepcis/encode/eclip/processed/design_matrix/train", full.names = T) %>%
  grep("extended.csv$", ., value=TRUE)
dt_all <- lapply(dir_list, fread) %>% rbindlist

rbp_subset <- c(proteins_of_choice, c("LARP4", "SMNDC1", "MTPAP", "FTO", "SUB1", "FXR1",
                                     "SF3B4", "YBX3", "LIN28B",
                                     "SFPQ", "HNRNPUL1", "SLTM", "RBFOX2", "HNRNPM", "HNRNPK")) # no effect
## rbp_subset <- c(proteins_of_choice, sample(dt_all[, unique(rbp)], 10))
features <- c("TSS", "polya", "exon_intron", "intron_exon", "start_codon", "stop_codon", 
              "gene_start", "gene_end")

features <- c("gene_start", "TSS", "start_codon",
              "exon_intron", "intron_exon",
              "stop_codon", "polya", "gene_end")

dtm_all <- melt(dt_all, id.vars=c("rbp", "binding_site"), measure.vars=features)
dtm_all_pl <- dtm_all[rbp %in% rbp_subset]
dtm_all_pl[, rbp := fct_relevel(rbp, rbp_subset)]
dtm_all_pl[, variable := fct_relevel(variable, features)]
#+ fig.height=15, fig.width=10

## TODO - provide this as 
ggplot(dtm_all_pl, aes(x = - sign(value) * log10(abs(value)), color = binding_site)) +
  geom_density() +
  facet_grid(rbp~variable) + 
  xlab("Signed log10 distance to landmark") + 
  theme(strip.text.y = element_text(angle=0))


