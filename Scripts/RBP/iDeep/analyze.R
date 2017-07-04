#'---
#' title: make the ideep positional plot
#' author: avsec
#' wb:
#'   input=[]
#'---
## analyze the results
library(cowplot)
source("Scripts/RBP/iDeep/functions.R")

dir_list <- list.dirs("~/github/iDeep/datasets/clip/") %>% grep("training_sample_0$", ., value=TRUE)

dirname_list <- sapply(dir_list, function(x) basename(gsub("/5000/training_sample_0$", "", x)))
names(dirname_list) <- NULL

dt_dist_list <- lapply(dir_list, function(cdir) fread(file.path(cdir, "positions.csv")))
names(dt_dist_list) <- dirname_list
dt_dist <- rbindlist(dt_dist_list, idcol="rbp")

dtm_dist <- dt_dist %>% melt(id.vars = c("rbp", "class"))
dtm_dist[, logval := sign(value) * log10(abs(value) + 1)]
dtm_dist[, rbp := dir2rbp(rbp)]
#+ fig.width=12, fig.height=10
qplot(x = logval, color = class, data=dtm_dist, geom="density") + facet_grid(rbp~variable) + 
  theme(strip.text.y = element_text(angle = 0))
ggsave("data/plots/RBP/iDeep/positional_preference.pdf")




## NSUN analysis

## TODO - create the torkler plot?
dt_nsun <- get_peak_centers(dir_list[30])

gr <- get_human_anno(version="hg19", use_mito=TRUE)
gr_genes <- gr[gr$type == "gene"]

gr_all <- get_human_anno(version="hg19", only_protein_coding=FALSE, use_mito=TRUE)
gr_genes_all <- gr_all[gr_all$type == "gene"]

dir_list <- dir_list[grepl("training", dir_list)]

#' 
#' ## Only protein coding
qpl <- tork_plot_iclip(dir_list, gr_genes)
qpl + facet_wrap(~rbp, nrow = 3) + ggtitle("protein-coding genes")

#' 
##' ## all genes
qpl_all <- tork_plot_iclip(dir_list, gr_genes_all)
qpl_all + facet_wrap(~rbp, nrow = 3) + ggtitle("all genes")



## TODO - fix the bug regarding gene length
## TODO - create the mapping between dir and rbp name


