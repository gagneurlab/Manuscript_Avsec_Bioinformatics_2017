#'---
#' title: Config
#'---
## LocusZoom
## Lancet
## npg
library(RColorBrewer)
library(cowplot)
library(ggsci)
library(dtplyr)
library("scales") #show_col

mypal <- brewer.pal(10, "Paired")

## mypal = pal_locuszoom()(7)
col_concise_shallow <- mypal[7]
col_concise_deep <- mypal[8]
## col_relu_shallow <- mypal[3]
## col_relu_deep <- mypal[4]
col_relu_shallow <- mypal[9]
col_relu_deep <- mypal[10]
col_glmnet <- mypal[1]
col_branchpointer <- mypal[2]

## aliases
legend_off = theme(legend.position = "none")

save_plot_mul <- function(basename, formats, ...) {
  library(cowplot)
  for (format in formats) {
    save_plot(paste0(basename, ".", format), ...)
  }
}

white_facets <- theme(strip.background = element_rect(colour="white",
                                                      fill="white")) # remove tabs from facets

source_all("Scripts/Figures/func/")

## mypal = pal_locuszoom( alpha = 0.7)(10)
## mypal
## library("scales")
## show_col(mypal)
