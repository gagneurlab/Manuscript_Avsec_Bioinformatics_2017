#'---
#' title: Config
#'---
## LocusZoom
## Lancet
## npg
library(cowplot)
library(ggsci)
library(dtplyr)
mypal = pal_locuszoom()(7)
col_concise_shallow <- mypal[4]
col_concise_deep <- mypal[5]
col_relu_shallow <- mypal[3]
col_relu_deep <- mypal[1]
col_branchpointer <- mypal[2]
col_glmnet <- mypal[7]

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

source_all("Scripts/Concise/Paper/func")

## mypal = pal_locuszoom( alpha = 0.7)(10)
## mypal
## library("scales")
## show_col(mypal)
