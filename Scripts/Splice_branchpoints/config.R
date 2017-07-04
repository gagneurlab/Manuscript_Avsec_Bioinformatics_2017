#'---
#' title: Analyze branchpointer performance
#' author: Ziga Avsec
#'---
library(stringr)
library(data.table)
library(caret)
library(ggplot2)
library(cowplot)
library(reshape2)
library(PRROC)

theme_figure <- theme_bw()+ theme(text=element_text(size=10),legend.key.size=unit(0.2, "inches"),
                                  panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.border = element_blank(),
                                  axis.line.x = element_line(), 
                                  axis.line.y = element_line(), 
                                  panel.background = element_rect(colour = "black", size=1, fill=NA),
                                  plot.title = element_text(hjust = 0.5))

nt_cols <- c("#359646","#4D7ABE","#FAA859","#CB3634")
