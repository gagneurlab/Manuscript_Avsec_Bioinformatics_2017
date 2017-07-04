#'---
#' title: Concise position effect inference
#' author: Žiga Avsec
#' wb:
#'  input: ["data/Concise/Splice_branchpoints/processed/branchpointer/train/branchpoint_df_HCN.csv",
#'          "data/Concise/Splice_branchpoints/interpret/position/concise_shallow.csv"]
#'  output : ["data/plots/Concise/Splice_branchpoints/concise_shallow_position_part1.pdf",
#'            "data/plots/Concise/Splice_branchpoints/concise_shallow_position_part1.png",
#'            "data/plots/Concise/Splice_branchpoints/concise_shallow_position_part2.pdf",
#'            "data/plots/Concise/Splice_branchpoints/concise_shallow_position_part2.png"]
#'---
#'
#' ## Conclusions
#'
#' - Fitting works perfectly
#'
#' ## TODO
#' 
#' - [ ] have a common scale for color and size (and display just one legend)
#'
#' 
#' ---------------
opts_chunk$set(echo=FALSE, cache=F, results = 'hide', messages = FALSE)
options(width=140)
#+ 
library(cowplot)
library(seqLogo)
library(gglogo)
plt_dir <- "./data/plots/Concise/Splice_branchpoints/"
dir.create(plt_dir, showWarnings = FALSE)
#+ cache=TRUE
source_all(dir = "Scripts/Concise/Splice_branchpoints/plot_functions")

#+ fig.width=8, fig.height=8

## TODO - add legend at some point
pl_l <- get_position_plots()
plt1 <- with(pl_l,
             plot_grid(pl_d2, m_pl_d2,
                       pl_d1, m_pl_d1,
                       pl_ppts, m_pl_ppts,
                       pl_pptl, m_pl_pptl,
                       ncol=2))
plt1
#+ fig.width=8, fig.height=10
plt2 <- with(pl_l, plot_grid(pl_c1, m_pl_c1,
                             pl_c2, m_pl_c2,
                             pl_c3, m_pl_c3,
                             pl_c4, m_pl_c4,
                             pl_c5, m_pl_c5,
                             ncol=2))
plt2

#+
save_plot(file.path(plt_dir, "concise_shallow_position_part1.pdf"), plt1,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 4,
          base_height = 2,
          base_aspect_ratio = 2,
          )

save_plot(file.path(plt_dir, "concise_shallow_position_part1.png"), plt1,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 4,
          base_height = 2,
          base_aspect_ratio = 2,
          )

save_plot(file.path(plt_dir, "concise_shallow_position_part2.pdf"), plt2,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 5,
          base_height = 2,
          base_aspect_ratio = 2,
          )

save_plot(file.path(plt_dir, "concise_shallow_position_part2.png"), plt2,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 5,
          base_height = 2,
          base_aspect_ratio = 2,
          ## base_height=3.5,
          )

a <- 1
