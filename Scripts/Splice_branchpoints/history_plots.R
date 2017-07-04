#'---
#' title: Concise position effect inference
#' author: Žiga Avsec
#' wb:
#'  input: ["data/Concise/Splice_branchpoints/trials/train_history/gam_vs_relu.csv"]
#'---
#'
#'
#' 
#' ---------------
opts_chunk$set(echo=FALSE, cache=F, results = 'hide', messages = FALSE)
options(width=140)
#+ 
library(seqLogo)
library(gglogo)
source_all(dir = "Scripts/Concise/Splice_branchpoints/plot_functions")
#' 
ret <- gam_vs_relu_loss_curves_plot(min_epoch = 5, use_metrics = c("accuracy", "f1", "loss"))
legend_off = theme(legend.position = "none")
with(ret, plot_grid(get_legend(deep+ theme(legend.position = "top") + guides(alpha="none")),
                    plot_grid(deep  + legend_off,
                              shallow  + legend_off),
                    ncol = 1,
                    rel_heights = c(0.05, 1)))


#+ fig.height = 3
ret <- gam_vs_relu_loss_curves_plot()
yl <- ylab("binary crossentropy loss (validation)")
legend_off = theme(legend.position = "none")
with(ret, plot_grid(get_legend(deep+ theme(legend.position = "top") + guides(alpha="none")),
                    plot_grid(deep + ylim(c(NA, 0.13)) + legend_off + yl,
                              shallow + ylim(c(NA, 0.18)) + legend_off + yl),
                    ncol = 1,
                    rel_heights = c(0.05, 1)))
a = 1
