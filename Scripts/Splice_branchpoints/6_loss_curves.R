#'---
#' title: Gam vs relu
#' author: Žiga Avsec
#' wb:
#'  input: ["data/Concise/Splice_branchpoints/trials/train_history/gam_vs_relu.csv",
#'          "data/Concise/Splice_branchpoints/trials/df/deep_gam.csv"]
#'---
#'
#'
#' 
#' ---------------
opts_chunk$set(echo=FALSE, cache=F, results = 'hide', messages = FALSE)
options(width=140)
#+ 
source_all(dir = "Scripts/Concise/Splice_branchpoints/plot_functions")
library(seqLogo)
library(gglogo)
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
pl_loss_curves <- with(ret, plot_grid(get_legend(deep+ theme(legend.position = "top") +
                                                   guides(alpha="none")),
                                      plot_grid(deep + ylim(c(NA, 0.13)) + legend_off + yl,
                                                shallow + ylim(c(NA, 0.18)) + legend_off + yl),
                                      ncol = 1,
                                      rel_heights = c(0.05, 1)))
pl_loss_curves

dt_eval <- get_trials_df_eval()
pl <- ggplot(dt_eval, aes(x = type, y = eval.auprc)) +
  geom_boxplot(aes(color = type)) +
  geom_jitter(aes(color = type), alpha = 0.1)
pl_auprc_deep <- pl %+% dt_eval[depth == "deep"] + ggtitle("deep") 
pl_auprc_shallow <- pl %+% dt_eval[depth == "shallow"] + ggtitle("shallow")

perf_boxplot <- plot_grid(lgd = get_legend(pl_auprc_deep + legend_top),
                          plot_grid(pl_auprc_deep + legend_off,
                                    pl_auprc_shallow + legend_off,
                                    ncol=2),
                          ncol=1, rel_heights = c(0.05, 1)
                          )
perf_boxplot

