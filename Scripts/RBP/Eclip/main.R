#'---
#' title: Methods performance
#' author: Žiga Avsec
#' wb:
#'   input: ["data/encode/eclip/processed/predictions/UPF1/kmer_glmnet-no_position.csv"]
#'---
##' ## Goal
##'
##' - generate AUC-comparing plot
##' 
##' ## Conclusions
##'
##' 
##' ## TODO's
##'
##' - fix the axis descriptions
##' - distinguish positional vs non-positional features
##' - show the performances in a boxplot:
##'   - use bootstrap to compare the performances
##' 
##' ---------------------------------------------------------------
##' 
##+ set_workdir, echo=F
library(knitr)
library(rmarkdown)
library(cowplot)
library(gridExtra)
opts_knit$set(root.dir = getwd())
opts_chunk$set(echo=FALSE, cache=F, results = 'hide', messages = FALSE)
options(width=140)
## 
plt_dir <- "data/plots/RBP/Eclip/"
dir.create(plt_dir, showWarnings = FALSE, recursive = TRUE)
source_all("Scripts/RBP/Eclip/plot_functions")
source("Scripts/RBP/Eclip/config.R")
##+

#' 
#' ## Get the markdown tables
roc_curves <- get_roc_curves()
pr_curves <- get_pr_curves()
roc_pl <- plot_roc(roc_curves)
pr_pl <- plot_pr(pr_curves)
#' ## Plot
#+ fig.width = 12, fig.height=6


dt_pred <- get_predictions()

lgd <- get_legend(roc_pl + legend_top)

plt <- plot_grid(roc_pl + legend_off,
                 pr_pl + legend_off,
                 lgd,
                 labels = c("a", "b"),
                 ncol=1,
                 rel_heights = c(1,1,0.1))
plt

dt_pred_auc <- dt_pred[, .(auc = bootstrap_metric(y_true, y_pred, metric_auc,
                                                  n_bootstrap=50, return_summary=FALSE)),
                       by = .(rbp, method)]

dt_pred_auprc <- dt_pred[, .(auprc = bootstrap_metric(y_true, y_pred, metric_auprc,
                                                      n_bootstrap=50, return_summary=FALSE)),
                       by = .(rbp, method)]


#+ fig.width = 12
ggplot(dt_pred_auc, aes(x = method, color = method, y = auc)) +
  geom_boxplot() +
  geom_jitter(alpha=0.2) +
  facet_grid(.~rbp) +
  legend_top + 
  ylim(c(0.85, NA)) +
  tilt_xlab

#+ fig.width = 12
ggplot(dt_pred_auprc, aes(x = method, color = method, y = auprc)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.2) +
  facet_grid(.~rbp) +
  legend_top + 
  tilt_xlab

ggplot(dt_pred[, .(auprc = metric_auprc(y_true, y_pred)),
               by = .(rbp, method)],
       aes(x = rbp, y = auprc, color = method)) +
  legend_top + 
  geom_point(position=position_dodge(0.5))



save_plot(file.path(plt_dir, "roc_pr.pdf"), plt,
          ncol = 6, # we're saving a grid plot of 2 columns
          nrow = 2,
          base_aspect_ratio = .8,
          base_height=3,
          # each individual subplot should have an aspect ratio of 1.3
          )

save_plot(file.path(plt_dir, "roc_pr.png"), plt,
          ncol = 6, # we're saving a grid plot of 2 columns
          nrow = 2,
          base_aspect_ratio = .8,
          base_height=3,
          )

y <- 1
