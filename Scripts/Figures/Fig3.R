#'---
#' title: Create figure 2
#' author: Ziga Avsec
#' wb:
#'  input: ["data/Concise/Splice_branchpoints/processed/pr_roc/roc_curves.csv",
#'          "data/Concise/Splice_branchpoints/processed/pr_roc/pr_curves.csv",
#'          "data/Concise/Splice_branchpoints/processed/pr_roc/bootstrap_auc.csv",
#'          "data/Concise/Splice_branchpoints/processed/pr_roc/bootstrap_auprc.csv"]
#'  output: ["data/plots/Concise/Paper/fig3_roc_pr.pdf",
#'           "data/plots/Concise/Paper/fig3_roc_pr.png"]
#'---
#'
#' 
opts_chunk$set(echo=FALSE, cache=F, results = 'hide', messages = FALSE)
options(width=140)
#+ 
plt_dir <- "~/projects-work/spline_trans/plots/subplots/fig3/"
dir.create(plt_dir, showWarnings = FALSE)
library(cowplot)
library(forcats)
source("Scripts/Concise/Paper/config.R")
# Get the data
load_dir <- "data/Concise/Splice_branchpoints/processed/pr_roc/"
roc_curves <- fread(file.path(load_dir, "roc_curves.csv"))
pr_curves <- fread(file.path(load_dir, "pr_curves.csv"))
#' 
roc_curves <- roc_curves[!grepl("relu", method) ]
pr_curves <- pr_curves[!grepl("relu", method) ]

library(gridExtra)
roc_curves[, method_perf := fct_reorder(method_perf, -auROC)]
pr_curves[, method_perf := fct_reorder(method_perf, -auPR)]
## set methods order by performance

roc_curves <- roc_curves[method != "mgcv_gam"]
pr_curves <- pr_curves[method != "mgcv_gam"]

## Colors
cols <- c(col_concise_deep, col_branchpointer,
          col_concise_shallow, col_glmnet)
#+
## Method, AUC
plt_ROC=ggplot(roc_curves, aes(x = FPR,y = TPR, col = method_perf)) + 
  geom_line() + 
  ## scale_color_manual(values=nt_cols) +
  legend_top_ver +
  xlab("False positive rate") + 
  ylab("True positive rate") + 
  ## ggtitle(label="ROC") +
  scale_color_manual(values=cols) +
  guides(color=guide_legend(title="Method, AUC")) + 
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(colour = "black", size=.8, linetype=1, inherit.blank = TRUE),
    legend.justification = c(1, 0), legend.position = c(1, 0), 
    )

#' 
#' ### Precision figures
pr_curves <- pr_curves[!(Recall == 0 & Precision == 0)] # Weird first point

plt_PR <- (plt_ROC+ aes(x = Recall,y = Precision, col = method_perf) +
             theme(legend.justification = c(0, 0), legend.position = c(0, 0)) ) %+% pr_curves

plt <- plot_grid(plt_ROC, plt_PR, ncol=1)
plt


save_plot_mul(file.path(plt_dir, "roc_pr"), c("pdf", "png"), plt,
              ncol = 1, # we're saving a grid plot of 2 columns
              nrow = 2,
              base_height=3.5
          # each individual subplot should have an aspect ratio of 1.3
              )

y <- 1

## --------------------------------------------
## boxplots
dtp_b_auc <- fread(file.path(load_dir, "bootstrap_auc.csv"))
dtp_b_auprc <- fread(file.path(load_dir, "bootstrap_auprc.csv"))

level_order <- c("glmnet",
                 ## "mgcv_gam",
                 ## "concise_shallow_relu",
                 "concise_shallow",
                 "branchpointer",
                 ## "concise_deep_relu",
                 "concise_deep")
## filter
dtp_b_auc <- dtp_b_auc[method %in% level_order]
dtp_b_auprc <- dtp_b_auprc[method %in% level_order]


## update color- legend
cols <- c(col_glmnet,
          ## col_relu_shallow,
          col_concise_shallow,
          col_branchpointer,
          ## col_relu_deep,
          col_concise_deep)
          ## mypal[7])

dtp_b_auc[, method := fct_relevel(method, level_order)]
dtp_b_auprc[, method := fct_relevel(method, level_order)]

## append _gam to concise_shallow

dtp_b_auc %>% setnames("method", "Method")
dtp_b_auprc %>% setnames("method", "Method")

q_auc <- ggplot(dtp_b_auc, aes(x = Method, y = auc)) +
  geom_boxplot(aes(color=Method)) + geom_jitter(aes(color=Method), alpha=0.1, size=0.5) +
  scale_color_manual(labels=c("glmnet", "NN w/ ST shallow",
                              "branchpointer", "NN w/ ST deep"),
                     values=cols) + 
  geom_signif(comparisons = list(## c("concise_shallow_relu", "concise_shallow"),
                                 ## c("concise_deep_relu", "concise_deep"),
                                 c("branchpointer", "concise_deep")
                                 ),
              hjust = 0.8,
              map_signif_level = TRUE,
              step_increase = 0.07) +
  ylab("Area under ROC curve") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
q_auprc <- (q_auc + aes(y = auprc) + ylab("Area under Precision-Recall curve")) %+% dtp_b_auprc

## put method on the top
lgd <- get_legend(q_auc +  legend_top)

plt_boxplot <- plot_grid(plot_grid(q_auc + legend_off, q_auprc + legend_off, align="v"),
                       lgd, ncol=1, rel_heights=c(1, 0.05))
plt_boxplot

save_plot_mul(file.path(plt_dir, "roc_pr_boxplot"), c("pdf", "png"), plt_boxplot,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 1,
          base_height=3.5,
          base_aspect_ratio = 0.82
          # each individual subplot should have an aspect ratio of 1.3
          )
## add spline functions
## ------------------------------------------------------------------
## Plot 1
## requirements
source_all(dir = "Scripts/Concise/Splice_branchpoints/plot_functions")
col_scheme_deep <- scale_colour_manual(values = c(col_concise_deep, col_relu_deep))
col_scheme_shallow <- scale_colour_manual(values = c(col_concise_shallow,
                                                     col_relu_shallow))
## 

dtt <- tidy_position()

## restrict regions
dtt <- dtt[!(variable == "canon_hit1" & position > 75)]
dtt <- dtt[!(grepl("canon_hit", variable) & position > 150)]
dtt <- dtt[!(grepl("canon_hit", variable) & position > 150)]
dtt <- dtt[!(variable == "dist1" & position > 3000)]
dtt <- dtt[!grepl("relu", method)]
## compute outliers
dtt[, outlier := p > 0.2 & method == "measured"]
dtt[outlier == TRUE, p := 0.2]

dtt

plm <- ggplot(dtt[method == "measured" & primary],
              aes(position, p, alpha = N, color = outlier)) +
  ## geom_hline(yintercept = 0.05, lty="dashed", alpha = 0.2) + 
  geom_point() +
  scale_colour_manual(values=c("black", "red")) + 
  facet_grid(~variable_name, scales = "free") +
  legend_off +
  xlab(NULL) +
  ylab("Marginal effect") +
  theme(strip.background = element_rect(colour="white", fill="white")) # remove tabs from facets

pli <- ggplot(dtt[method != "measured" & primary],
              aes(position, p, color = method)) +
  ## geom_hline(yintercept = 0.5, lty="dashed", alpha = 0.2)+
  geom_line() +
  facet_grid(~variable_name, scales = "free") +
  col_scheme_shallow + 
  legend_off +
  theme(strip.text.x = element_blank(),
        strip.background = element_blank()) +
  xlab("Distance") +
  ylab("Inferred effect") 

## main figure - position
plt_pos_effect <- plot_grid(plm, pli, ncol=1, align="v")

## supplementary figure - position
plt_pos_effect_supp <- plot_grid(plm %+% dtt[method == "measured" & !primary],
                                 pli %+% dtt[method != "measured" & !primary],
                                 ncol=1, align="v")


save_plot(file.path(plt_dir, "pos_effects.pdf"), plt_pos_effect,
          ncol=2, nrow=1, base_height=3.5, base_aspect_ratio=0.7,
          dpi=600)

save_plot_mul(file.path(plt_dir, "fig3_sup"), c("pdf", "png", "eps"), plt_pos_effect_supp,
          ncol=7, nrow=1, base_height=3.5, base_aspect_ratio=0.7,
          dpi=600)

## save_plot(file.path(plt_dir, "pos_effects.png"), plt_pos_effect,
##           ncol=4, nrow=1, base_height=3.5, base_aspect_ratio=0.7,
##           dpi=600)
