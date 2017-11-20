#'---
#' title: Create figure 2
#' author: Ziga Avsec
#' wb:
#'  input: ["data/Splice_branchpoints/processed/pr_roc/roc_curves.csv",
#'          "data/Splice_branchpoints/processed/pr_roc/pr_curves.csv",
#'          "data/Splice_branchpoints/processed/pr_roc/bootstrap_auc.csv",
#'          "data/Splice_branchpoints/processed/pr_roc/bootstrap_auprc.csv"]
#'  output: ["data/plots/fig3_roc_pr.pdf",
#'           "data/plots/fig3_roc_pr.png"]
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
source("Scripts/Figures/config.R")
# Get the data
load_dir <- "data/Splice_branchpoints/processed/pr_roc/"
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


cols <- c(col_glmnet,
          ## col_relu_shallow,
          col_concise_shallow,
          col_branchpointer,
          ## col_relu_deep,
          col_concise_deep)
          ## mypal[7])

dtp_b_auc[, method := fct_relevel(method, level_order)]
dtp_b_auprc[, method := fct_relevel(method, level_order)]


dtp_b_auc %>% setnames("method", "Method")
dtp_b_auprc %>% setnames("method", "Method")

## TODO - use smaller point-size
q_auc <- ggplot(dtp_b_auc, aes(x = Method, y = auc)) +
  geom_boxplot(aes(color=Method)) + geom_jitter(aes(color=Method), alpha=0.1, size=0.5) +
  scale_color_manual(labels=c("glmnet", "NN w/ ST shallow",
                              "branchpointer", "NN w/ ST deep"),
                     values=cols) + 
  ## TODO - turn off x-labels
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

## TODO - put method on the top
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
## ------------------------------------------------------------------
## Plot 1
## requirements
source_all(dir = "Scripts/Splice_branchpoints/plot_functions")
col_scheme_deep <- scale_colour_manual(values = c(col_concise_deep, col_relu_deep))
col_scheme_shallow <- scale_colour_manual(values = c(col_concise_shallow,
                                                     col_relu_shallow))
## 

dtt <- tidy_position()

## TODO - update this part -
## - add position
## - add motif

## --------------------------------------------
## Motif plots
basedir = "data/Splice_branchpoints/interpret/sequence/shallow_spline_trans"
df <- lapply(list.files(basedir, full.names=TRUE), function(f) {
  dt <- fread(f, header=FALSE)
  dt[, method := gsub("\\.txt", "", gsub("pwm_", "", basename(f)))]
  return(dt)
}) %>% rbindlist
DNA <- c("A", "C", "G", "T")
setnames(df, c(DNA, "method"))

library(ggseqlogo)
pwm_list = list("Branchpoint-centered PWM (Data)" = t(as.matrix(df[method == "data"][, DNA, with = F])),
                "Inferred PWM" = t(as.matrix(df[method == "model_transformed"][, DNA, with = F])),
                "Filter weights" = t(as.matrix(df[method == "model_raw"][, DNA, with = F])))

## TODO - get background probabilty
background <- c(0.2053,  0.254 ,  0.2001,  0.3407)
## TODO - create one model_raw_pwm
mpwm <- exp(pwm_list[["Filter weights"]]) * background
## normalize
pwm_list[["Model PWM (Inferred)"]] <- t(t(mpwm) / colSums(mpwm))

pl_raw <- ggseqlogo(pwm_list[3], method='custom', seq_type='dna') + ylab("Weight size")
pl_raw_pwm <- ggseqlogo(pwm_list[4], seq_type='dna')

pl2 <- ggseqlogo(pwm_list[c(4, 1)], seq_type='dna')
pl2

pl <- ggseqlogo(pwm_list[c(1)])
## pl2 <- ggseqlogo(pwm_list[c(2)], method="probability")
## pl1 <- ggseqlogo(pwm_list[c(1)], method="probability")
## pl_out <- plot_grid(pl_raw, pl2, pl1, nrow=1)
pl_out <- cowplot::plot_grid(pl_raw_pwm, pl, nrow=1, rel_widths = c(1,1))
pl_out

save_plot(file.path(plt_dir, "motif.pdf"), pl_out,
          ncol = 2, base_height=1.5, base_aspect_ratio=2)

save_plot(file.path(plt_dir, "motif2.pdf"), pl2,
          ncol = 2, base_height=1.6, base_aspect_ratio=1.6)

ggseqlogo(t(as.matrix(df[method == "data"][, DNA, with = F])))
## --------------------------------------------
basedir = "data/Splice_branchpoints/interpret/positions/shallow_spline_trans"
dtpos_all<- fread(file.path(basedir, "all_positions.csv"))
setnames(dtpos_all, "outlier", "Outlier")
dtpos_all[, primary := feature %in% c("dist2", "ppt_start", "canon_hit1", "canon_hit2")]
dtpos_all[,V1 := NULL]
dtpos_all[,V2 := NULL]
dtpos_all[,V3 := NULL]

## TODO - add variable name 
dtpos_all[, feature_name := DIST_FEATURES[feature]]
dtpos_all[type =="measured", Frac := N/sum(N), by = feature]

dtpos_all[, feature_name := fct_relevel(feature_name, DIST_FEATURES[c("dist2", "ppt_start")])]

library(plyr)

plm <- ggplot(dtpos_all[primary == TRUE],
       aes(x = x, y=y, alpha = Frac)) + 
  geom_point(aes(group=1, color = "Data"),
             data= function(x)  x[type == "measured" & !Outlier]) + 
  geom_point(#aes(group=1, color = "Outlier"),
             color = "red",
             ## guide="none",
             guide=FALSE,
             data= function(x)  x[type == "measured" & Outlier]) + 
  geom_line(aes(group=1, color = "Inferred"),
            data= function(x)  x[type == "inferred"],
            size=1, alpha=1, show_guide = TRUE) + 
  geom_line(aes(group=1, color = "Predicted"),
            data=function(x)  x[type == "predicted"],
            size=1, alpha=.3) +
  scale_color_manual(values=c("black", col_concise_shallow, "blue"),
                     breaks=c("Data", "Inferred", "Predicted"),
                     name=NULL,
                     guide = guide_legend(order=1,
                                          override.aes = list(
                       linetype = c(rep("blank", 1), "solid", "solid"),
                       shape = c(rep(16, 1), NA, NA)))) + 
  scale_alpha_continuous(guide=guide_legend(override.aes = list(linetype="blank",
                                                                shape=16,
                                                                color="black")))+
  ylab("Branchpoint log-odds") + 
  xlab(NULL) + 
  facet_wrap(~feature_name, ncol=2, scales="free_x") + 
  theme_cowplot() + 
  theme(## legend.justification = c(1, 0),
        ## legend.position = c(1, 0),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal") + 
 ## theme(legend.position = "bottom") +
  theme(strip.background = element_rect(colour="white", fill="white")) # remove tabs from facets
plm

## main figure - position

## supplementary figure - position
plt_pos_effect_supp <- (plm + facet_wrap(~feature_name, ncol=3, scales="free_x")) %+%
  dtpos_all[primary==FALSE]

save_plot(file.path(plt_dir, "pos_effects.pdf"), plm,
          ncol=2, nrow=2, base_height=2, base_aspect_ratio=1.5,
          dpi=600)

save_plot_mul(file.path(plt_dir, "fig3_sup"), c("pdf", "png", "eps"), plt_pos_effect_supp,
          ncol=3, nrow=2, base_height=2, base_aspect_ratio=1.5,
          dpi=600)

## save_plot(file.path(plt_dir, "pos_effects.png"), plt_pos_effect,
##           ncol=4, nrow=1, base_height=3.5, base_aspect_ratio=0.7,
##           dpi=600)
