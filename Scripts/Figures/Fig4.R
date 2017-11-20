#'---
#' title: Create figure 4
#' author: Žiga Avsec
#' wb:
#'  input: ["data/Splice_branchpoints/trials/train_history/gam_vs_relu.csv"]
#'  output: ["data/plots/fig4.pdf",
#'           "data/plots/fig4.png"]
#'---
#' 
knitr::opts_chunk$set(echo=FALSE, cache=F, results = 'hide', messages = FALSE)
options(width=140)
#+ 
plt_dir <- "data/plots"
# plt_dir <- "~/projects-work/spline_trans/plots/"
dir.create(plt_dir, showWarnings = FALSE)
library(cowplot)
source("Scripts/Figures/config.R")
## 
source_all(dir = "Scripts/Splice_branchpoints/plot_functions")

col_scheme_deep <- scale_colour_manual(values = c(col_relu_deep, col_concise_deep),
                                       labels = c("PL", "Spline"),
                                       name="Transformation")
col_scheme_shallow <- scale_colour_manual(values = c(col_relu_shallow,
                                                     col_concise_shallow),
                                       labels = c("PL", "Spline"),
                                       name="Transformation")
col_scheme_all <- scale_colour_manual(values = c(col_relu_shallow, col_concise_shallow,
                                                 col_relu_deep, col_concise_deep))
rbps <- c("UPF1", "PUM2", "DDX3X", "NKRF", "TARDBP", "SUGP2")
## --------------------------------------------

## 1. two plots with relu vs gam comparision (same scatterplot as in fig2)

dtm <- get_metric_comparison(1e-4)
df1 <- dtm[metric == "auprc"]


in_top <- function(x, n = 20) {
  return(rank(- x) <= n)
}
df1[, add_label := (subtask %in% rbps) |  !signif.gam__relu|
        (task == "eCLIP" & in_top(gam- relu, 15)) |
        (task == "CLIP" & in_top(gam- relu, 10))
  , by = task] # TOP 20 from 
df1[, from_before := subtask %in% rbps & task != "CLIP"]



dummy <- data.table(task = c("eCLIP", "eCLIP",
                             "CLIP", "CLIP",
                             "branchpoint", "branchpoint"),
                    x = c(0.33, 0.95,
                          0.5, 0.97,
                          0.57, 0.68)
                    )

df1[, task_name := ifelse(task == "CLIP", "CLIP (iDeep)", task)]
df1[, task_name := fct_relevel(task_name, c("eCLIP", "CLIP (iDeep)"))]
dummy[, task_name := ifelse(task == "CLIP", "CLIP (iDeep)", task)]
dummy[, task_name := fct_relevel(task_name, c("eCLIP", "CLIP (iDeep)"))]



plt_test_gam_vs_relu <- ggplot(df1,
                               aes(x = relu,
                                   y = gam,
                                   color = !signif.gam__relu))+#factor(sign(gam - relu) * signif.gam__relu, levels = c(1, 0, -1)))) +
  geom_abline(alpha=0.2) +
  geom_point(size=2, alpha=0.7) +
  ## geom_text_repel(data = df1[add_label==TRUE ], aes(label=subtask, color = from_before),
  ##                 box.padding = unit(0.4, "lines"),
  ##                 point.padding = unit(0.4, "lines"),
  ##                 alpha=0.5) +
  facet_wrap(~task_name, scales="free") + 
  scale_color_manual(values=c("black", "grey", "blue", mypal[3], "black", mypal[1] , "black", "red"), # mypal[7] - alternatively for grey
                     ## limits = c(3, 4, 5),
                     guide = guide_legend("Wilcoxon test"),
                     drop=FALSE,
                     labels = c(bquote(list(p[adj] < 10^-4)),
                                "Not significant",
                                bquote(list(p[adj] < 0.0001, bar(x) > bar(y))),
                                "", "")) + 
  ## xlim(c(0.29, .97)) + 
  ## ylim(c(0.29, .97)) + 
  ## xlim(c(0.68, 1)) +
  ## ylim(c(0.68, 1)) +
  geom_blank(aes(x=x, y=x, color=TRUE), data=dummy) + 
  ylab("auPR with spline\ntransformation of dist.") +
  xlab("auPR with piecewise lin. transformation of dist.") + 
  legend_top_ver +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0))
  ## white_facets
plt_test_gam_vs_relu



## --------------------------------------------
## loss-curves
dtm_w_mean <- gam_vs_relu_loss_curves_dt()
type_hash <- c("relu"="NN w/ PLT", "gam"="NN w/ ST")

## TODO - set col-order first
## dtm_w_mean[, method := as.factor_keep_order(paste(type_hash[as.character(type)], depth))]
## dtm_w_mean[, method := fct_rev(method)]

loss_deep <- ggplot(dtm_w_mean[metric %in% c("loss") &
                                 epoch >= 0 & epoch < 50 &
                                   dataset == "validation"][depth == "deep"],
                    aes(x = epoch, y = value, 
                 color = type, 
                 alpha = is_mean,
                 group = interaction(tid, dataset, type),
                 )) + 
  geom_line() + 
  theme_bw()
loss_shallow <- loss_deep %+%
  dtm_w_mean[metric %in% "loss" & epoch >= 0 & dataset == "validation"][depth == "shallow"]


## (loss_deep + aes(color = method, alpha=NULL, group=NULL) + col_scheme_all
##    ) %+% dtm_w_mean

yl <- ylab("Binary crossentropy \nloss (validation)")
## yl <- ylab("Binary crossentropy loss")
xl <- xlab("Epoch")

title_style <- theme(plot.title = element_text(face="plain",  size=12))

lgd <- get_legend(loss_shallow + scale_color_manual(name="Method",
                                                   values = c(col_relu_shallow,
                                                              col_concise_shallow,
                                                              col_relu_deep,
                                                              col_concise_deep),
                                                   limits = c("shallow_relu",
                                                              "shallow_gam",
                                                              "deep_relu",
                                                              "deep_gam")) +
                    theme_cowplot() + 
                    legend_top + 
                      guides(alpha="none"))

yl_loss <- ggdraw() + draw_label("Binary crossentropy \nloss (validation)", angle=90)
xl_label  <- ggdraw() + draw_label("Epoch")

loss_plt <- plot_grid(plot_grid(yl_loss,
                               loss_shallow + ylim(c(NA, 0.18)) + theme_cowplot() +
                                 legend_off + ylab(NULL) + xl + ggtitle("shallow") +
                                   xlab(NULL) + 
                                   title_style +
                                   col_scheme_shallow,
                               loss_deep + ylim(c(NA, 0.13)) + theme_cowplot() +
                                 theme(legend.justification = c(1, 1),
                                       legend.position = c(1, 1)) + 
                                 scale_alpha_discrete(guide = FALSE) +
                                 xl + ggtitle("deep") +
                                   ylab(NULL) + 
                                   title_style +
                                   xlab(NULL)+
                                   col_scheme_deep,
                               nrow=1,
                               rel_widths = c(0.3, 1, 1)
                               ),
                     xl_label,
                     ncol=1,
                     rel_heights=c(1,0.1)
                     )
a <- 1

## --------------------------------------------
## boxplot
dt_eval <- get_trials_df_eval()
dt_eval[, type := fct_relevel(type, c("relu", "gam"))]
dt_eval[, type := fct_recode(type, "PL"="relu", "Spline"="gam")]
pl <- ggplot(dt_eval, aes(x = type, y = eval.auprc)) +
  geom_boxplot(aes(color = type)) +
  geom_jitter(aes(color = type), alpha = 0.1) +
  ## xlab("Position transformation") +
  xlab(NULL)

## TODO 
pl2 <- ggplot(dt_eval, aes(x = type, y = eval.auprc)) +
  geom_boxplot(aes(color = type)) +
  geom_jitter(aes(color = type), alpha = 0.1) +
  ## xlab("Position transformation") +
  xlab(NULL)

pl_auprc_deep <- pl %+% dt_eval[depth == "deep"] + ggtitle("deep")  + title_style 
pl_auprc_shallow <- pl %+% dt_eval[depth == "shallow"] + ggtitle("shallow") + title_style


yl_box <- ggdraw() + draw_label("auPR (validation)", angle=90)
xl_box_label  <- ggdraw() + draw_label("Distance transformation type")

perf_boxplot <- plot_grid(yl_box,
                          pl_auprc_shallow + legend_off+
                            ylab(NULL) +
                            col_scheme_shallow,
                          pl_auprc_deep + legend_off +
                            ylab(NULL) + 
                            ## ylab("Area under\nprecision-recall curve") +
                            col_scheme_deep,
                          nrow=1,
                          rel_widths = c(0.15, 1, 1))
perf_boxplot <- plot_grid(perf_boxplot,
                          xl_box_label,
                          ncol=1,
                          rel_heights = c(1,0.1))

#+ fig.width = 12, fig.height = 8

plt_global <- plot_grid(plt_test_gam_vs_relu,
                        plot_grid(perf_boxplot, loss_plt,
                                  ncol=2,
                                  labels = c("d", "e"),
                                  rel_widths=c(1,1.2),
                                  align="h"
                                  ),
                        labels = "a",
                        rel_widths = c(2.5, 1),
                        ncol=1
                        ## lgd
                        )
plt_global

save_plot_mul(file.path(plt_dir, "fig4"), c("png", "pdf", "eps"), plt_global,
          ncol=3, nrow=2, base_height=2.8, base_aspect_ratio=1,
          dpi=600)

## create S3...
## save_plot(file.path(plt_dir, "S_fig4.pdf"), plt1_o,
##           ncol=1, nrow=1, base_height=8, base_aspect_ratio=1.7/2,
##           dpi=600)
## save_plot(file.path(plt_dir, "S_fig4.png"), plt1_o,
##           ncol=1, nrow=1, base_height=8, base_aspect_ratio=1.7/2,
##           dpi=600)
