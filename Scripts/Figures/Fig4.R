#'---
#' title: Create figure 3
#' author: Žiga Avsec
#' wb:
#'  input: ["data/Concise/Splice_branchpoints/trials/train_history/gam_vs_relu.csv"]
#'  output: ["data/plots/Concise/Paper/fig4.pdf",
#'           "data/plots/Concise/Paper/fig4.png"]
#'---
#' 
opts_chunk$set(echo=FALSE, cache=F, results = 'hide', messages = FALSE)
options(width=140)
#+ 
plt_dir <- "~/projects-work/spline_trans/plots/"
dir.create(plt_dir, showWarnings = FALSE)
library(cowplot)
source("Scripts/Concise/Paper/config.R")
## 
source_all(dir = "Scripts/Concise/Splice_branchpoints/plot_functions")

col_scheme_deep <- scale_colour_manual(values = c(col_concise_deep, col_relu_deep))
col_scheme_shallow <- scale_colour_manual(values = c(col_concise_shallow,
                                                     col_relu_shallow))
rbps <- c("UPF1", "PUM2", "DDX3X", "NKRF", "TARDBP", "SUGP2")
## --------------------------------------------

## 1. two plots with relu vs gam comparision (same scatterplot as in fig2)

dtm <- get_metric_comparison(1e-4)
df1 <- dtm[metric == "auprc"]


in_top <- function(x, n = 20) {
  return(rank(- x) <= n)
}
df1[, add_label := (subtask %in% rbps) |  !signif.gam__relu|
        (task == "eClip" & in_top(gam- relu, 15)) |
        (task == "iClip" & in_top(gam- relu, 10))
  , by = task] # TOP 20 from 
df1[, from_before := subtask %in% rbps & task != "iClip"]

df1[, task_name := ifelse(task == "iClip", "iClip (iDeep)", task)]
df1[, task_name := fct_relevel(task_name, c("eClip", "iClip (iDeep)"))]

plt_test_gam_vs_relu <- ggplot(df1,
                               aes(x = relu,
                                   y = gam,
                                   color = !signif.gam__relu))+#factor(sign(gam - relu) * signif.gam__relu, levels = c(1, 0, -1)))) +
  geom_abline(alpha=0.2) +
  geom_point(size=2, alpha=0.5) +
  ## geom_text_repel(data = df1[add_label==TRUE ], aes(label=subtask, color = from_before),
  ##                 box.padding = unit(0.4, "lines"),
  ##                 point.padding = unit(0.4, "lines"),
  ##                 alpha=0.5) +
  facet_wrap(~task_name, scales="free") + 
  scale_color_manual(values=c("black", "grey", mypal[3], "black", mypal[1] , "black", "red"), # mypal[7] - alternatively for grey
                     ## limits = c(3, 4, 5),
                     guide = guide_legend("Wilcoxon test"),
                     drop=FALSE,
                     labels = c(bquote(list(p[adj] < 10^-4)),
                                "Not significant",
                                bquote(list(p[adj] < 0.0001, bar(x) > bar(y))),
                                "", "")) + 
  xlim(c(0.29, .97)) + 
  ylim(c(0.29, .97)) + 
  ## xlim(c(0.68, 1)) +
  ## ylim(c(0.68, 1)) +
  ylab("auPR with spline\ntransformation of pos.") +
  xlab("auPR with piece-wise linear transformation of pos.") + 
  legend_top_ver +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0))
  ## white_facets
plt_test_gam_vs_relu



## --------------------------------------------
## loss-curves
ret <- gam_vs_relu_loss_curves_plot(facets=NULL)
yl <- ylab("Binary crossentropy \nloss (validation)")
## yl <- ylab("Binary crossentropy loss")
xl <- xlab("Epoch")

title_style <- theme(plot.title = element_text(face="plain",  size=12))

lgd <- get_legend(ret$shallow + scale_color_manual(name="Method",
                                                   values = c(col_concise_shallow,
                                                              col_relu_shallow,
                                                              col_concise_deep,
                                                              col_relu_deep),
                                                   limits = c("shallow_gam",
                                                              "shallow_relu",
                                                              "deep_gam",
                                                              "deep_relu")) +
                    theme_cowplot() + 
                    legend_top + 
                      guides(alpha="none"))

yl_loss <- ggdraw() + draw_label("Binary crossentropy \nloss (validation)", angle=90)

loss_plt <- with(ret, plot_grid(yl_loss,
                                shallow + ylim(c(NA, 0.18)) + theme_cowplot() +
                                  legend_off + ylab(NULL) + xl + ggtitle("shallow") +
                                    title_style + 
                                    col_scheme_shallow,
                                deep + ylim(c(NA, 0.13)) + theme_cowplot() +
                                  legend_off +  xl + ggtitle("deep") +
                                    ylab(NULL) + 
                                    title_style + 
                                    col_scheme_deep,
                                nrow=1,
                                rel_widths = c(0.2, 1, 1)
                                ))
a <- 1

## --------------------------------------------
## boxplot
dt_eval <- get_trials_df_eval()
dt_eval[, type := fct_relevel(type, c("relu", "gam"))]
pl <- ggplot(dt_eval, aes(x = type, y = eval.auprc)) +
  geom_boxplot(aes(color = type)) +
  geom_jitter(aes(color = type), alpha = 0.1) +
  xlab(NULL)

pl_auprc_deep <- pl %+% dt_eval[depth == "deep"] + ggtitle("deep")  + title_style
pl_auprc_shallow <- pl %+% dt_eval[depth == "shallow"] + ggtitle("shallow") + title_style


yl_loss <- ggdraw() + draw_label("auPR (validation)", angle=90)
perf_boxplot <- plot_grid(yl_loss,
                          pl_auprc_shallow + legend_off+
                            ylab(NULL) +
                            col_scheme_shallow,
                          pl_auprc_deep + legend_off +
                            ylab(NULL) + 
                            ## ylab("Area under\nprecision-recall curve") +
                            col_scheme_deep,
                          nrow=1,
                          rel_widths = c(0.1, 1, 1))

#+ fig.width = 12, fig.height = 8

plt_global <- plot_grid(plt_test_gam_vs_relu,
                        plot_grid(perf_boxplot, loss_plt,
                                  ncol=2,
                                  labels = c("b", "c"),
                                  rel_heights=c(1.2,1),
                                  align="h"
                                  ),
                        labels = "a",
                        rel_widths = c(2.5, 1),
                        ncol=1
                        ## lgd
                        )
plt_global

save_plot_mul(file.path(plt_dir, "fig4"), c("png", "pdf", "eps"), plt_global,
          ncol=3, nrow=2, base_height=3.5, base_aspect_ratio=1,
          dpi=600)

## create S3...
## save_plot(file.path(plt_dir, "S_fig4.pdf"), plt1_o,
##           ncol=1, nrow=1, base_height=8, base_aspect_ratio=1.7/2,
##           dpi=600)
## save_plot(file.path(plt_dir, "S_fig4.png"), plt1_o,
##           ncol=1, nrow=1, base_height=8, base_aspect_ratio=1.7/2,
##           dpi=600)
