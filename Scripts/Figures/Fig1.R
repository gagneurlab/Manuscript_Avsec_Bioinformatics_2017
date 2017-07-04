#'---
#' title: Create figure 2
#' author: Ziga Avsec
#' wb:
#'  input: ["data/Concise/Splice_branchpoints/processed/pr_roc/roc_curves.csv",
#'          "data/Concise/Splice_branchpoints/processed/pr_roc/pr_curves.csv"]
#'  output: ["data/plots/Concise/Paper/fig2_roc_pr.pdf",
#'           "data/plots/Concise/Paper/fig2_roc_pr.png"]
#'---


plt_dir <- "./data/plots/Concise/Paper/"
dir.create(plt_dir, showWarnings = FALSE)

save_dir <- "data/Concise/Splice_branchpoints/processed/pr_roc/"
roc_curves <- fread(file.path(save_dir, "roc_curves.csv"))
pr_curves <- fread(file.path(save_dir, "pr_curves.csv"))
source("Scripts/Concise/Paper/config.R")
library(gridExtra)
#+
Figure1C=ggplot(roc_curves, aes(x = FPR,y = TPR, col = method_perf)) + 
  geom_line() + 
  ## scale_color_manual(values=nt_cols) +
  legend_top_ver + 
  ggtitle(label="ROC") +
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(colour = "black", size=.8, linetype=1, inherit.blank = TRUE),
    legend.justification = c(1, 0), legend.position = c(1, 0), 
    )

#' 
#' ### Precision figures
pr_curves <- pr_curves[!(Recall == 0 & Precision == 0)] # Weird first point

Figure1D=ggplot(pr_curves, aes(x = Recall,y = Precision, col = method_perf)) + 
  geom_line() + 
  ## scale_color_manual(values=nt_cols) +  
  legend_top_ver + 
  ggtitle(label="Precision Recall") +
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border =  element_rect(colour = "black", size=.8, linetype=1, inherit.blank = TRUE),
    legend.justification = c(0, 0), legend.position = c(0, 0), 
    )

plt <- plot_grid(Figure1C, Figure1D, labels = c("a", "b"))
plt


save_plot(file.path(plt_dir, "fig2_roc_pr.pdf"), plt,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 1,
          base_height=3.5,
          # each individual subplot should have an aspect ratio of 1.3
          )

save_plot(file.path(plt_dir, "fig2_roc_pr.png"), plt,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 1,
          base_height=3.5,
          )

y <- 1
