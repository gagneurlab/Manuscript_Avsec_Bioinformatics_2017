#'---
#' title: Analyze branchpointer performance
#' author: Ziga Avsec
#' wb:
#'  input: ["data/Concise/Splice_branchpoints/test_predictions/branchpointer.csv",
#'          "data/Concise/Splice_branchpoints/test_predictions/glmnet.csv",
#'          "data/Concise/Splice_branchpoints/test_predictions/concise_deep.csv",
#'          "data/Concise/Splice_branchpoints/test_predictions/concise_shallow.csv"]
#'  output: ["data/plots/Concise/Splice_branchpoints/roc_pr.pdf",
#'           "data/plots/Concise/Splice_branchpoints/roc_pr.png",
#'           "data/Concise/Splice_branchpoints/processed/pr_roc/roc_curves.csv",
#'           "data/Concise/Splice_branchpoints/processed/pr_roc/pr_curves.csv",
#'           "data/Concise/Splice_branchpoints/processed/pr_roc/bootstrap_auc.csv",
#'           "data/Concise/Splice_branchpoints/processed/pr_roc/bootstrap_auprc.csv"]
#'---
opts_chunk$set(echo=FALSE, cache=F, results = 'hide', messages = FALSE)
options(width=140)
#+ 
plt_dir <- "./data/plots/Concise/Splice_branchpoints/"
dir.create(plt_dir, showWarnings = FALSE)
source("Scripts/Concise/Splice_branchpoints/config.R")
library(cowplot)
#'
#' ## Conclusions
#'
#' - we can reproduce the performance figures of branchpointer
#' - Metrics from the paper:
#'     - AUC: 0.941
#'     - auPR: 0.617
#' 
#'----------------------
read_csv_pred <- function(file) {
  dt <- fread(file)
  if ("V1" %in% names(dt)) {
    dt[, V1 := NULL]
  }
  dt[, method := gsub("\\.csv$", "", basename(file))]  
  return(dt)  
}

pred_csv2dt_roc <- function(file) {
  dt <- read_csv_pred(file)
  roc <- dt[, PRROC::roc.curve(y_pred[y_true=="HC"], y_pred[y_true=="NEG"], curve = TRUE)]
  roc_curves <- data.table(roc$curve)
  colnames(roc_curves) <- c("FPR","TPR","cut")
  roc_curves[, method := unique(dt$method)]
  roc_curves[, auROC := roc$auc]
  roc_curves[, method_perf := paste0(method, ", ", sprintf("%.3f", auROC))]
  return(roc_curves)
}
pred_csv2dt_pr<- function(file) {
  dt <- read_csv_pred(file)
  pr <- dt[, PRROC::pr.curve(y_pred[y_true=="HC"],
                             y_pred[y_true=="NEG"], curve = TRUE)]
  pr_curves <- data.table(pr$curve)
  colnames(pr_curves) <- c("Recall","Precision", "cut")
  pr_curves[, method := unique(dt$method)]
  pr_curves[, auPR := pr$auc.integral]
  pr_curves[, method_perf := paste0(method, ", ", sprintf("%.3f", auPR))]
  return(pr_curves)
}

files <- list.files("data/Concise/Splice_branchpoints/test_predictions/", full.names=T)
files <- files[!grepl("gam_0|_mul_", files)]

## ---------------------------------
dtp <- lapply(files, read_csv_pred) %>% rbindlist(fill=TRUE)
dtp[, id := NULL]


dtp[ , method] %>% unique
dtp[, method := fct_recode(method,
                           concise_deep_relu = "deep_relu_0" ,
                           concise_shallow_relu = "shallow_relu_0"
                           )]
dtp[, method := fct_relevel(method, c("glmnet",
                                      "concise_shallow_relu", "concise_shallow",
                                      "branchpointer",
                                      "concise_deep_relu", "concise_deep"
                                      ))]
## TODO - bootstrap
set.seed(42)
dtp_b_auc <- dtp[, bootstrap_measure(y_true=="HC", y_pred,
                                     cem_auc, 200, return_summary=FALSE),
                 by = method]
dtp_b_auc %>% setnames("V1", "auc")
dtp_b_auprc <- dtp[, bootstrap_measure(y_true=="HC", y_pred,
                                       cem_auprc, 200, return_summary=FALSE),
                 by = method]
dtp_b_auprc %>% setnames("V1", "auprc")

library(ggsignif)
q_auc <- ggplot(dtp_b_auc, aes(x = method, y = auc)) +
  geom_boxplot() + geom_jitter(alpha=0.1) +
  geom_signif(comparisons = list(c("concise_shallow_relu", "concise_shallow"),
                                 c("concise_deep_relu", "concise_deep"),
                                 c("branchpointer", "concise_deep")
                                 ),
              step_increase = 0.05)

q_auprc <- ggplot(dtp_b_auprc, aes(x = method, y = auprc)) +
  geom_boxplot() + geom_jitter(alpha=0.1) +
  geom_signif(comparisons = list(c("concise_shallow_relu", "concise_shallow"),
                                 c("concise_deep_relu", "concise_deep"),
                                 c("branchpointer", "concise_deep")
                                 ),
              step_increase = 0.05)
#+ fig.width=6
plot_grid(q_auc, q_auprc, ncol=2)
## ---------------------------------

roc_curves <- lapply(files, pred_csv2dt_roc) %>% rbindlist
pr_curves <- lapply(files, pred_csv2dt_pr) %>% rbindlist

save_dir <- "data/Concise/Splice_branchpoints/processed/pr_roc/"
write_csv(roc_curves, file.path(save_dir, "roc_curves.csv"))
write_csv(pr_curves, file.path(save_dir, "pr_curves.csv"))

write_csv(dtp_b_auc, file.path(save_dir, "bootstrap_auc.csv"))
write_csv(dtp_b_auprc, file.path(save_dir, "bootstrap_auprc.csv"))


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


save_plot(file.path(plt_dir, "roc_pr.pdf"), plt,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 1,
          base_height=3.5,
          # each individual subplot should have an aspect ratio of 1.3
          )

save_plot(file.path(plt_dir, "roc_pr.png"), plt,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 1,
          base_height=3.5,
          )

y <- 1
