#'---
#' title: Bootstrap performance
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
source("Scripts/Concise/Splice_branchpoints/config.R")
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

roc_curves <- lapply(files, pred_csv2dt_roc) %>% rbindlist
pr_curves <- lapply(files, pred_csv2dt_pr) %>% rbindlist

save_dir <- "data/Concise/Splice_branchpoints/processed/pr_roc/"
write_csv(roc_curves, file.path(save_dir, "roc_curves.csv"))
write_csv(pr_curves, file.path(save_dir, "pr_curves.csv"))

write_csv(dtp_b_auc, file.path(save_dir, "bootstrap_auc.csv"))
write_csv(dtp_b_auprc, file.path(save_dir, "bootstrap_auprc.csv"))
