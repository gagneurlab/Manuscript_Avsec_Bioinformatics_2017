#'---
#' title: helper
#'---
RBP_LIST <- c("UPF1", "PUM2", "DDX3X", "NKRF", "TARDBP", "SUGP2")

library(PRROC)
read_csv_pred <- function(file) {
  dt <- fread(file)
  if ("V1" %in% colnames(dt)) {
    dt[, V1 := NULL]
  }
  dt[, method := gsub("\\.csv$", "", basename(file))]  
  dt[, y_true := as.numeric(y_true)]
  return(dt)  
}

pred_csv2dt_roc <- function(file) {
  dt <- read_csv_pred(file)
  roc <- dt[, PRROC::roc.curve(y_pred[y_true==1], y_pred[y_true==0], curve = TRUE)]
  roc_curves <- data.table(roc$curve)
  colnames(roc_curves) <- c("FPR","TPR","cut")
  roc_curves[, method := unique(dt$method)]
  roc_curves[, auROC := roc$auc]
  roc_curves[, method_perf := paste0(method, ", ", sprintf("%.3f", auROC))]
  return(roc_curves)
}
pred_csv2dt_pr<- function(file) {
  dt <- read_csv_pred(file)
  if (uniqueN(dt$y_pred) == 1) {
    warning("only one predicted value for file: ", file)
    return(NULL)
  }

  pr <- dt[, PRROC::pr.curve(y_pred[y_true==1],
                             y_pred[y_true==0], curve = TRUE)]
  pr_curves <- data.table(pr$curve)
  colnames(pr_curves) <- c("Recall","Precision", "cut")    

  pr_curves[, method := unique(dt$method)]
  pr_curves[, auPR := pr$auc.integral]
  pr_curves[, method_perf := paste0(method, ", ", sprintf("%.3f", auPR))]
  return(pr_curves)
}


get_roc_curves <- function() {
  BASE_DIR <- "data/eclip/processed/predictions"
  rbp_dirlist <- list.dirs(BASE_DIR, recursive = FALSE, full.names=T)
  roc_curves <- lapply(rbp_dirlist, function(dir) {
    rbp <- basename(dir)
    dt <- lapply(list.files(dir, full.names=T), pred_csv2dt_roc) %>%
      rbindlist(use.names = TRUE)
    dt[, rbp := rbp]
    return(dt)
  }) %>% rbindlist(use.names = TRUE)
  roc_curves[, rbp := fct_relevel(rbp, RBP_LIST)]
  return(roc_curves)
}


get_pr_curves <- function() {
  BASE_DIR <- "data/eclip/processed/predictions"
  rbp_dirlist <- list.dirs(BASE_DIR, recursive = FALSE, full.names=T)
  pr_curves <- lapply(rbp_dirlist, function(dir) {
    rbp <- basename(dir)
    dt <- lapply(list.files(dir, full.names=T), pred_csv2dt_pr) %>%
      rbindlist(use.names = TRUE)
    dt[, rbp := rbp]
    return(dt)
  }) %>% rbindlist(use.names = TRUE)  
  pr_curves[, rbp := fct_relevel(rbp, RBP_LIST)]
  return(pr_curves)
}

get_predictions <- function() {
  BASE_DIR <- "data/eclip/processed/predictions"
  rbp_dirlist <- list.dirs(BASE_DIR, recursive = FALSE, full.names=T)
  pr_curves <- lapply(rbp_dirlist, function(dir) {
    rbp <- basename(dir)
    dt <- lapply(list.files(dir, full.names=T), read_csv_pred) %>%
      rbindlist(use.names = TRUE)
    dt[, rbp := rbp]
    return(dt)
  }) %>% rbindlist(use.names = TRUE)  
  pr_curves[, rbp := fct_relevel(rbp, RBP_LIST)]
  return(pr_curves)
}


plot_roc <- function(roc_curves) {
  Figure1C=ggplot(roc_curves, aes(x = FPR,y = TPR, col = method)) + 
    geom_line() + 
    ## scale_color_manual(values=nt_cols) +
    ## legend_top_ver + 
    ggtitle(label="ROC") +
    facet_grid(.~rbp) + 
    theme(
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      panel.border =  element_rect(colour = "black", size=.8, linetype=1, inherit.blank = TRUE),
      ## legend.justification = c(1, 0), legend.position = c(1, 0), 
      )
  return(Figure1C)
}


plot_pr <- function(pr_curves) {
  Figure1D=ggplot(pr_curves, aes(x = Recall,y = Precision, col = method)) + 
    geom_line() + 
    ## scale_color_manual(values=nt_cols) +  
    ## legend_top_ver + 
    facet_grid(.~rbp) + 
    ggtitle(label="Precision Recall") +
    theme(
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      panel.border =  element_rect(colour = "black", size=.8, linetype=1, inherit.blank = TRUE),
      ## legend.justification = c(0, 0), legend.position = c(0, 0), 
      )
  return(Figure1D)
}
