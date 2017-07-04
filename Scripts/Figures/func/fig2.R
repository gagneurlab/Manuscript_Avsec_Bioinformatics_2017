## Idea - for each plot, you need two essential components: data
## (returned by some data function and the plot itself)

get_dtb <- function() {
  fix_names <- function(x) {
    x %>% str_replace("_2", "") %>% str_replace("scalar_", "") %>% str_replace("position", "pos")
  }
  ##+
  ## deepnn
  message("Get deepnn model")
  dtb_auc <- fread("data/encode/eclip/processed/predictions/bootstrap_auc.csv")
  dtb_auprc <- fread("data/encode/eclip/processed/predictions/bootstrap_auprc.csv")
  dtb_auc[, V1 := NULL]
  dtb_auprc[, V1 := NULL]
  dtb_auc %>% setnames("method", "Method")
  dtb_auprc %>% setnames("method", "Method")
  dtb_auc[, Method := fix_names(Method)]
  dtb_auprc[, Method := fix_names(Method)]
  dtb_auc[, Method] %>% unique

  rbps <- c("UPF1", "PUM2", "DDX3X", "NKRF", "TARDBP", "SUGP2")
  dtb_auc[, rbp := fct_relevel(rbp, rbps)]
  dtb_auprc[, rbp := fct_relevel(rbp, rbps)]

  ## glmnet
  message("Get glmnet")
  dt_pred <- get_predictions()
  dt_pred <- dt_pred[grepl("kmer", method)]
  dt_pred[, method := str_replace(method, "_glmnet", "-glmnet")]
  dt_pred[, method := str_replace(method, "-no_position", "")]
  dt_pred[, method := str_replace(method, "position", "pos")]
  dt_pred[, method := str_replace(method, "-w", "")]
  set.seed(42)
  message("bootstrap auc")
  dt_pred_auc <- dt_pred[, .(auc = bootstrap_metric(y_true, y_pred, cem_auc,
                                                    n_bootstrap=200, return_summary=FALSE)),
                         by = .(rbp, method)]
  message("bootstrap auprc")
  dt_pred_auprc <- dt_pred[, .(auprc = bootstrap_metric(y_true, y_pred, cem_auprc,
                                                        n_bootstrap=200, return_summary=FALSE)),
                           by = .(rbp, method)]
  message("reformat")
  dt_pred_auc %>% setnames("method", "Method")
  dt_pred_auprc %>% setnames("method", "Method")

  dtb_auc <- rbindlist(list(dt_pred_auc, dtb_auc), use.names=TRUE, fill=TRUE)
  dtb_auprc <- rbindlist(list(dt_pred_auprc, dtb_auprc), use.names=TRUE, fill=TRUE)

  ## TODO - from here 
  dtb_auc[, Method] %>% unique


  dtb_auc[, Method := fct_relevel(Method, c("kmer-glmnet", "kmer-glmnet_pos"))]
  dtb_auprc[, Method := fct_relevel(Method, c("kmer-glmnet", "kmer-glmnet_pos"))]

  ## TODO - rename
  dtb_auc[, Method := fct_recode(Method, "DeepNN_pos_spline-t"="DeepNN_pos_gam")]
  dtb_auprc[, Method := fct_recode(Method, "DeepNN_pos_spline-t"="DeepNN_pos_gam")]
  message("Done!")
  return(list(dtb_auc=dtb_auc, dtb_auprc=dtb_auprc))
}
