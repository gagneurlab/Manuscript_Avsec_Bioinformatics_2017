## Idea - for each plot, you need two essential components: data
## (returned by some data function and the plot itself)

source("Scripts/Figures/func/pr_roc.R")

get_glmnet <- function(cache=TRUE) {
  message("Get glmnet")
  cache_file <- "data/eclip/processed/predictions/boostrap_glmnet_auc_auprc.csv"
  if (cache==TRUE & file.exists(cache_file)) {
    return(fread(cache_file))
  }
  dt_pred <- get_predictions()
  dt_pred <- dt_pred[grepl("kmer", method)]
  dt_pred[, method := str_replace(method, "_glmnet", "-glmnet")]
  dt_pred[, method := str_replace(method, "-no_position", "")]
  dt_pred[, method := str_replace(method, "position", "pos")]
  dt_pred[, method := str_replace(method, "-w", "")]
  set.seed(42)
  message("bootstrap auc & auprc")
  dtb <- dt_pred[, bootstrap_metric(y_true, y_pred, 
                                    metric_fn = list("auc"= cem_auc,
                                                     "auprc"= cem_auprc),
                                    n_bootstrap=200, return_summary=FALSE),
                 by = .(rbp, method)]
  dtb <- melt(dtb, measure.vars=c("auc", "auprc"),
              variable.name="metric")
  ## dtb[, task := "eClip_subset"]
  dtb %>% setnames("method", "Method")
  if (cache==TRUE) {
    write_csv(dtb, cache_file)
  }
  return(dtb)
}

get_dtb_eclip_subset <- function(cache=TRUE) {
  ## deepnn
  dtb_clip <- get_eClip_metric(cache=cache, which="subset")
  dtb_glmnet <- get_glmnet(cache=cache)

  dtb_clip %>% setnames("subtask", "rbp")
  dtb_clip[ ,task := NULL]

  dtb <- rbindlist(list(dtb_clip, dtb_glmnet), use.names=TRUE)
  message("Get deepnn model")
  rbps <- c("UPF1", "PUM2", "DDX3X", "NKRF", "TARDBP", "SUGP2")
  dtb[, rbp := fct_relevel(rbp, rbps)]

  dtb <- dtb[!grepl("relu", as.character(Method))]
  ## re-naming
  mname_hash <- c("kmer-glmnet" = "kmer-glmnet",
                  "kmer-glmnet_pos" = "kmer-glmnet_pos",
                  "no_pos" = "DeepNN",
                  "gam" = "DeepNN_pos_spline-t"
                  )
  hash_name <- function(x,h) {
    fct_relevel(h[as.character(x)], h)
  }
  dtb[, Method := hash_name(Method, mname_hash)]

  message("Done!")
  return(dtb)
}
