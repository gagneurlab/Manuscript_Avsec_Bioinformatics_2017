### Figure 3
## 
## - [x] update the data
## - [x] create a function to read in each individual dataset
##   - one for each dataset
##   - choose which combindation of inputs to use (useful for the task afterwards)
##   - 
## - [x] create a function to merge them all
##   - have a nice, clean representation
## - [x] boostrap all-together in parallel (user bplapply)
##   - look for parallel split-apply-combine paradigm in data.table
## - [x] save the results - (one big table)
##   - Put the previous code into one script (as it will probably be a long-running job)
## 	- Use snake-make?
## - [x] create a function to load the cached results

## 1. Predictive performance for iDeep
library(dtplyr)
DIR_ROOT <- "data/eclip/"

## data format 3.
get_iCLIP_metrics <- function() {
  dti <- fread("data/eclip/processed/predictions/iClip-iDeep_boostrap_Amin_3.csv")
  dti[, V1 := NULL]
  dti %>% setnames("method", "Method")
  dti[, Method := str_replace(Method, "_scaler", "")]
  dti[, Method := str_replace(Method, "_position", "_pos")]
  dti[,i := 1:.N , by = .(Method, rbp)]
  dti[, metric := "auc"]
  setnames(dti, "auc", "value")

  dti_auprc<- fread("data/eclip/processed/predictions/iClip-iDeep_boostrap_Amin_auprc_3.csv")
  dti_auprc[, V1 := NULL]
  dti_auprc %>% setnames("method", "Method")
  dti_auprc[, Method := str_replace(Method, "_scaler", "")]
  dti_auprc[, Method := str_replace(Method, "_position", "_pos")]
  dti_auprc[,i := 1:.N , by = .(Method, rbp)]
  dti_auprc[, metric := "auprc"]
  setnames(dti_auprc, "auprc", "value")
  dtim <- rbindlist(list(dti, dti_auprc))

  dtim <- dtim %>% rename(subtask=rbp) %>% mutate(i=NULL, task="iCLIP")
  dtim[, Method := gsub("iDeep_pos_", "", Method)]
  dtim[, Method := gsub("iDeep", "no_pos", Method)]
  dtim <- dtim %>% setcolorder_first(c("task", "subtask", "Method", "metric", "value"))
  return(dtim)
}

## TODO - do the same for eCLIP with 2 features

## 2. Predictive performance for eClip
get_eClip_metric <- function(n_cores=4, n_bootstrap=200, cache=TRUE, which="ext") {
  ## ugly but ok...
  if (which == "ext") {
    cache_file <- "data/eclip/processed/predictions/boostrap_auc_auprc_ext.csv"
    exp_list <- c("DeepNN_ext",
                  "DeepNN_scalar_position_ext_gam",
                  "DeepNN_scalar_position_ext_relu")
    name_hash <- function(exp) {
      return(ifelse(exp == "DeepNN_ext", "no_pos", gsub("DeepNN_scalar_position_ext_", "", exp)))
    }
    base_from <- "DeepNN_scalar_position_relu"
    rbp_list <- list.dirs("data/eclip/processed/predictions",
                          recursive=FALSE, full.names = FALSE)
  } else if (which == "subset") {
    cache_file <- "data/eclip/processed/predictions/boostrap_auc_auprc_2.csv"
    exp_list <- c("DeepNN_2",
                  "DeepNN_scalar_position_gam_2",
                  "DeepNN_scalar_position_relu_2")    
    name_hash <- function(exp) {
      exp <- gsub("_2", "", exp)
      return(ifelse(exp == "DeepNN", "no_pos", gsub("DeepNN_scalar_position_", "", exp)))
    }
    rbp_list <- c("UPF1", "PUM2", "DDX3X", "NKRF", "TARDBP", "SUGP2")
  }


  get_eClip_metric_single <- function(rbp, exp, n_bootstrap = 200) {
    path <- file.path(DIR_ROOT, "/processed/predictions/", rbp, paste0(exp, ".csv"))
    dt <- fread(path)
    dt[, V1 := NULL]
    dt[, subtask := rbp]
    dt[, task := "eClip"]
    dt[, Method := name_hash(exp)]

    ## boostrap and get the metric
    dt <- dt[, bootstrap_metric(y_true, y_pred,
                                metric_fn = list("auc"= cem_auc,
                                                 "auprc"= cem_auprc),
                                return_summary=FALSE,
                                n_bootstrap=n_bootstrap),
             by = .(task, subtask, Method)]

    dt <- melt(dt, measure.vars=c("auc", "auprc"),
               variable.name="metric")
    return(dt)
  }

  if (cache==TRUE & file.exists(cache_file)) {
    dt <- fread(cache_file)
    dt[, task := gsub("eClip", "eCLIP", task)]
    return(dt)
  }


  dt <- pbmclapply(rbp_list, function(rbp) {
    lapply(exp_list, function(exp)
      get_eClip_metric_single(rbp, exp, n_bootstrap)) %>%
      rbindlist
  }, mc.cores=n_cores) %>% rbindlist
  if (cache==TRUE) {
    write_csv(dt, cache_file)
  }
  return(dt)
}
## write_csv(dte, "data/eclip/processed/predictions/boostrap_auc_auprc_ext.csv")

## 3. Predictive performance for branchpoint prediction task 
get_BP_metric <- function() {
  load_dir <- "data/Concise/Splice_branchpoints/processed/pr_roc/"
  dtp_b_auc <- fread(file.path(load_dir, "bootstrap_auc.csv"))
  dtp_b_auprc <- fread(file.path(load_dir, "bootstrap_auprc.csv"))
  dtp_b_auc <- dtp_b_auc %>% setnames("auc", "value") %>% mutate(metric="auc")
  dtp_b_auprc <- dtp_b_auprc %>% setnames("auprc", "value") %>% mutate(metric="auprc")

  dtp <- rbindlist(list(dtp_b_auc, dtp_b_auprc))
  dtp[, task := "branchpoint"]
  ## TODO - rename 

  dtp <- dtp[grepl("concise", method)]
  dtp[, method := gsub("concise_", "", method)]

  dtp <- dtp %>% separate(method, c("subtask", "Method"), sep="_", extra="merge", fill="right")
  dtp[, Method := ifelse(is.na(Method), "gam", Method)]

  
  dtp <- dtp %>% setcolorder_first(c("task", "subtask", "Method", "metric", "value"))
  ## TODO - restructure the table
  return(dtp)
}

## data format 1. (not applicable here)
## y_pred, y_true, method (no_pos, gam, relu), task (eCLip, iCLIP, BP), subtask (RBP, deepBP)

## data format:
## method, task, subtask, metric (auc, auprc), value (bootstrapped samples)
get_all_metric <- function() {
  dt_iclip <- get_iCLIP_metrics()

  dt_eclip <- get_eClip_metric(n_cores=7)

  dt_bp <- get_BP_metric()
  dt <- rbindlist(list(dt_iclip, dt_eclip, dt_bp))
  dt[, i := 1:.N, by = .(task, subtask, metric, Method)]
  return(dt)
}

## always compare no_pos <-> gam, relu <-> gam -> compare p-value and average score


## data format 4.
## method, task, subtask, metric, gam, relu, no_pos, p.gam__relu, p.gam__no_pos, signif.gam__relu (bool), signif.gam__no_pos(bool)
get_metric_comparison <- function(signif_threshold=0.05) {
  dt <- get_all_metric()
  ## data format 3.
  ## method, task, subtask, metric, gam, relu, no_pos
  dtc <- dt %>% spread(Method, value)
  dtc

  dtct <- dtc[, .(gam = mean(gam), relu=mean(relu), no_pos=mean(no_pos),
                  p.gam__relu = wilcox.test(gam, relu)$p.value,
                  p.gam__no_pos = if (task == "branchpoint") as.numeric(NA) else wilcox.test(gam, no_pos)$p.value
                  ), by = .(task, subtask, metric)]

  dtct[, signif.gam__relu := p.adjust(p.gam__relu, method="bonferroni") < signif_threshold]
  dtct[, signif.gam__no_pos:= p.adjust(p.gam__no_pos, method="bonferroni") < signif_threshold]
  return(dtct)
}
