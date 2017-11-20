#'---
#' title: scripts
#' author: Žiga Avsec
#'---
# 
#  input: ["data/Splice_branchpoints/trials/train_history/gam_vs_relu.csv"]

## example:
## 
## ret <- gam_vs_relu_loss_curves_plot()
## yl <- ylab("binary crossentropy loss (validation)")
## legend_off = theme(legend.position = "none")
## with(ret, plot_grid(get_legend(deep+ theme(legend.position = "top") + guides(alpha="none")),
##                     plot_grid(deep + ylim(c(NA, 0.13)) + legend_off + yl,
##                               shallow + ylim(c(NA, 0.18)) + legend_off + yl),
##                     ncol = 1,
##                     rel_heights = c(0.05, 1)))

gam_vs_relu_loss_curves_plot <- function(top_N=10, use_metrics = c("loss"),
                                       facets = ~ depth, 
                                       min_epoch = 0,
                                       deep_max_epoch=50) {
  library(cowplot)
  EXP_DIR = "data/Splice_branchpoints/"
  dt <- fread(file.path(EXP_DIR, "trials/train_history/gam_vs_relu.csv"))  
  dt[, V1 := NULL]
  top_dt <- dt[, .(tid, trial)] %>% unique %>% .[, .SD[1:top_N], by = trial]
  metrics <- names(dt) %>% grep("^val_", ., value = TRUE) %>% gsub("^val_", "", .)
  all_metrics <- paste0(c("", "val_"), metrics)
  # make tidy data
  dt_use = merge(dt, top_dt, by = c("trial", "tid"))
  dtm = melt(dt_use, id.vars = c("tid", "trial", "epoch"))
  dtm <- dtm[trial != "shallow_gam_mul"]
  dtm[, dataset := ifelse(grepl("^val_", variable), "validation", "train")]
  dtm[, metric := gsub("^val_", "", variable)]
  dtm <- separate(dtm, trial, c("depth", "type"))

  dtm_mean = dtm[, .(value = mean(value, na.rm=TRUE)),  by = .(metric, depth, type, epoch, dataset)]
  dtm_mean[, is_mean := TRUE]
  dtm[, is_mean := FALSE]
  dtm_w_mean = rbindlist(list(dtm_mean, dtm), use.names = TRUE, fill = TRUE)
  dtm_w_mean[, tid := ifelse(is.na(tid), 0, tid)]

  dtm_w_mean[, type := fct_relevel(type, "relu")]
  a = ggplot(dtm_w_mean[metric %in% use_metrics &
                          epoch >= min_epoch & epoch < deep_max_epoch &
                          dataset == "validation"][depth == "deep"],
             aes(x = epoch, y = value, 
                 color = type, 
                 alpha = is_mean,
                 group = interaction(tid, dataset, type),
                 )) + 
    geom_line() + 
    theme_bw()
  if (!is.null(facets)) {
    a <- a + facet_grid(facets)
  }
  b = a %+% dtm_w_mean[metric %in% use_metrics &  epoch >= min_epoch & dataset == "validation"][depth == "shallow"]
  return(list(deep = a, shallow = b))
}

gam_vs_relu_loss_curves_dt <- function(top_N=10) {
  library(cowplot)
  EXP_DIR = "data/Splice_branchpoints/"
  dt <- fread(file.path(EXP_DIR, "trials/train_history/gam_vs_relu.csv"))  
  dt[, V1 := NULL]
  top_dt <- dt[, .(tid, trial)] %>% unique %>% .[, .SD[1:top_N], by = trial]
  metrics <- names(dt) %>% grep("^val_", ., value = TRUE) %>% gsub("^val_", "", .)
  all_metrics <- paste0(c("", "val_"), metrics)
  # make tidy data
  dt_use = merge(dt, top_dt, by = c("trial", "tid"))
  dtm = melt(dt_use, id.vars = c("tid", "trial", "epoch"))
  dtm <- dtm[trial != "shallow_gam_mul"]
  dtm[, dataset := ifelse(grepl("^val_", variable), "validation", "train")]
  dtm[, metric := gsub("^val_", "", variable)]
  dtm <- separate(dtm, trial, c("depth", "type"))

  dtm_mean = dtm[, .(value = mean(value, na.rm=TRUE)),  by = .(metric, depth, type, epoch, dataset)]
  dtm_mean[, is_mean := TRUE]
  dtm[, is_mean := FALSE]
  dtm_w_mean = rbindlist(list(dtm_mean, dtm), use.names = TRUE, fill = TRUE)
  dtm_w_mean[, tid := ifelse(is.na(tid), 0, tid)]

  dtm_w_mean[, type := fct_relevel(type, "relu")]
  return(dtm_w_mean)
}
