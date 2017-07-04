#'---
#' title: Trials
#'---

#   input: ["data/Concise/Splice_branchpoints/trials/df/deep_dam.csv"]
# 

get_trials_df_eval <- function(all_trials = c("deep_gam", "deep_relu",
                                              "shallow_gam", "shallow_relu")) {
  dt_eval <- lapply(all_trials, function(trial) {
      EXP_DIR <- "data/Concise/Splice_branchpoints/"
      dt <- fread(file.path(EXP_DIR, "/trials/df", paste0(trial, ".csv")))
      dt[, trial := trial]
      dt[, V1 := NULL]
      dt <- dt[, c("tid", "trial",
                   grep("^eval", colnames(dt), value=TRUE),
                   grep("^time", colnames(dt), value=TRUE)),
               with = F]
      return(dt)
  }) %>% rbindlist(use.names=TRUE, fill=TRUE) %>% separate(trial, c("depth", "type"))
  return(dt_eval)
}
