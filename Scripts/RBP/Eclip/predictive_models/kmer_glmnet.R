#!/usr/bin/env Rscript
#'---
#' title: kmer_glmnet
#'---
##' usage script.R rbp_name

## train glmnet k-mer model

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) stop("One argument required")
rbp_name <- args[1]
position_type <- args[2]

## parse positional argument
if (position_type == "w_position") {
  with_position <- TRUE
} else if (position_type == "no_position") {
  with_position <- FALSE
} else {
  stop("position_type unrecognized")
}


## dummy variables
## rbp_name <- c("UPF1")
## k <- 3
## alpha <- 0

flog.info(paste0("Read in all the data. RBP: ", rbp_name))
eclip_load_modelling_data(rbp_name)     #get all required data


## add position
if (isTRUE(with_position)) {
  flog.info(paste0("Generate X_TSS and X_polya"))
  N_BASES <- 10
  dt_train[, dataset := "train"]
  dt_valid[, dataset := "valid"]
  dt_test[, dataset := "test"]
  dt <- rbindlist(list(dt_train, dt_valid, dt_test), use.names=TRUE)
  ## TSS
  dt[, TSS_distance := log10(TSS_distance + 1)]
  stopifnot(all(!is.na(dt$TSS_distance)))
  gs_TSS <- get_gam_splines(dt$TSS_distance, n_bases=N_BASES)
  X_TSS <- gs_TSS$X %>% Matrix(sparse=TRUE)
  colnames(X_TSS) <- paste0("TSS_spline", 1:N_BASES)
  ## polya
  dt[, polya_distance := log10(-polya_distance+1)]
  stopifnot(all(!is.na(dt$polya_distance)))
  gs_polya <- get_gam_splines(dt$polya_distance, n_bases=N_BASES)
  X_polya <- gs_polya$X %>% Matrix(sparse=TRUE)
  colnames(X_polya) <- paste0("polya_spline", 1:N_BASES)

  ## check that the order didn't change
  stopifnot(all(dt[dataset == "train", seq] == dt_train[, seq]),
            all(dt[dataset == "valid", seq] == dt_valid[, seq]),
            all(dt[dataset == "test", seq] == dt_test[, seq]))
}

## CONFIG
this_script <- "kmer_glmnet.R"
k_vec <- 6:7                      #it breaks at 8-mers on 31 gigs of RAM...
## alpha_vec <- seq(.5, 1, by= 0.25)
alpha_vec <- 0.5
##----
flog.info("Run the model")
list_dt_kmer <- lapply(k_vec, function(k) {
  gc()
  flog.info("Generate kmers train")
  X_train <- dt_train[, get_kmers_X(get(info$sequence), k)] %>% Matrix(sparse = TRUE)
  flog.info("Generate kmers valid")
  X_valid <- dt_valid[, get_kmers_X(get(info$sequence), k)] %>% Matrix(sparse = TRUE)


  if (isTRUE(with_position)) {
    ## append positional columns
    X_train <- cbind(X_train, X_TSS[dt$dataset == "train",], X_polya[dt$dataset == "train",])
    X_valid <- cbind(X_valid, X_TSS[dt$dataset == "valid",], X_polya[dt$dataset == "valid",])
  }

  X_train_valid <- rbind(X_train, X_valid)

  flog.info("Start fitting")
  lapply(alpha_vec, function(alpha) {
    message("k = ", k, ", alpha = ", alpha)
    res <- tryCatch({
      tic()
      cv_fit <- cv.glmnet(X_train, y_train,
                          family = "binomial",
                          alpha = alpha,
                          standardize = FALSE, intercept = TRUE,
                          type.measure = "auc", nfolds = 10)
      auc_valid <- auc(y_valid, predict(cv_fit, newx = X_valid, s = "lambda.min", type = "response")[,1])
      time = toc()
      elapsed = time$toc - time$tic
      res <- data.table(auc= auc_valid,
                        auc_cv_train = cv.glmnet_best_perf(cv_fit)[1],
                        auc_cv_train_sd = cv.glmnet_best_perf(cv_fit)[2],
                        lambda = cv_fit$lambda.min,
                        alpha = alpha,
                        k = k,
                        training_time = elapsed
                        )
      gc()
      ## get global model
      overall_model <- glmnet(x = X_train_valid,
                              y = y_train_valid,
                              family = "binomial",
                              standardize = FALSE, intercept = TRUE,
                              alpha = res$alpha,
                              lambda = res$lambda)
      return(list(res = res, final_model = overall_model))
    }, error = function(e) {
      flog.error(e)
      return(NULL)
    }
    )
    return(res)
  })
}) %>% unlist(recursive = FALSE)

## remove NULL values
failed <- sapply(list_dt_kmer, is.null)
if (sum(failed)>0) flog.warn(paste0("Number of failed jobs: ",
                                       sum(failed), "/", length(failed)))
list_dt_kmer <- list_dt_kmer[!failed]
flog.info("Predict for test data")
## best model params
## remove NONE
dt_kmer <- list_dt_kmer %>% lapply(function(x) x$res) %>% rbindlist
best_model <- list_dt_kmer[[which.max(dt_kmer$auc)]]
model_obj <- best_model$final_model
param <- best_model$res %>% as.list %>% .[c("k", "alpha", "lambda")]
execution_time <- best_model$res$training_time
valid_performance <- best_model$res %>% as.list %>% .[["auc"]]
## test performance
X_test <- dt_test[, get_kmers_X(get(info$sequence), param$k)] %>% Matrix(sparse = TRUE) 
## add positional features
if (isTRUE(with_position)) {
  X_test <- cbind(X_test, X_TSS[dt$dataset == "test",], X_polya[dt$dataset == "test",])
}
y_test <- dt_test[[info$response]]
y_test_pred <- predict(model_obj, newx = X_test, s = "lambda.min", type = "response") %>% .[,1]
dt_pred_best <- data.table(y_true = y_test, y_pred = y_test_pred)
test_performance <- auc(as.integer(y_test), y_test_pred)
test_performance_roc <- roc(as.integer(y_test), y_test_pred)

all_kmers <- list(
  description = "All kmers (glmnet)",
  script = this_script,
  best_model = list(
    param = param,
    valid_performance = valid_performance,
    test_performance = test_performance,
    test_performance_roc = test_performance_roc,
    dt_pred = dt_pred_best,
    execution_time = execution_time,
    model_obj = model_obj
  ),
  validation_performance_table = dt_kmer
)

## save the model 
flog.info("Save")
saveRDS(all_kmers, eclip_modelling_get_output_file(this_script, rbp_name, position_type=position_type))
write_csv(dt_pred_best, eclip_modelling_get_output_file_csv(this_script,
                                                            rbp_name,
                                                            position_type))
flog.info("Done")
