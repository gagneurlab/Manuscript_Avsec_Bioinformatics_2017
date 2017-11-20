#'---
#' title: Simple glmnet model
#' author: Žiga Avsec
#' wb:
#'   input: ["data/Splice_branchpoints/processed/branchpointer/train/filteredDescr.csv",
#'            "data/Splice_branchpoints/processed/branchpointer/test/filteredDescr.csv"]
#'   output: ["data/Splice_branchpoints/test_predictions/glmnet.csv"]
#'---
dt <- fread("data/Splice_branchpoints/processed/branchpointer/train/filteredDescr.csv")
dt_test <- fread("data/Splice_branchpoints/processed/branchpointer/train/filteredDescr.csv")

dt[, V1 := NULL]
dt_test[, V1 := NULL]

dt %>% head(3)

x2splines <- function(dt, spline_features, n_bases) {
  lapply(spline_features, function(feat) {
    if (is.matrix(dt)) {
      xvec <- dt[, spline_features]
    } else {
      xvec <- dt[[spline_features]]
    }

    gs <- get_gam_splines(xvec, n_bases=n_bases)
    X <- gs$X
    print(paste0("n_bases: ", n_bases))
    print(paste0("n_cols: ", ncol(X)))
    colnames(X) <- paste0(feat, 1:n_bases)
    X <- X %>% Matrix(sparse=TRUE)
    return(X)
  }) %>% Reduce(cBind, .)
}


branchpoint_data <- function(use_positions=TRUE, n_bases=10) {
  ## TODO - update
  dt <- fread("localdata/Splice_branchpoints/processed/branchpointer/train/filteredDescr.csv")
  dt_test <- fread("localdata/Splice_branchpoints/processed/branchpointer/train/filteredDescr.csv")
  dt[, V1 := NULL]
  dt_test[, V1 := NULL]

  create_x <- function(x, use_positions, n_bases) {
    x_seq <- x[, grepl("^seq", colnames(x))]
    if (!isTRUE(use_positions)) {
      return(x_seq)
    }
    pos_features <- colnames(x)[grepl("^seq", colnames(x))]
    x_spl <- x2splines(x, pos_features, n_bases)
    return(cBind(x_seq, x_spl))
  }
  x_train <- model.matrix(set~. - 1, dt)
  x_test <- model.matrix(set~.-1, dt_test)
  y_train <- dt[, as.integer(set == "HC")]
  y_test <- dt_test[, as.integer(set == "HC")]

  x_train <- create_x(x_train, use_positions, n_bases)
  x_test <- create_x(x_test, use_positions, n_bases)    

  ## center and scale - shouldn't do for the other features
  pp <- preProcess(x_train, method = c("center", "scale"))
  x_train <- predict(pp, x_train)
  x_test <- predict(pp, x_test)

  x_train <- as(x_train, "sparseMatrix")
  x_test <- as(x_test, "sparseMatrix")


  return(
    list(train = list(x = x_train,
                      y = y_train),
         test = list(x = x_test,
                     y = y_test)
       ))
}

x_train <- model.matrix(set~. - 1, dt)
x_test <- model.matrix(set~.-1, dt_test)
y_train <- dt[, as.factor(set)]
y_test <- dt_test[, as.factor(set)]

## center and scale - shouldn't do for the other features
pp <- preProcess(x_train, method = c("center", "scale"))
x_train <- predict(pp, x_train)
x_test <- predict(pp, x_test)

x_train <- as(x_train, "sparseMatrix")
x_test <- as(x_test, "sparseMatrix")

library(doMC)
registerDoMC(cores=5)
fit <- cv.glmnet(x_train, as.integer(y_train == "HC"),
                 family = "binomial", alpha = .5, nfolds=5)
plot(fit)
## lambda=0 if for sure the right choice (large dataset)
fit <- glmnet(x_train, as.integer(y_train == "HC"), family = "binomial", alpha = .5, lambda = 0)
plot(fit)
y_test_pred <- predict(fit, x_test)[,1]
## AUC is .90
confusionMatrix(as.integer(y_test_pred>0), as.integer(y_test=="HC"))

## get the ID
dt_pred <- data.table(y_true = y_test, y_pred = y_test_pred)
write_csv(dt_pred, "/s/project/deepcis/Splice_branchpoints/test_predictions/glmnet.csv")

