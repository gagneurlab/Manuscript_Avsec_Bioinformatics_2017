cv.glmnet_best_perf <- function(gfit) {
  i_min <- which(gfit$lambda == gfit$lambda.min)
  rnames <- c(paste0(gfit$name), paste0(gfit$name,"_sd"))
  result <- c(gfit$cvm[i_min],gfit$cvm[i_min])
  names(result) <- rnames
  return(result)
}

cv.glmnet_nonzero_coef <- function(gfit) {
  lasso_coef <- coef.cv.glmnet(gfit)
  bcoef <- lasso_coef[lasso_coef[,1] != 0,]
  return(bcoef[order(-bcoef)])
}


tidy.glmnet <- function(fit) {
  fit %>% coef %>% as.matrix %>% data.table(keep.rownames = TRUE) %>% setnames(c("term", "beta"))
}

tidy.cv.glmnet <- tidy.glmnet

##'
##' Jun's function
### Compute explained variance of linear model from the fit result
# method: "var", "mad"
variance_explained_lm_fit <- function(fit, method='var'){
  y_pred <- predict(fit)
  variance_explained(residuals(fit) + y_pred, y_pred, method = method)
}

variance_explained_glmnet_no_cv <- function(fit, x, y, method='var'){
  variance_explained(y, predict(fit, newx = x, s = "lambda.min")[, 1], method = method)
}

## fit needs to be run with keep = TRUE
glmnet_predict_dt <- function(fit, x, y) {
  ## 1. get out of fold predictions
  lambda <- fit$lambda.min
  alpha <- fit$glmnet.fit$call$alpha
  id <- 1:nrow(x)
  if (is.null(alpha)) alpha <- 1
  ## fit in CV - returns y, and y_pred
  ## generate folds
  fids <- fit$foldid
  if (is.null(fids)) stop("no foldid information. glmnet has to be trained with argument: keep = TRUE")
  dt <- lapply(unique(fids), function(i) {
    ## train/test set

    X_train <- x[fids != i,]
    X_test <- x[fids == i,]
    y_train <- y[fids != i]
    y_test <- y[fids == i]
    id_test <- id[fids == i]

    ## glmnet fit&predict
    fit <- glmnet(x = X_train, y = y_train, alpha = alpha, lambda = lambda)
    y_pred <- predict(fit, newx = X_test)

    return(data.table(y_test = y_test,
                      y_pred = y_pred[, 1],
                      fold_id = i,
                      row_id = id_test))
  }) %>% rbindlist  %>% .[order(row_id)] %>% .[, row_id := NULL]
}

variance_explained_glmnet <- function(fit, x, y, method = 'var', by_folds = FALSE) {
  dt <- glmnet_predict_dt(fit, x, y)
  if (isTRUE(by_folds)) {
    dt <- dt[, .(variance_explained = variance_explained(y_test, y_pred, method = method)), by = fold_id]
    dt <- dt[, mean_variance_explained := mean(variance_explained)]
    return(dt)
  } else {
  return(variance_explained(dt$y_test, dt$y_pred, method = method))
  }
}

mse_glmnet <- function(fit, x, y, by_folds = FALSE) {
  dt <- glmnet_predict_dt(fit, x, y)
  if (isTRUE(by_folds)) {
    dt <- dt[, .(mse = mse(y_test, y_pred)), by = fold_id]
    dt <- dt[, mean_mse := mean(mse)]
    return(dt)
  } else {
  return(mse(dt$y_test, dt$y_pred))
  }
}


fit_glmnet_best_perf <- function(x, y, alpha, ...) {
  glmnet(x, y,
         alpha = alpha,
         lambda=cv.glmnet(x, y, alpha = alpha, ...)$lambda.min)  
}

