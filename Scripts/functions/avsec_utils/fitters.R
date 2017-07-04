##' fit a linear model in cross-validation
##' @param formula lm formula
##' @param data containing all the variables specified in formula
##' @param folds nested list of test, train indicies
fit_lm_cv <- function(formula, data, folds, ...) {

  ## run the cross-validation
  cv_model <- lapply(folds, function(fold) {
    X_train <- data[fold$train,]
    rownames(X_train) <- rownames(data)[fold$train]
    X_test <- data[fold$test,]
    rownames(X_test) <- rownames(data)[fold$test]
    
    fit <- lm(formula, data = X_train)
    y_predict <- predict(fit, newdata = X_test)
    y <- model.response(model.frame(formula, data = X_test))
    return(list(
      mse_test_acc = mse(y_predict, y),
      coef = coef(fit),
      fit = fit,
      id = names(y_predict),
      y = y,
      y_predict = y_predict
    ))
  })

  ## add fold-names
  names(cv_model) <- paste0("fold.", 1:length(cv_model))

  final_model_fit <- lm(formula = formula, data = data, ...)
  final_model <- list(
    mse_train_acc = mse(predict(final_model_fit), data$hlt),
    coef = coef(final_model_fit),
    fit = final_model_fit
  )

  return(list(cv_model = cv_model,
              final_model = final_model))
}

extract_cv_performance <- function(fit_lm_cv_obj) {
  mse_cv <- sapply(fit_lm_cv_obj$cv_model, function(x) x$mse_test_acc)
  mse_mean <- mean(mse_cv)
  mse_sd <- sd(mse_cv)
  return(c("mean" = mse_mean, "sd" = mse_sd))
}

get_predict_table_lm <- function(lm_fit) {
  y <- model.response(lm_fit$model)
  y_predict <- predict(lm_fit)
  id <- names(y_predict)

  dt <- data.table(id = id,
                   y = y,
                   y_predict = y_predict
                   )
  return(dt)
}

##' get the prediction table from fit_lm_cv_obj in cross validation
##' @return data.table with features:
##' id - row name
##' y - response variable
##' y_predict - out-of-fold prediction of y
##' y_predict_in_fold - prediction of y without using cross-validation
get_cv_predict_table_lm <- function(fit_lm_cv_obj) {
  
  if (all(is.null(names(fit_lm_cv_obj$cv_model)))) {
    names(fit_lm_cv_obj$cv_model) <- paste0("fold.", 1:length(fit_lm_cv_obj$cv_model))
  }

  
  dt <- rbindlist(lapply(names(fit_lm_cv_obj$cv_model), function(fold) {
    x <- fit_lm_cv_obj$cv_model[[fold]]
    dt <- as.data.table(x[c("y", "y_predict", "id")])
    dt[, fold := fold]
    return(dt)
  }
    ))

  dt_overall <- get_predict_table_lm(fit_lm_cv_obj$final_model$fit)
  dt_overall[, y := NULL]
  setnames(dt_overall, "y_predict", "y_predict_in_fold")
  dt_final <- merge(dt, dt_overall, by = "id", all = TRUE)

  stopifnot(!any(is.na(dt_final)))

  return(dt_final)
}

############################################
## helpers

##' Compute the mean-squared error between two vectors
mse <- function(x,y) {
  return(mean((x-y)^2, na.rm = TRUE))
}

