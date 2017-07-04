## Evaluation metrics

cem_auc <- function(y_true, y_pred) {
  roc <- PRROC::roc.curve(y_pred[y_true==1], y_pred[y_true==0], curve = FALSE)
  return(roc$auc)
}


cem_auprc <- function(y_true, y_pred) {
  pr <- PRROC::pr.curve(y_pred[y_true==1], y_pred[y_true==0], curve = FALSE)
  return(pr$auc.integral)
}

##' Compute variance explained
##'
##' @param String. Which method to use. Can be "var" or "mad".
metric_variance_explained <- function(y, y_pred, method='var'){
  resid <- y - y_pred
  if(method == 'var'){
    # residual variance
    var_resid <- var(resid, na.rm=TRUE)
    # y response variance
    var_y <- var(y, na.rm=TRUE)
  }else{
    if(method == 'mad'){
      var_resid <- mad(resid, na.rm=TRUE)^2
      var_y <- mad(y, na.rm=TRUE)^2
    }
  }
  explained <- 1 - var_resid / var_y
  return(explained)
}
variance_explained <- metric_variance_explained

##' Bootstrap a variable
bootstrap_metric <- function(y_true, y_pred, metric_fn, n_bootstrap = 1000,
                              return_summary=TRUE,
                              ...)  {
  stopifnot(length(y_true) == length(y_pred), !(isTRUE(is.list(metric_fn)) && return_summary))
  m <- sapply(1:n_bootstrap, function(i) {
    ivec <- sample(1:length(y_true), replace = TRUE)
    if (is.list(metric_fn)) {
      return(as.data.table(lapply(metric_fn, function(fn) fn(y_true[ivec], y_pred[ivec], ...))))
    } else {
      return(metric_fn(y_true[ivec], y_pred[ivec], ...))
    }
  }, simplify=!is.list(metric_fn))
  if (is.list(metric_fn)) {
    return(rbindlist(m))
  }

  if (isTRUE(return_summary)) {
    return(list(mean = mean(m), sd = sd(m)))
  } else {
    return(m)
  }
}

bootstrap_measure <- bootstrap_metric

##'
##' plot the accuracy
get_performance <- function(accuracy) {
  lapply(accuracy, function(ac) {
    dt_train <- data.table(step = ac$step_history, loss = ac$train_acc_history, mode = "train")
    dt_test <- data.table(step = ac$step_history, loss = ac$val_acc_history, mode = "test")
    dt <- rbindlist(list(dt_train, dt_test))
    return(dt)
  }) %>% df_list2dt("fold")
}

##' Plot the accuracy
##' @param dc Deepcis in cross-validation
##' @param ylim - ylim of the curve
plot_accuracy <- function(dc, ylim = NULL) {
  accuracy <- get_accuracy(dc)
  dt <- get_performance(accuracy)  
  loss <- dt[, loss]
  if (is.null(ylim)) ylim <- c(min(loss), quantile(loss, 0.9))
  qplot(step, loss, color = mode, data = dt, facets = fold~., xlab = "epoch", ylab = "mse") + geom_line() + coord_cartesian(ylim=ylim) + 
    ggtitle("Accuracy plot")
}



## qplot(y = t(ac$loss.history), xlab = "step", ylab = "loss on the training batch", alpha = I(0.3))

mse2ermse <- function(mse) {
  10^sqrt(mse)
}



##' Convert the mean squared error to the coefficient of variation
##'
##' @param mse Mean squared error
##' @param base Logarithm base on which the mse is measured. Typically 10 or exp(1)
##' @return Coefficient of variation.
mse2cv <- function(mse, base) {
  sqrt(exp(log(base) * mse) - 1)
}

