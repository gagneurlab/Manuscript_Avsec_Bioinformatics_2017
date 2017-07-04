
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



## TODO - implement
## mse_na.omit <- function(y, y_pred) {
##   dto <- data.table(y = c(y), y_pred = c(y_pred)) %>% na.omit
##   mse(dto$y, dto$y_pred)
## }

##' Compute variance explained
##'
##' @param String. Which method to use. Can be "var" or "mad".
variance_explained <- function(y, y_pred, method='var'){
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


##' Bootstrap a variable
bootstrap_measure <- function(y_true, y_pred, measure_fun, n_bootstrap = 1000,
                              return_summary=TRUE,
                              ...)  {
  stopifnot(length(y_true) == length(y_pred))
  m <- sapply(1:n_bootstrap, function(i) {
    ivec <- sample(1:length(y_true), replace = TRUE)
    return(measure_fun(y_true[ivec], y_pred[ivec], ...))
  })
  if (isTRUE(return_summary)) {
    return(list(mean = mean(m), sd = sd(m)))
  } else {
    return(m)
  }
}
