## get gam splines from mgcv

##' Get the basis for B-splines using the mgcv package.
##'
##' It uses equal spacing to distribute the splines.
##'
##' @param x A vector of bases at which we want to evaluate the splines.
##' @param n_bases Number of basis functions created
##' @param spline_order Spline order. 3 for cubic, 2 for quadratic.
##' @param add_intercepd Should we add the intercept term to the design matrix + penalty matrix?
##' @return A list of elements:
##' - X  design matrix
##' - S  penalty matrix (I think mgcv takes just some differences and doesn't compute the integral)
##' - knots Positions of the spline knots.
##'
##' @examples
##' from = 0
##' to = 100
##' obj <- get_gam_splines(x = from:to, n_bases = 20, spline_order = 2, add_intercept = TRUE)
##' matplot(from:to, obj$X, type = "l")
get_gam_splines <- function(x = 0:10, n_bases = 10, spline_order = 3, add_intercept = TRUE) {

  dummy_dt = data.frame(x = x)
  sm <- mgcv::smoothCon(mgcv::s(x, k=n_bases, bs = "ps", m = spline_order - 1), data=dummy_dt, scale.penalty = FALSE, knots=NULL)[[1]]

  X <- sm$X
  S <- sm$S[[1]]

  if (add_intercept == TRUE) {
    X <- cbind(1, X)
    ## don't penalize the intercept term
    S <- cbind(0, rbind(0, S))
  }
  
  return(list(
    X = X,
    S = S,
    knots = sm$knots,
    x = x
    ## args = as.list(match.call()) ## get the function call arguments
  ))
}
