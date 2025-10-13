#' Make predictions from a "rome" object.
#'
#' Similar to other predict methods, this function predicts fitted values,
#' coefficients, and number of non-zero coefficients from a fitted \code{"rome"} object.
#' @param object Fitted \code{"rome"} model object.
#' @param X Matrix of new values for \code{x} at which predictions are to be made.
#' @param lambda Value(s) of the penalty parameter \code{lambda} at which predictions are required.
#' @param type Type of prediction: "response", "coefficients", or "nvars".
#' @param ... Additional arguments passed to underlying methods.
#' @method predict cv.rome
#' @exportS3Method predict cv.rome
predict.cv.rome <- function(object, 
                         X, lambda = c("lambda.1se","lambda.min"), 
                         type=c("response","coefficients","nvars"),...) {
  type = match.arg(type)
  if (is.character(lambda)) {
    lambda = match.arg(lambda)
    lambda = object[[lambda]]
  } else if(!is.numeric(lambda)) {
    stop("Invalid lambda")
  }
  predict(object$fit, X, lambda = lambda, type = type, ...)
}

#' Extract coefficients from a "cv.rome" object.
#'
#' @param object Fitted \code{"cv.rome"} model object.
#' @param lambda Value(s) of the penalty parameter \code{lambda}.
#' @param ... Additional arguments passed to underlying methods.
#' @method coef cv.rome
#' @exportS3Method coef cv.rome
coef.cv.rome <- function(object, lambda=c("lambda.lse","lambda.min"),...) {
  if (is.character(lambda)) {
    lambda = match.arg(lambda)
    lambda = object[[lambda]]
  } else if(!is.numeric(lambda)) {
    stop("Invalid lambda")
  }
  coef(object$fit, lambda = lambda, ...)
}

