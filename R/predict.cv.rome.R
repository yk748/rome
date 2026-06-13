#' Make predictions from a \code{"cv.rome"} object.
#'
#' @param object Fitted \code{"cv.rome"} object.
#' @param X Matrix of new predictor values.
#' @param lambda Either a numeric value or one of \code{"lambda.min"} or
#'   \code{"lambda.1se"}. Default is \code{"lambda.1se"}.
#' @param type Type of prediction: \code{"response"}, \code{"coefficients"},
#'   or \code{"nvars"}.
#' @param ... Additional arguments passed to \code{predict.rome}.
#' @method predict cv.rome
#' @exportS3Method predict cv.rome
predict.cv.rome <- function(object,
                            X,
                            lambda = c("lambda.1se", "lambda.min"),
                            type   = c("response", "coefficients", "nvars"),
                            ...) {
  type <- match.arg(type)
  
  if (is.character(lambda)) {
    lambda <- match.arg(lambda)
    lambda <- object[[lambda]]
  } else if (!is.numeric(lambda)) {
    stop("'lambda' must be numeric or one of 'lambda.min' / 'lambda.1se'.")
  }
  
  predict.rome(object$fit, X = X, lambda = lambda, type = type, ...)
}

#' Extract coefficients from a \code{"cv.rome"} object.
#'
#' @param object Fitted \code{"cv.rome"} object.
#' @param lambda Either a numeric value or one of \code{"lambda.min"} or
#'   \code{"lambda.1se"}. Default is \code{"lambda.1se"}.
#' @param ... Additional arguments passed to \code{coef.rome}.
#' @method coef cv.rome
#' @exportS3Method coef cv.rome
coef.cv.rome <- function(object,
                         lambda = c("lambda.1se", "lambda.min"),
                         ...) {
  if (is.character(lambda)) {
    lambda <- match.arg(lambda)
    lambda <- object[[lambda]]
  } else if (!is.numeric(lambda)) {
    stop("'lambda' must be numeric or one of 'lambda.min' / 'lambda.1se'.")
  }
  
  coef.rome(object$fit, lambda = lambda, ...)
}
