#' Make predictions from a \code{"rome"} object.
#'
#' Predicts fitted values, coefficients, or number of nonzero coefficients
#' from a fitted \code{"rome"} object.
#'
#' @param object Fitted \code{"rome"} model object.
#' @param X Matrix of new values for \code{x}. Not used for
#'   \code{type = "coefficients"} or \code{type = "nvars"}.
#' @param lambda Value(s) of \code{lambda} at which predictions are required.
#'   Default is the full sequence used to fit the model.
#' @param type Type of prediction: \code{"response"}, \code{"coefficients"},
#'   or \code{"nvars"}.
#' @param exact If \code{TRUE}, refit the model at the supplied \code{lambda}
#'   values rather than interpolating. Requires the original \code{x} and
#'   \code{y} to be supplied via \code{...}.
#' @param ... Additional arguments (e.g. \code{x}, \code{y}) required when
#'   \code{exact = TRUE}.
#' @method predict rome
#' @export
predict.rome <- function(object,
                         X,
                         lambda,
                         type  = c("response", "coefficients", "nvars"),
                         exact = FALSE,
                         ...) {
  type <- match.arg(type)
  
  if (missing(X) && type == "response")
    stop("You need to supply a value for 'X' when type = 'response'.")
  
  beta <- coef.rome(object, lambda = lambda, exact = exact, ...)
  
  switch(type,
         coefficients = beta,
         nvars = if (is.matrix(beta)) apply(beta != 0, 2, sum) else sum(beta != 0),
         response     = X %*% beta
  )
}


#' Extract coefficients from a \code{"rome"} object.
#'
#' @param object Fitted \code{"rome"} model object.
#' @param lambda Value(s) of \code{lambda}. If missing, the full coefficient
#'   matrix is returned.
#' @param exact If \code{TRUE}, refit the model at the supplied \code{lambda}
#'   values. Requires original data via \code{...}.
#' @param ... Passed to the underlying fitter when \code{exact = TRUE}.
#' @method coef rome
#' @export
coef.rome <- function(object, lambda, exact = FALSE, ...) {
  
  if (missing(lambda)) return(object$beta)
  
  ls <- object$lambda
  
  if (exact) {
    ind <- match(lambda, ls, 0L)
    if (any(ind == 0L)) {
      ls     <- unique(sort(c(lambda, ls), decreasing = TRUE))
      object <- update(object, lambda = ls)
      ind    <- match(lambda, ls)
    }
    return(object$beta[, ind, drop = FALSE])
  }
  
  # Linear interpolation for lambda values not in the fitted sequence
  lambda <- pmin(pmax(lambda, min(ls)), max(ls))
  ind    <- approx(ls, seq_along(ls), lambda)$y
  left   <- floor(ind)
  right  <- ceiling(ind)
  weight <- ind %% 1
  
  beta <- (1 - weight) * object$beta[, left, drop = FALSE] +
    weight  * object$beta[, right, drop = FALSE]
  
  if (length(lambda) > 1) colnames(beta) <- round(lambda, 4)
  beta
}