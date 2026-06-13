#' Plot the cross-validation curve from a \code{"cv.rome"} object.
#'
#' @param x Fitted \code{"cv.rome"} object.
#' @param nvars Logical; if \code{TRUE} (default), show number of nonzero
#'    coefficients on the top axis.
#' @param ... Additional graphical arguments passed to \code{plot}.
#'
#' @method plot cv.rome
#' @importFrom graphics abline axis points segments
#' @export
plot.cv.rome <- function(x, nvars = TRUE, ...) {
  
  fit  <- x
  l    <- log(fit$lambda)
  cve  <- fit$cve
  L.cv <- cve - fit$cvse
  U.cv <- cve + fit$cvse
  
  ylab <- switch(fit$type.measure,
                 deviance = "Deviance",
                 mse      = "Mean Squared Error",
                 mae      = "Mean Absolute Error"
  )
  ylim <- range(c(L.cv, U.cv))
  
  plot.args <- list(
    x    = l, y = cve, ylim = ylim,
    xlab = expression(log(lambda)), ylab = ylab,
    type = "n", xlim = rev(range(l)), las = 1
  )
  new.args <- list(...)
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  
  points(l, cve, col = "red", pch = 19, cex = 0.5)
  error.bars(l, U.cv, L.cv, width = 0.01, col = "darkgrey")
  abline(v = log(fit$lambda_min), lty = 3)
  abline(v = log(fit$lambda_1se), lty = 3)
  
  if (nvars) {
    nv <- predict.rome(fit$fit, lambda = fit$lambda, type = "nvars")
    axis(3, at = l, labels = nv, tick = FALSE, line = -0.5)
  }
  
  invisible(fit)
}


#' Extract coefficients from a rome object
#'
#' @param object Fitted \code{"rome"} object.
#' @param lambda Optional values of penalty parameter \code{lambda} at which coefficients are requested.
#' @param exact Logical. If \code{TRUE}, exact coefficients are returned by refitting if needed.
#' @param ... Additional arguments passed to methods.
#'
#' @method coef rome
#' @export
coef.rome <- function(object, lambda, exact = FALSE, ...) {
  if (missing(lambda)) {
    beta <- object$beta
  } else if (exact) {
    ls <- object$lambda
    ind <- match(lambda, ls, 0)
    if (any(ind == 0)) {
      ls <- unique(rev(sort(c(lambda,ls))))
      object <- update(object, lambda=ls)
      ind <- match(lambda, ls)
    }
    beta <- object$beta[, ind]
  } else {
    ls <- object$lambda
    lambda[lambda>max(ls)] <- max(ls)
    lambda[lambda<min(ls)] <- min(ls)
    ind <- approx(ls, seq(ls), lambda)$y
    left <- floor(ind)
    right <- ceiling(ind)
    weight <- ind %% 1
    beta <- (1-weight)*object$beta[,left] + weight*object$beta[,right]
    if (length(lambda) > 1) {
      colnames(beta) <- round(lambda, 4)
    }
  }
  return(beta)
}


lambda.interp <- function(lambda, s) {
  k <- length(lambda)
  if (k == 1L) {
    nums  <- length(s)
    left  <- rep(1L, nums)
    right <- left
    sfrac <- rep(1, nums)
  } else {
    sfrac  <- (lambda[1] - s) / (lambda[1] - lambda[k])
    lam    <- (lambda[1] - lambda) / (lambda[1] - lambda[k])
    sfrac  <- pmin(pmax(sfrac, min(lam)), max(lam))
    coord  <- approx(lam, seq_along(lam), sfrac)$y
    left   <- floor(coord)
    right  <- ceiling(coord)
    sfrac  <- (sfrac - lam[right]) / (lam[left] - lam[right])
    sfrac[left == right] <- 1
    sfrac[abs(lam[left] - lam[right]) < .Machine$double.eps] <- 1
  }
  list(left = left, right = right, frac = sfrac)
}


check.dots <- function(object, ..., 
                       need = c("x", "y", "weights", "offset", "penalty.factor",
                                "lower.limits", "upper.limits"), 
                       error_start = "used coef.rome() or predict.rome() with `exact=TRUE`",
                       error_end = " in order to safely rerun rome",
                       prefix = NULL) {
  if (is.null(need)) {
    return(invisible())
  }
  
  # extract the function options we need from the object's call
  thiscall <- object$call
  ncall <- names(thiscall)[-1]
  w <- match(ncall, need, 0)
  need <- need[w]
  if (length(need) == 0) {
    return(invisible())
  }
  
  # check that ... indeed has those function options
  if (!is.null(prefix)) {
    need <- paste0(prefix, need)
  }
  nargs <- names(list(...))
  w <- match(need, nargs, 0) > 0
  if(!all(w)) {
    margs <- need[!w]
    stop(paste(error_start, 
               "so must in addition supply original argument(s) ",
               paste(margs,collapse=" and "), 
               error_end), call.=FALSE)
  }
  invisible()
}

error.bars <- function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x,        upper, x,        lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  invisible(range(upper, lower))
}
