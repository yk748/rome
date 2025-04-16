error.bars <- function(x,upper,lower,width=0.02,...){
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}

#' plot the cross-validation curve produced by cv.rome
#'
#' A plot is produced, and nothing is returned.
#'
#' @aliases plot.cv.rome
#' @param x fitted \code{"cv.rome"} object
#' @seealso \code{rome} and \code{cv.rome}.
#'
#' @method plot cv.rome
#' @export
plot.cv.rome <- function(fit,...){
  l <- fit$lambda
  l <- -log(l)
  xlab <- expression(-log(lambda))
  
  L.cve <- fit$cve - fit$cvse
  U.cve <- fit$cve + fit$cvse
  y <- fit$cve
  L <- L.cve
  U <- U.cve
  ylab <- "Deviance"
  
  ylim <- range(c(L, U))
  ind <- ((U-L)/diff(ylim) > 1e-3)
  plot.args <- list(x=l, y=y, ylim=ylim, xlab=xlab, ylab=ylab, 
                    type="n", xlim=rev(range(l)), las=1)
  new.args <- list(...)
  if (length(new.args)) {
    plot.args[names(new.args)] <- new.args
  }
  do.call("plot", plot.args)
  points(l, y, col="red", pch=19, cex=.5)
  beta <- fit$fit$beta
  
  error.bars(-log(fit$lambda),U.cve,L.cve,width=0.01,col="darkgrey")
  abline(v=log(fit$lambda_min),lty=3)
  abline(v=log(fit$lambda_1se),lty=3)
  
  nv <- predict.rome(fit$fit, lambda = fit$lambda, type = "nonzero")
  axis(3, at=l, labels=nv, tick=FALSE, line=-0.5)
}


#' Extract coefficients from a rome object
#'
#' @method coef rome
#' @rdname predict.rome
#' @export
#' @export coef.rome
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


lambda.interp <- function(lambda,s){
  
  if(length(lambda)==1){# degenerate case of only one lambda
    nums <- length(s)
    left <- rep(1,nums)
    right <- left
    sfrac <- rep(1,nums)
  }
  else{
    k <- length(lambda)
    sfrac <- (lambda[1]-s)/(lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
    sfrac[sfrac < min(lambda)] <- min(lambda)
    sfrac[sfrac > max(lambda)] <- max(lambda)
    coord <- approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sfrac-lambda[right])/(lambda[left] - lambda[right])
    sfrac[left==right] <- 1
    sfrac[abs(lambda[left]-lambda[right])<.Machine$double.eps] <- 1
    
  }
  list(left=left,right=right,frac=sfrac)
}


check_dots <- function(object, ..., 
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