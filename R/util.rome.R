error.bars <- function(x,upper,lower,width=0.02,...){
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}

plot_cv <- function(fit,...){
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
  
  nv <- predict_ecd(fit$fit, lambda = fit$lambda, type = "nvars")
  axis(3, at=l, labels=nv, tick=FALSE, line=-0.5)
}

coef_rome <- function(object, lambda, exact = FALSE, ...) {
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