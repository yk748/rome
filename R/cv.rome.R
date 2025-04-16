#' Cross-validation a penalized robust regression
#'
#' Does k-fold cross-validation for rome, produces a plot, and returns a
#' value for \code{lambda}
#' 
#' @param y response variable.
#' @param x input matrix, of dimension nobs x nvars; 
#' @param weights observation weights. Can be total counts if responses are
#' proportion matrices. Default is 1 for each observation
#' @param delta threshold for Huber loss function. The default is 1.5. 
#' @param lambda a sequence of \code{lambda}. It can be provided by users. The default is \code{NULL}.
#' @param nfolds number of folds - default is 10.
#' @param foldid an optional vector of values between 1 and \code{nfolds} identifying what fold each observation is in. If supplied, \code{nfolds} can be missing.
#' @param preprocess a logical flag whether the data should be preprocessed. Currently, only \code{none} (no preprocessing) is provided.
#' @param centering a logical flag whether the data should be centered. If so, the data is shifted by their mean levels. The default is \code{FALSE}.
#' @param max.iter a maximum iteration for each update of a variable at a fixed \code{lambda}. The default is 100. 
#' @param tol tolerance for convergence, used for stopping criteria. The each update of a variable is terminated when a ell2 norm for difference between two solutions of consecutive iteration is less than tolerance. The default is 1e-6.
#' @param intercept a logical flag whether the intercept is added. If so, the dimension of variable increases by 1. The default is \code{FALSE}.
#' @author Younghoon Kim \cr Maintainer: Younghoon Kim
#' \email{yk748@cornell.edu}
#' @references Kim, Y. Loh, PL. S. Basu (2025)
#' \emph{Exact Coordinate Descent for High-Dimensional Regularized Robust M-Estimators, ??, Vol. ??(??), ??-??},
#' \doi{??}.\cr
#' @keywords models regression
#' @export cv.rome
#' 

cv.rome <- function(x,y,weights=NULL,delta=NULL,lambda=NULL,
                    nfolds=10,foldid=NULL,
                    preprocess = c("none"), 
                    centering = FALSE,
                    max.iter = 100, tolerance = 1e-6, 
                    intercept = FALSE){
  

  if (!is.null(lambda) && length(lambda) < 2){
    stop("Need more than one value of lambda for cv_rome")
  }
  
  if (is.null(foldid)){
    foldid <- ceiling(sample(1:n)/n*nfolds)
  } else {
    nfolds <- max(foldid)
  }
  
  fit <- rome(x, y, weights=weights, delta=delta, lambda=lambda, 
              max.iter=max.iter, tolerance=tolerance)
  
  n <- length(y)
  E <- matrix(NA, nrow = n, ncol = length(fit$lambda))
  for (i in 1:nfolds) {
    cat("CV fold #",i,"begins\n")
    fit_i_fold <- cv.ftn(i,x,y,foldid,delta=fit$delta)
    E[foldid == i, 1:fit_i_fold$nl] <- fit_i_fold$pe
  }
  
  idx <- which(apply(is.finite(E), 2, all))
  E <- E[,idx]
  lambda <- fit$lambda
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, sd) / sqrt(n)
  idx_min <- which.min(cve)
  idx_1se <- min(which(cve < cve[idx_min]+cvse[idx_min]))
  val <- list(cve = cve, cvse = cvse,
              lambda = lambda, fit=fit, 
              lambda_1se = lambda[idx_1se], lambda_min = lambda[idx_min])
  structure(val, class="cv.rome")
}


cv.ftn <- function(i,x,y,foldid,delta) {
  x_tr <- x[foldid != i,,drop = FALSE]
  y_tr <- y[foldid != i]
  x_ts <- x[foldid == i,,drop = FALSE]
  y_ts <- y[foldid == i]
  
  fit_i <- rome(x=x_tr,y=y_tr,delta=delta)
  
  y_hat <- matrix(predict.rome(fit_i,newx=x_ts,type="response"),length(y_ts))
  
  # Currently, only deviance can be provided
  r_hat <- y_ts - y_hat
  val <- array(NA,dim=dim(r_hat))
  for (r in 1:dim(r_hat)[1]){
    for (s in 1:dim(r_hat)[2]){
      val[r,s] <- ifelse(r_hat[r,s] <= delta,
                         r_hat[r,s]^2/2,
                         delta*abs(r_hat[r,s]) - delta^2/2)
    }
  }
  return(list(pe=val,nl=length(fit_i$lambda)))
}
