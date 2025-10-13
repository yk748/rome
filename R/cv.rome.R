#' Cross-validation a penalized robust regression
#'
#' Does k-fold cross-validation for rome, produces a plot, and returns a
#' value for \code{lambda}
#' 
#' @param y response variable.
#' @param x input matrix, of dimension nobs x nvars; 
#' @param adapt a logical flag whether the adaptive (weighting) scheme is used. Default is \code{FALSE}.
#' @param weights observation weights. Can be total counts if responses are
#' proportion matrices. Default is 1 for each observation
#' @param delta threshold for Huber loss function. The default is 0.5. 
#' @param nlambda Length of a sequence of \code{lambda}. The default length is 100.
#' @param lambda.min A ratio for minimum \code{lambda} from its maximum value. The default is 0.05.
#' @param lambda a sequence of \code{lambda}. It can be provided by users.
#' @param nfolds number of folds. Default is 10.
#' @param foldid an optional vector of values between 1 and \code{nfolds} identifying what fold each observation is in. If supplied, \code{nfolds} can be missing.
#' @param preprocess a logical flag whether the data should be preprocessed. Currently, only \code{none} (no preprocessing) is provided.
#' @param centering a logical flag whether the data should be centered. If so, the data is shifted by their mean levels. The default is \code{FALSE}.
#' @param max.iter a maximum iteration for each update of a variable at a fixed \code{lambda}. The default is 100. 
#' @param eps tolerance for convergence, used for stopping criteria. If the screen rules are used, it will be used for scaling the null deviance. If no screen rule is used, each update of a variable is terminated when a ell2 norm for difference between two solutions of consecutive iteration is less than tolerance. The default is 1e-4.
#' @param intercept a logical flag whether the intercept is added. If so, the dimension of variable increases by 1. The default is \code{FALSE}.
#' @param screen a logical flag whether the adaptive strong screen rule, strong screen rule, or no screen rule (\code{none}) is applied.
#' @param dfmax Upper bound for the number of nonzero coefficients. The algorithm exits and 
#' returns a partial path if \code{dfmax} is reached. Useful for very large dimensions.
#' @param trace Logical flag to ask the message for each progress. Default is \code{FALSE}.
#' @author Younghoon Kim \cr Maintainer: Younghoon Kim
#' \email{yk748@cornell.edu}
#' @references Kim, Y. Loh, PL. S. Basu (2025)
#' \emph{Exact Coordinate Descent for High-Dimensional Regularized Huber Regression, Preprint}
#' @keywords models regression
#' @export cv.rome
#' 

cv.rome <- function(x,y, FUN=c("rome.adaptive"),
                    ncores = 1, nfolds=10, foldid=NULL, 
                    type.measure = c("deviance", "mse", "mae"),...){
  
  FUN <- get(match.arg(FUN))
  type.measure <- match.arg(type.measure)
  n <- length(y)
  if (is.null(foldid)){
    foldid <- ceiling(sample(1:n)/n*nfolds)
  } else {
    nfolds <- max(foldid)
  }
  
  fit <- FUN(y=y, x=x, ...) 
  
  cv.args <- list()
  cv.args$lambda <- fit$lambda 
  cv.args$delta <- fit$delta
  measure.args <- list(delta=fit$delta, type.measure = type.measure)
  
  parallel <- FALSE
  if (ncores > 1) {
    max.cores <- detectCores()
    if (ncores > max.cores) {
      cat("The number of cores specified (", ncores, ") is larger than
          the number of avaiable cores (", max.cores, "), so", 
          max.cores, "cores are used.", "\n")
      ncores <- max.cores
    }
    cluster <- makeCluster(ncores)
    if (!("cluster" %in% class(cluster))) {
      stop("Cluster is not of class 'cluster'; see ?makeCluster")
    }
    parallel <- TRUE
    cat("Start parallel computing for cross-validation...")
    clusterExport(cluster, c("foldid", "y", "x", "cv.args", "measure.args"),
                  envir=environment())
    clusterCall(cluster, function() require(rome))
    fold.results <- parLapply(cl = cluster, X = 1:nfolds, fun=cv.ftn, y=y, x=x,
                              foldid=foldid, cv.args=cv.args, measure.args=measure.args)
    stopCluster(cluster)
  }
  
  E <- matrix(NA, nrow=length(y), ncol=length(cv.args$lambda))
  for (i in 1:nfolds) {
    if (parallel) {
      fit_i_fold <- fold.results[[i]]
    } else {
      cat("CV fold #",i,"begins\n")
      fit_i_fold <- cv.ftn(i,x,y,foldid,cv.args,measure.args,FUN)
    }
    E[foldid == i, 1:fit_i_fold$nl] <- fit_i_fold$pe
  }
  
  
  idx <- which(apply(is.finite(E), 2, all))
  E <- E[,idx]
  lambda <- fit$lambda
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, sd) / sqrt(n)
  idx_min <- which.min(cve)
  idx_1se <- min(which(cve < cve[idx_min]+cvse[idx_min]))
  val <- list(cve = cve, cvse = cvse, type.measure = type.measure,
              lambda = lambda, fit=fit, 
              lambda_1se = lambda[idx_1se], lambda_min = lambda[idx_min])
  structure(val, class="cv.rome")
}


cv.ftn <- function(i,x,y,foldid,cv.args,measure.args,FUN) {
  cv.args$x <- x[foldid != i,,drop = FALSE]
  cv.args$y <- y[foldid != i]
  x_ts <- x[foldid == i,,drop = FALSE]
  y_ts <- y[foldid == i]
  
  # fit_i <- do.call(FUN, cv.args)
  fit_i <- do.call("rome.adaptive", cv.args)
  
  y_hat <- matrix(predict.rome(fit_i,x_ts,type="response"),length(y_ts))
  
  list(pe = measure.rome(y_ts, y_hat, measure.args), 
       nl = length(fit_i$lambda))
}

measure.rome <- function(y, y_hat, measure.args) {
  
  r <- y - y_hat
  type.measure <- measure.args$type.measure
  
  if (type.measure == "deviance") {
    
    delta <- measure.args$delta
    val <- huber.loss(r, delta)
    
  } else if (type.measure == "mse") {
    val <- r^2
  } else {
    val <- abs(r)
  }
  return(val)
}


huber.loss <- function(v, thres){
  return( ifelse(abs(v) <= thres,v^2/2,thres*abs(v) - thres^2/2) )
}

