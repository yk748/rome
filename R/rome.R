#' fit a penalized robust regression
#'
#' Fit a penalized robust regression (Huber loss) via exact coordinate descent algorithm.
#' 
#' @param y response variable.
#' @param x input matrix, of dimension nobs x nvars; 
#' @param weights observation weights. Can be total counts if responses are
#' proportion matrices. Default is 1 for each observation
#' @param delta threshold for Huber loss function. The default is 1.5. 
#' @param nlambda Length of a sequence of \code{lambda}. The default length is 100.
#' @param lambda.min A ratio for minimum \code{lambda} from its maximum value. The default is 0.05.
#' @lambda a sequence of \code{lambda}. It can be provided by users.
#' @preprocess a logical flag whether the data should be preprocessed. Currently, only \code{none} (no preprocessing) is provided.
#' @centering a logical flag whether the data should be centered. If so, the data is shifted by their mean levels. The default is \code{FALSE}.
#' @max.iter a maximum iteration for each update of a variable at a fixed \code{lambda}. The default is 100. 
#' @tol tolerance for convergence, used for stopping criteria. The each update of a variable is terminated when a ell2 norm for difference between two solutions of consecutive iteration is less than tolerance. The default is 1e-6.
#' @intercept a logical flag whether the intercept is added. If so, the dimension of variable increases by 1. The default is \code{FALSE}.
#' @author Younghoon Kim \cr Maintainer: Younghoon Kim
#' \email{yk748@cornell.edu}
#' @references Kim, Y. Loh, PL. S. Basu (2025)
#' \emph{Exact Coordinate Descent for High-Dimensional Regularized Robust M-Estimators, ??, Vol. ??(??), ??-??},
#' \doi{??}.\cr
#' @keywords models regression
#' @import utils
#' @import stats
#' @import mvtnorm
#' @useDynLib rome ecd_huber_
#' @export rome
#' 
#' @examples
rome <- function(x, y, method = c("huber"), weights=NULL,
                        delta = NULL, 
                        nlambda=100, lambda.min = 0.05, lambda = NULL, 
                        preprocess = c("none"), 
                        centering = FALSE,
                        max.iter = 100, tolerance = 1e-6, 
                        intercept = FALSE){
  # ------------------------------------------- #
  # Needs to be deleted:
  # dyn.load("D:/High-dimensional time series/Exact coordinate descent/rome/src/rome.dll")
  
  # ------------------------------------------- #
  # Input check:
  method <- match.arg(method)
  preprocess <- match.arg(preprocess)
  # screen <- match.arg(screen)
  if (missing(lambda) && nlambda < 2) {
    stop("nlambda (number of lambdas) should be at least 2")
  }
  call <- match.call()
  
  # ------------------------------------------- #
  # Include a column for an intercept & weights & threshold
  intcpt <- 0
  n <- dim(x)[1]
  if (intercept == TRUE){
    XX <- cbind(rep(1,n), x)
    intcpt <- 1
  }else{
    XX <- x
  }
  p <- dim(XX)[2]
  
  if (is.null(weights)){
    ww <- rep(1,n)
  }else{
    ww <- weights
  }
  
  if (is.null(delta)){
    del <- 1.5
  }else{
    del <- delta
  }
  
  # ------------------------------------------- #
  # centering
  shift <- 0
  if (centering == TRUE){
    if (method == "huber") {
      if(gamma > sd(y)) {
        shift <- mean(y)
      }else {
        shift <- median(y)
      }
    } 
  }
  yy <- y - shift
  
  # ------------------------------------------- #
  # Flag for user-supplied lambdas
  user <- 0
  if (is.null(lambda)) {
    lambda <- double(nlambda)
  } else {
    nlambda <- length(lambda)
    user <- 1
  }
  
  # ------------------------------------------- #
  # Flags for preprocessing and screening
  ppflag <- switch(preprocess, standardize = 1L, rescale = 2L, none = 0L)
  # scrflag <- switch(screen, SR = 1L, none = 0L)
  
  # ------------------------------------------- #
  # Fitting
  if (method == "huber") {
    fit <- .C(ecd_huber_, 
              double(p*nlambda), 
              integer(nlambda), 
              as.double(lambda), 
              as.double(XX), 
              as.double(yy), 
              as.double(ww), 
              as.double(del), 
              as.double(tolerance), 
              as.double(lambda.min), 
              as.integer(nlambda),  
              as.integer(n), 
              as.integer(p), 
              as.integer(ppflag), 
              as.integer(max.iter),
              as.integer(user) )
  } 
  
  beta <- matrix(fit[[1]], nrow=p)
  iter <- fit[[2]]
  lambda <- fit[[3]]
  
  # ------------------------------------------- #
  # Get back intercept, if exists
  if (intercept == TRUE){
    beta[1,] <- beta[1,] + shift
  }
  
  # ------------------------------------------- #
  # Names
  vnames <- colnames(x)
  if (intercept == TRUE){
    if (is.null(vnames)) {
      vnames <- paste0("V",seq(p-1))
    }
    vnames <- c("(Intercept)", vnames)
  }else{
    if (is.null(vnames)) {
      vnames <- paste0("V",seq(p))
    }
  }
  dimnames(beta) <- list(vnames, paste0("L", 1:length(lambda)))
  
  # ------------------------------------------- #
  # Output
  structure(list(call = call,
                 beta = beta,
                 iter = iter,
                 lambda = lambda,
                 delta = delta,
                 method = method),
            class = "rome")
  
}
