#' fit a penalized adaptive robust regression
#'
#' Fit a penalized weighted  robust regression (Huber loss) via exact coordinate descent algorithm. The weighted are either chosen by users or computed internally.
#' 
#' @param y response variable.
#' @param x input matrix, of dimension nobs x nvars; 
#' @param adapt a logical flag whether the adaptive (weighting) scheme is used. Default is \code{FALSE}.
#' @param weights observation weights. Can be total counts if responses are
#' proportion matrices. Default is 1 for each observation
#' @param delta threshold for Huber loss function. The default is 1.5. 
#' @param nlambda Length of a sequence of \code{lambda}. The default length is 100.
#' @param lambda.min A ratio for minimum \code{lambda} from its maximum value. The default is 0.05.
#' @param lambda a sequence of \code{lambda}. It can be provided by users.
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
#' \emph{Exact Coordinate Descent for High-Dimensional Regularized Robust M-Estimators, ??, Vol. ??(??), ??-??},
#' \doi{??}.\cr
#' @keywords models regression
#' @import utils
#' @import stats
#' @import mvtnorm
#' @useDynLib rome ecd_huber_adaptive_active_
#' @useDynLib rome ecd_huber_active_
#' @useDynLib rome ecd_huber_noscreen_
#' @export rome_adaptive
#' 
rome_adaptive <- function(y, x, adapt = FALSE,
                          method = c("huber"), weights=NULL,
                        delta = NULL, 
                        nlambda=100, lambda.min = 0.05, lambda = NULL, 
                        preprocess = c("none"), 
                        centering = FALSE,
                        max.iter = 1e5, eps = 1e-4, 
                        intercept = FALSE, screen=c("adaptive","strong","none"),
                        dfmax = ncol(x)+1,trace = FALSE){

  # ------------------------------------------- #
  # Input check:
  method <- match.arg(method)
  preprocess <- match.arg(preprocess)
  screen <- match.arg(screen)
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
  
  if (!is.null(weights)){
    ww <- weights
  }else {
    ww <- rep(1,n)
  }
  
  if (is.null(delta)){
    del <- 0.5
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
  scrflag <- switch(screen, adaptive = 1L, strong = 2L, none = 0L)
  
  # ------------------------------------------- #
  # Fitting
  if (method == "huber") {
    if (adapt & scrflag !=0){
      fit <- .C("ecd_huber_adaptive_active_", 
                double(p*nlambda), 
                integer(nlambda), 
                as.double(lambda), 
                integer(1),
                integer(1),
                as.double(XX), 
                as.double(yy), 
                as.double(ww), 
                as.double(del), 
                as.double(eps), 
                as.double(lambda.min), 
                as.integer(nlambda),  
                as.integer(n), 
                as.integer(p), 
                as.integer(ppflag), 
                as.integer(dfmax), 
                as.integer(max.iter),
                as.integer(user),
                as.integer(scrflag),
                as.integer(c(trace)*1))
      
      beta <- matrix(fit[[1]], nrow=p, byrow=FALSE)
      iter <- fit[[2]]
      lambda <- fit[[3]]
      saturated <- fit[[4]]
      nv <- fit[[5]]
      # Eliminate saturated lambda values
      ind <- !is.na(iter)
      beta <- beta[, ind]
      iter <- iter[ind]
      lambda <- lambda[ind]
      
      delta <- fit[[9]]
      
      
    }else if (!adapt & scrflag !=0){
      fit <- .C("ecd_huber_active_", 
                double(p*nlambda), 
                integer(nlambda), 
                as.double(lambda), 
                integer(1),
                integer(1),
                as.double(XX), 
                as.double(yy), 
                as.double(del), 
                as.double(eps), 
                as.double(lambda.min), 
                as.integer(nlambda),  
                as.integer(n), 
                as.integer(p), 
                as.integer(ppflag), 
                as.integer(dfmax), 
                as.integer(max.iter),
                as.integer(user),
                as.integer(scrflag),
                as.integer(c(trace)*1)) 
      
      beta <- matrix(fit[[1]], nrow=p, byrow=FALSE)
      iter <- fit[[2]]
      lambda <- fit[[3]]
      saturated <- fit[[4]]
      nv <- fit[[5]]
      # Eliminate saturated lambda values
      ind <- !is.na(iter)
      beta <- beta[, ind]
      iter <- iter[ind]
      lambda <- lambda[ind]
      
      delta <- fit[[9]]
      
    }else{
      fit <- .C("ecd_huber_noscreen_", 
                double(p*nlambda), 
                integer(nlambda), 
                as.double(lambda), 
                as.double(XX), 
                as.double(yy), 
                as.double(ww),
                as.double(del), 
                as.double(eps), 
                as.double(lambda.min), 
                as.integer(nlambda),  
                as.integer(n), 
                as.integer(p), 
                as.integer(ppflag), 
                as.integer(max.iter),
                as.integer(user),
                as.integer(c(trace)*1)) 
    }
    
    beta <- matrix(fit[[1]], nrow=p, byrow=FALSE)
    iter <- fit[[2]]
    lambda <- fit[[3]]
    delta <- fit[[7]]
    
    saturated <- FALSE
    nv <- NULL
  } 
  
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
                 method = method,
                 adapt = adapt,
                 screen = screen,
                 weights = ww,
                 dim = dim(XX),
                 saturated = saturated,
                 nv = nv,
                 intercept = ifelse(intcpt,TRUE,FALSE)),
            class = "rome")
  
}
