#' Fit a penalized adaptive robust regression
#'
#' Fit a penalized weighted robust regression (Huber loss) via exact coordinate
#' descent algorithm. The weights are either chosen by users or computed internally.
#'
#' @param y response variable.
#' @param x input matrix, of dimension nobs x nvars.
#' @param adapt a logical flag whether the adaptive (weighting) scheme is used. Default is \code{FALSE}.
#' @param method Choice of loss functions. Currently, \code{huber} is only provided.
#' @param weights observation weights. Default is 1 for each observation.
#' @param alpha elastic-net mixing parameter, between 0 and 1. \code{alpha = 1}
#'   is the Lasso (default), \code{alpha = 0} is Ridge. Values in between give
#'   the elastic-net.
#' @param delta threshold for Huber loss function. The default is 0.5.
#' @param nlambda Length of a sequence of \code{lambda}. The default is 100.
#' @param lambda.min A ratio for minimum \code{lambda} from its maximum value. The default is 0.05.
#' @param lambda a sequence of \code{lambda}. It can be provided by users.
#' @param preprocess preprocessing flag. One of \code{"none"} (default),
#'   \code{"standardize"} (center and scale by sd), or \code{"rescale"}
#'   (shift by min, scale by range).
#' @param centering a logical flag whether the data should be centered. Default is \code{FALSE}.
#' @param max.iter maximum iterations per coordinate update. Default is 1e5.
#' @param eps tolerance for convergence. Default is 1e-4.
#' @param intercept a logical flag whether the intercept is added. Default is \code{FALSE}.
#' @param screen screening rule: \code{adaptive}, \code{strong}, or \code{none}.
#' @param KKT a logical flag whether KKT checking is applied. Default is \code{TRUE}.
#' @param dfmax Upper bound for the number of nonzero coefficients. Default is \code{ncol(x)+1}.
#' @param trace Logical flag to print progress messages. Default is \code{FALSE}.
#' @author Younghoon Kim \cr Maintainer: Younghoon Kim
#' \email{yk748@cornell.edu}
#' @references Kim, Y. Loh, PL. S. Basu (2026)
#' \emph{Exact Coordinate Descent for High-Dimensional Regularized Huber Regression, Preprint}
#' @keywords models regression
#' @import utils
#' @import stats
#' @import mvtnorm
#' @useDynLib rome ecd_huber_adaptive_active_
#' @useDynLib rome ecd_huber_benchmarktest_
#' @useDynLib rome ecd_huber_kkttest_
#' @export rome.adaptive

rome.adaptive <- function(y,x,
                          adapt = FALSE,
                          method = c("huber"),
                          weights = NULL,
                          alpha = 1,
                          delta = NULL,
                          nlambda = 100, lambda.min = 0.05, lambda = NULL,
                          preprocess = c("none", "standardize", "rescale"),
                          centering = FALSE,
                          max.iter = 1e5,
                          eps = 1e-4,
                          intercept = FALSE,
                          screen = c("adaptive", "strong", "none"),
                          KKT = TRUE,
                          dfmax = ncol(x) + 1,
                          trace = FALSE) {
  
  # ------------------------------------------- #
  # Input checks
  method     <- match.arg(method)
  preprocess <- match.arg(preprocess)
  screen     <- match.arg(screen)
  
  if (missing(lambda) && nlambda < 2){
    stop("nlambda (number of lambdas) should be at least 2")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1){
    stop("alpha must be a single number in [0, 1]")
  }
  call <- match.call()
  
  # ------------------------------------------- #
  # Intercept column
  n <- nrow(x)
  intcpt <- 0L
  if (intercept) {
    XX <- cbind(rep(1, n), x)
    intcpt <- 1L
  } else {
    XX <- x
  }
  p <- ncol(XX)
  
  # ------------------------------------------- #
  # Weights
  ww <- if (!is.null(weights)) weights else rep(1, n)
  
  # ------------------------------------------- #
  # Huber threshold
  del <- if (is.null(delta)) 0.5 else delta
  
  # ------------------------------------------- #
  # Centering
  shift <- 0
  if (centering) {
    shift <- if (del > sd(y)) mean(y) else median(y)
  }
  yy <- y - shift
  
  # ------------------------------------------- #
  # Lambda sequence
  user <- 0L
  if (is.null(lambda)) {
    lambda <- double(nlambda)
  } else {
    nlambda <- length(lambda)
    user    <- 1L
  }
  
  # ------------------------------------------- #
  # Preprocessing / screening flags
  ppflag  <- switch(preprocess, standardize = 1L, rescale = 2L, none = 0L)
  scrflag <- switch(screen, adaptive = 1L, strong = 2L, none = 0L)
  
  # ------------------------------------------- #
  # Fitting — alpha is passed as the 2nd new argument in every .C() call.
  if (method == "huber") {
    
    if (adapt && scrflag != 0L && KKT) {
      # ---- adaptive weights + screening + KKT ----
      fit <- .C("ecd_huber_adaptive_active_",
                double(p * nlambda),        # [1]  beta path
                integer(nlambda),           # [2]  iter per lambda
                as.double(lambda),          # [3]  lambda sequence (in/out)
                integer(1),                 # [4]  saturated flag
                integer(1),                 # [5]  nv
                as.double(XX),              # [6]  design matrix
                as.double(yy),              # [7]  response
                as.double(ww),              # [8]  weights
                as.double(del),             # [9]  delta
                as.double(alpha),           # [10] alpha  (NEW)
                as.double(eps),             # [11] eps
                as.double(lambda.min),      # [12] lambda.min
                as.integer(nlambda),        # [13] nlambda
                as.integer(n),              # [14] n
                as.integer(p),              # [15] p
                as.integer(ppflag),         # [16] ppflag
                as.integer(dfmax),          # [17] dfmax
                as.integer(max.iter),       # [18] max.iter
                as.integer(user),           # [19] user
                as.integer(scrflag),        # [20] scrflag
                as.integer(trace * 1L))     # [21] trace
      
      beta      <- matrix(fit[[1]], nrow = p, byrow = FALSE)
      iter      <- fit[[2]]
      lambda    <- fit[[3]]
      saturated <- fit[[4]]
      nv        <- fit[[5]]
      
    } else if (!adapt && scrflag != 0L && KKT) {
      # ---- no adaptive weights + screening + KKT (benchmark) ----
      fit <- .C("ecd_huber_benchmarktest_",
                double(p * nlambda),        # [1]  beta path
                integer(nlambda),           # [2]  iter per lambda
                as.double(lambda),          # [3]  lambda sequence (in/out)
                integer(1),                 # [4]  saturated flag
                integer(1),                 # [5]  nv
                as.double(XX),              # [6]  design matrix
                as.double(yy),              # [7]  response
                as.double(del),             # [8]  delta
                as.double(alpha),           # [9]  alpha
                as.double(eps),             # [10] eps
                as.double(lambda.min),      # [11] lambda.min
                as.integer(nlambda),        # [12] nlambda
                as.integer(n),              # [13] n
                as.integer(p),              # [14] p
                as.integer(ppflag),         # [15] ppflag
                as.integer(dfmax),          # [16] dfmax
                as.integer(max.iter),       # [17] max.iter
                as.integer(user),           # [18] user
                as.integer(scrflag),        # [19] scrflag
                as.integer(trace * 1L))     # [20] trace
      
      beta      <- matrix(fit[[1]], nrow = p, byrow = FALSE)
      iter      <- fit[[2]]
      lambda    <- fit[[3]]
      saturated <- fit[[4]]
      nv        <- fit[[5]]
      
    } else if (!adapt && scrflag == 0L) {
      # ---- no screening: KKT on (ecd_huber_kkttest_ with kkt_flag=1)
      #                 or KKT off (ecd_huber_kkttest_ with kkt_flag=0) ----
      fit <- .C("ecd_huber_kkttest_",
                double(p * nlambda),        # [1]  beta path
                integer(nlambda),           # [2]  iter per lambda
                as.double(lambda),          # [3]  lambda sequence (in/out)
                as.double(XX),              # [4]  design matrix
                as.double(yy),              # [5]  response
                as.double(del),             # [6]  delta
                as.double(alpha),           # [7]  alpha
                as.double(eps),             # [8]  eps
                as.double(lambda.min),      # [9] lambda.min
                as.integer(nlambda),        # [10] nlambda
                as.integer(n),              # [11] n
                as.integer(p),              # [12] p
                as.integer(ppflag),         # [13] ppflag
                as.integer(max.iter),       # [14] max.iter
                as.integer(user),           # [15] user
                as.integer(trace * 1L),     # [16] trace
                as.integer(KKT * 1L))       # [17] kkt_flag
      
      beta      <- matrix(fit[[1]], nrow = p, byrow = FALSE)
      iter      <- fit[[2]]
      lambda    <- fit[[3]]
      saturated <- FALSE
      nv        <- NULL
    }
  }
  
  # ------------------------------------------- #
  # Eliminate saturated (non-converged) lambda values
  # (applies to screening branches only)
  if (!is.null(nv)) {
    ind    <- !is.na(iter)
    beta   <- beta[, ind, drop = FALSE]
    iter   <- iter[ind]
    lambda <- lambda[ind]
  }
  
  # ------------------------------------------- #
  # Restore intercept shift
  if (intercept) {
    beta[1, ] <- beta[1, ] + shift
  }
  
  # ------------------------------------------- #
  # Dimension names
  vnames <- colnames(x)
  if (intercept) {
    vnames <- c("(Intercept)", if (is.null(vnames)) paste0("V", seq(p - 1)) else vnames)
  } else {
    if (is.null(vnames)) {
      vnames <- paste0("V", seq(p))
    }
  }
  dimnames(beta) <- list(vnames, paste0("L", seq_len(ncol(beta))))
  
  # ------------------------------------------- #
  # Output
  structure(
    list(call      = call,
         beta      = beta,
         iter      = iter,
         lambda    = lambda,
         delta     = del,
         alpha     = alpha,
         method    = method,
         adapt     = adapt,
         screen    = screen,
         weights   = ww,
         dim       = dim(XX),
         saturated = saturated,
         nv        = nv,
         intercept = as.logical(intcpt)),
    class = "rome"
  )
}