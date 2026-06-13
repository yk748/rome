#' Cross-validation for a penalized robust regression
#'
#' Does k-fold cross-validation for rome, produces a plot, and returns a
#' value for \code{lambda}.
#'
#' @param y response variable.
#' @param x input matrix, of dimension nobs x nvars.
#' @param FUN fitting function to use. Default is \code{"rome.adaptive"}.
#' @param ncores number of cores for parallel CV. Default is 1.
#' @param nfolds number of folds. Default is 10.
#' @param foldid an optional vector of fold assignments (values 1 to nfolds).
#' @param alpha elastic-net mixing parameter in `[0,1]`. \code{alpha=1} is Lasso
#'   (default), \code{alpha=0} is Ridge. Passed explicitly to \code{FUN} and
#'   to each fold fit so that CV is always performed at the correct penalty.
#' @param preprocess preprocessing flag passed to \code{FUN}. One of
#'   \code{"none"} (default), \code{"standardize"}, or \code{"rescale"}.
#' @param type.measure loss for CV: \code{"deviance"} (Huber), \code{"mse"},
#'   or \code{"mae"}.
#' @param track Logical flag to print progress messages on cross-validation. Default is \code{FALSE}.
#' @param ... additional arguments passed to \code{FUN}.
#' @author Younghoon Kim \cr Maintainer: Younghoon Kim
#' \email{yk748@cornell.edu}
#' @references Kim, Y. Loh, PL. S. Basu (2026)
#' \emph{Exact Coordinate Descent for High-Dimensional Regularized Huber Regression, Preprint}
#' @keywords models regression
#' @export cv.rome
cv.rome <- function(x, y,
                    FUN = c("rome.adaptive"),
                    ncores = 1,
                    nfolds = 10,
                    foldid = NULL,
                    alpha = 1,
                    preprocess   = c("none", "standardize", "rescale"),
                    type.measure = c("deviance", "mse", "mae"),
                    track = FALSE,
                    ...) {
  
  FUN <- get(match.arg(FUN))
  type.measure <- match.arg(type.measure)
  preprocess   <- match.arg(preprocess)
  n <- length(y)
  
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1){
    stop("alpha must be a single number in [0, 1]")
  }
  
  # ------------------------------------------- #
  # Fold assignment
  if (is.null(foldid)) {
    foldid <- ceiling(sample(seq_len(n)) / n * nfolds)
  } else {
    nfolds <- max(foldid)
  }
  
  # ------------------------------------------- #
  # Full-data fit to obtain lambda sequence and delta
  fit <- FUN(y = y, x = x, alpha = alpha, preprocess = preprocess, ...)
  cv.args <- list(lambda = fit$lambda, delta = fit$delta, alpha = alpha, preprocess = preprocess)
  measure.args <- list(delta = fit$delta, type.measure = type.measure)
  
  # ------------------------------------------- #
  # CV — parallel or sequential
  parallel <- FALSE
  if (ncores > 1) {
    max.cores <- parallel::detectCores()
    if (ncores > max.cores) {
      if (track){
        message("Requested ", ncores, " cores but only ", max.cores,
                " available; using ", max.cores, ".") 
      }
      ncores <- max.cores
    }
    cluster <- parallel::makeCluster(ncores)
    parallel <- TRUE
    message("Starting parallel cross-validation...")
    parallel::clusterExport(cluster,
                            c("foldid", "y", "x", "cv.args", "measure.args"),
                            envir = environment())
    parallel::clusterCall(cluster, function() require(rome, quietly = TRUE))
    fold.results <- parallel::parLapply(
      cl = cluster,
      X = seq_len(nfolds),
      fun = cv.ftn,
      y = y, x = x, foldid = foldid,
      cv.args = cv.args, measure.args = measure.args, FUN = FUN
    )
    parallel::stopCluster(cluster)
  }
  
  # ------------------------------------------- #
  # Collect fold errors
  E <- matrix(NA, nrow = n, ncol = length(cv.args$lambda))
  for (i in seq_len(nfolds)) {
    if (parallel) {
      fit_i_fold <- fold.results[[i]]
    } else {
      if (track){
        message("CV fold #", i, " begins") 
      }
      fit_i_fold <- cv.ftn(i, x, y, foldid, cv.args, measure.args, FUN)
    }
    E[foldid == i, seq_len(fit_i_fold$nl)] <- fit_i_fold$pe
  }
  
  # ------------------------------------------- #
  # Summary
  idx    <- which(apply(is.finite(E), 2, all))
  E      <- E[, idx, drop = FALSE]
  lambda <- fit$lambda[idx]
  cve    <- colMeans(E)
  cvse   <- apply(E, 2, sd) / sqrt(n)
  
  idx_min <- which.min(cve)
  idx_1se <- min(which(cve <= cve[idx_min] + cvse[idx_min]))
  
  structure(
    list(cve          = cve,
         cvse         = cvse,
         type.measure = type.measure,
         lambda       = lambda,
         alpha        = alpha,
         fit          = fit,
         lambda_min   = lambda[idx_min],
         lambda_1se   = lambda[idx_1se]),
    class = "cv.rome"
  )
}


cv.ftn <- function(i, x, y, foldid, cv.args, measure.args, FUN) {
  
  cv.args$x <- x[foldid != i, , drop = FALSE]
  cv.args$y <- y[foldid != i]
  x_ts <- x[foldid == i, , drop = FALSE]
  y_ts <- y[foldid == i]
  
  fit_i <- do.call(FUN, cv.args)
  
  y_hat <- matrix(predict.rome(fit_i, X = x_ts, type = "response"),
                  nrow = length(y_ts))
  
  list(pe = measure.rome(y_ts, y_hat, measure.args),
       nl = length(fit_i$lambda))
}



measure.rome <- function(y, y_hat, measure.args) {
  r <- y - y_hat
  switch(measure.args$type.measure,
         deviance = huber.loss(r, measure.args$delta),
         mse = r^2,
         mae = abs(r)
  )
}


huber.loss <- function(v, thres) {
  ifelse(abs(v) <= thres, v^2 / 2, thres * abs(v) - thres^2 / 2)
}