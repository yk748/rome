#' Make predictions from a "rome" object.
#'
#' Similar to other predict methods, this functions predicts fitted values,
#' coefficients and number of non-zero coefficients from a fitted \code{"rome"} object.
#' @aliases coef.rome predict.rome 
#' @param object Fitted \code{"rome"} model object.
#' @param X Matrix of new values for \code{x} at which predictions are to be
#' made. Must be a matrix; can be sparse as in \code{Matrix} package. This
#' argument is not used for \code{type=c("coefficients","nonzero")}
#' @param s Value(s) of the penalty parameter \code{lambda} at which
#' predictions are required. Default is the entire sequence used to create the
#' model.
#' @param exact This argument is relevant only when predictions are made at
#' values of \code{s} (lambda) \emph{different} from those used in the fitting
#' of the original model. Not available for \code{"relaxed"} objects. If
#' \code{exact=FALSE} (default), then the predict function uses linear
#' interpolation to make predictions for values of \code{s} (lambda) that do
#' not coincide with those used in the fitting algorithm. While this is often a
#' good approximation, it can sometimes be a bit coarse.  With
#' \code{exact=TRUE}, these different values of \code{s} are merged (and
#' sorted) with \code{object$lambda}, and the model is refit before predictions
#' are made. In this case, it is required to supply the original data \code{x=}
#' and \code{y=} as additional named arguments to \code{predict()} or
#' \code{coef()}. 
#' @seealso \code{rome.adaptive} and \code{cv.rome}.
#' @keywords models regression
#' @method predict rome
#' @export
#' @export predict.rome
predict.rome <- function(object, 
                         X, lambda, 
                         type=c("response","coefficients","nvars"), exact = FALSE, ...) {
  type <- match.arg(type)
  if( missing(X) ){
    if(!match(type,c("coefficients","nonzero"),FALSE)){
      stop( "You need to supply a value for 'new X'" )
    }
  }
  
  beta <- coef.rome(object,lambda=lambda,exact=exact)
  if (is.matrix(beta)) {
    num_coef <- nrow(beta) 
  }else {
    num_coef <- length(beta)
  }
  
  if(type=="coefficients"){
    return(beta)
  } else if (type == "nvars") {
    if (is.matrix(beta)) {
      return(apply(beta!=0, 2, sum))
    }
    else {
      return(sum(beta!=0))
    }
  } else{
    return(X %*% beta)  
  }   
}

coef.rome <- function(object, lambda, exact = FALSE, ...) {
  if (missing(lambda)) {
    beta <- object$beta
    
  } else if (exact) {
    # augment the lambda sequence with the new values, and refit the model
    ls <- object$lambda
    ind <- match(lambda, ls, 0)
    if (any(ind == 0)) {
      ls <- unique(rev(sort(c(lambda,ls))))
      object <- update(object, lambda=ls)
      ind <- match(lambda, ls)
    }
    beta <- object$beta[, ind]
    
  } else {
    
    # use linear interpolation to estimate coefficients for supplied lambda
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

