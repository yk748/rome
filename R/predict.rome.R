#' make predictions from a "rome" object.
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
#' @seealso \code{rome} and \code{cv.rome}.
#' @keywords models regression
#' @method predict rome
#' @export
#' @export predict.rome
predict.rome <- function(object, newx, s=NULL, 
                         type=c("response","coefficients","nonzero"), 
                         exact = FALSE, ...) {
  type=match.arg(type)
  if( missing(newx) ){
    if(!match(type,c("coefficients","nonzero"),FALSE)){
      stop( "You need to supply a value for 'newx'" )
    }
  }
  
  if( exact&&(!is.null(s)) ){
    lambda <- object$lambda
    which <- match(s,lambda,FALSE)
    if(!all(which>0)){
      lambda <- unique(rev(sort(c(s,lambda))))
      check_dots(object,...)
      object <- update(object,lambda=lambda,...)
    }
  }
  
  beta <- coef.rome(object,exact = exact)
  nbeta <- object$beta
  
  if(!is.null(s)){
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL,NULL)
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda,s)
    
    nbeta <- nbeta[,lamlist$left,drop=FALSE]%*%Diagonal(x=lamlist$frac) +nbeta[,lamlist$right,drop=FALSE]%*%Diagonal(x=1-lamlist$frac)
    namess <- names(s)
    if(is.null(namess))namess <- paste0("s",seq(along=s))
    dimnames(nbeta) <- list(vnames,namess)
  }
  
  if(type=="coefficients"){
    return(nbeta)
  }
  if(type == "nonzero") {
    if (is.matrix(beta)) {
      return(apply(beta!=0, 2, sum))
    }
    else {
      return(sum(beta!=0))
    }
  }
  
  ###Check on newx
  if(inherits(newx, "sparseMatrix")){
    newx <- as(newx,"dMatrix")
  }
  dx <- dim(newx); p <- object$dim[2]
  if(is.null(dx)){
    newx <- matrix(newx,1,byrow=TRUE)
  }
  if(ncol(newx) != p){
    stop(paste0("The number of variables in newx must be ",p))
  }
  nfit <- as.matrix(newx %*% beta)
  nfit
}



