predict_rome <- function(object, X, lambda, 
                         type=c("response","coefficients","nvars"), 
                         exact = FALSE, ...) {
  type=match.arg(type)
  if (missing(X) && type == "response") {
    stop("Need to supply 'X'")
  }
  beta <- coef_ecd(object,exact = exact)
  
  if (type == "nvars") {
    if (is.matrix(beta)) {
      return(apply(beta!=0, 1, sum))
    }
    else {
      return(sum(beta!=0))
    }
  }
  if (type == "response") {
    return(X %*% t(beta)) 
  }
}