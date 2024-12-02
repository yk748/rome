#' Exact Coordinate Descent algorithm for fitting a \deqn{\ell_1}-penalized adaptive Huber regression model
#'
#'
#' Fit solution paths for \deqn{\ell_1}-penalized adaptive Huber regression over a grid of values for the regularization parameter lambda.
#'
#' @param y Response variable.
#' @param x Input matrix, of dimension nobs x nvar;
#' @param lambda Tuning parameter.
#' @param delta Parameter for pivoting.
#' @param w Weight vectors. Default is 1 for each observation. 
#' @param max.iter Maximum number of iterations. Default is 100000.
#' @param tol Tolerance parameter. Default is 1e-7.
#' @return An object with S3 class \code{"huber"}.
#' @author Younghoon Kim, Sumanta Basu, Po-Ling Loh \cr Maintainer: Younghoon Kim
#' \email{yk748@cornell.edu}
#' @examples
#' 
#' set.seed(123)
#' lambda <- 0.05; 
#' delta <- 1.5
#' n <- 100
#' p <- 20
#' x <- matrix(rnorm(n*p,0,1),n,p)
#' e <- rnorm(n,0,1)
#' b <- c(10,-10,rep(0,p-2))
#' y <- x %*% b + e
#' fit <- ecd.huber(y,x,delta,lambda)
#' 
#' @export ecd.huber
#' @useDynLib rome, .registration=TRUE
#' @importFrom Rcpp sourceCpp evalCpp
#' 
ecd.huber <- function(y,x,delta,lambda,w=NULL,max.iter=1e6,tol=1e-7){
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  if (is.null(w)){
    w <- rep(1,n)
  }
  
  beta.cd <- matrix(NA,max.iter,p)
  for (iter in 1:max.iter){
    
    if (iter == 1){
      b.old <- b.new <- rep(0,p)
      beta.cd[1,] <- b.old
    }
    
    for (j in 1:p){
      part_grad <- .C("pg_Huber_reg", as.double(y),as.double(x[,-j]),as.double(w),
                      as.double(b.old[-j]), as.integer(n), as.integer(p), as.double(delta))
      # part_grad <-pg.Huber.reg(b.old[-j],n,y,x[,-j],w,delta)
      if (abs(part_grad) <= lambda){
        b.new[j] <- 0
        next
      }else{
        rj <- (y - x[,-j]%*%b.old[-j])/x[,j]
        
        fit <- .C("cd_huber", as.double(x[,j]),as.double(w),as.double(rj),
                  as.integer(n),as.double(lambda))
        
        b.new[j] <- matrix(fit[[1]],nrow = p)
      }
    }
    
    # print(round(b.new,3))
    beta.cd[iter,] <- b.new
    if (norm(b.new-b.old,"2") < tol){
      # print("converged!, from normal condition check")
      beta.cd <- beta.cd[1:iter,]
      return(list(beta_hat = beta_optimize, iter = iter))
    }else if (iter >= 3){
      if( norm(beta.cd[iter,] - beta.cd[(iter-2),],"2") < tol){
        # print("converged!, from oscillating condition check")
        beta.cd <- beta.cd[1:iter,]
        return(list(beta_hat = beta.cd, iter = iter))
      }
    }else{
      b.old <- b.new
    }
  }
  
  structure(list(beta=beta.cd,iter=iter),class="huber")
}