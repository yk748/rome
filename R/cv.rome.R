cv_rome <- function(y,x,weights=NULL,delta=NULL,n_folds=10){
  
  fold_id <- ceiling(sample(1:n)/n*n_folds)
  
  fit <- rome(x, y, weights=weigths, delta=delta)
  
  n <- length(y)
  E <- matrix(NA, nrow = n, ncol = length(fit$lambda))
  for (i in 1:n_folds) {
    cat("CV fold #",i,"begins\n")
    fit_i_fold <- cv_ftn(i,y,x,fold_id,delta=fit$delta)
    E[fold_id == i, 1:fit_i_fold$nl] <- fit_i_fold$pe
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


cv_ftn <- function(i,y,x,fold_id,delta) {
  x_tr <- x[fold_id != i,,drop = FALSE]
  y_tr <- y[fold_id != i]
  x_ts <- x[fold_id == i,,drop = FALSE]
  y_ts <- y[fold_id == i]
  
  fit_i <- rome(x=x_tr,y=y_tr)
  
  y_hat <- matrix(predict_rome(fit_i,X=x_ts,type="response"),length(y_ts))
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
