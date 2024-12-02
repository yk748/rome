# Exact coordinate descent for Huber, univariate:
fun_cd <- function(x,w,r,n,delta,lambda){
  
  thres <- delta/abs(w*x)
  
  dat <- array(NA, c(2*n, 6))
  colnames(dat) <- c("r_i", "i", "r pm delta/|x_j|", "B", "Slope", "f'(c)")
  dat[,1] <- c(r,r) # r_i.
  dat[,2] <- c(rep(1:n, 2)) # indices i of r_i.
  dat[,3] <- c(r-thres, r+thres) # kinks; r_i pm r pm delta/|x_j|.
  dat[,4] <- c(w^2*x^2/n,-w^2*x^2/n) # increment when passing the kinks with respect to c.
  dat <- dat[order(dat[,3]),] # rearrange data in the order of kinks.
  dat[,5] <- cumsum(dat[,4]) # This tell us about the slope.
  
  S <- -delta*sum(abs(w*x))/n
  dat[1,6] <- S + ifelse(dat[1,3] >= 0, lambda, 0) + ifelse(dat[1,3] <= 0, -lambda,0)
  for (i in 1:(nrow(dat)-1)){
    S <- S + dat[i,5]*(dat[i+1,3]-dat[i,3])
    dat[i+1,6] <- S + ifelse(dat[(i+1),3] >= 0, lambda, 0) + ifelse(dat[(i+1),3] <= 0, -lambda,0)
    
    if (i>1){
      if ((dat[(i-1),6] <0) & (dat[i,6] >= 0)){
        x1 <- dat[i-1,3]; 
        x2 <- dat[i,3]
        y0 <- dat[i-1,6]; 
        y1 <- dat[i,6]
        out <- x1 - y0*(x2-x1)/(y1-y0)
      }
    }
  }
  return(out)
}


# partial gradient of penalized Huber regression:
pg.Huber.reg <- function(beta,n,y,x,w,delta){
  sum(mapply(i=1:n,function(i) fun.dL(w[i]*(y[i] - x[i,] %*% beta),delta)))/n
}

# Huber loss:
fun.L <- function(v, delta){
  if (abs(v) <= delta){
    v^2/2
  }else{
    delta*abs(v) - delta^2/2
  }
}

# Huber loss gradient:
fun.dL <- function(v, delta){
  if (abs(v) <= delta){
    v
  }else{
    delta*sign(v)
  }
}