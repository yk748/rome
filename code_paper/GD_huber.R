
########################################################################
# Huber loss function
huber_loss <- function(v, thresh){
  ifelse(abs(v)<=thresh,0.5*v^2,thresh*(abs(v)-0.5*thresh))
}

########################################################################
# Derivative (gradient) of Huber loss
huber_grad <- function(v, thresh){
  ifelse(abs(v)<=thresh,v,thresh*sign(v))
}

########################################################################
# Soft-thresholding operator (L1)
soft_threshold <- function(v, thres) {
  pmax(0, v - thres) - pmax(0, -v - thres)
}

########################################################################
# Composite Gradient Descent for L1-penalized Huber regression
GD_huber_l1 <- function(X, y, delta = 1, lambda = 0.1,
                        step_size = -1, 
                        max_iter = 10000,
                        tol = 1e-4, verbose = FALSE) {
  if (!is.matrix(X)) stop("X must be a numeric matrix")
  if (!is.numeric(y)) stop("y must be numeric")
  
  dyn.load("GD_huber.so")
  .Call("comp_grad_huber_l1", X, y,
        as.numeric(delta), as.numeric(lambda),
        as.numeric(step_size), as.integer(max_iter),
        as.numeric(tol), as.integer(verbose))
}
