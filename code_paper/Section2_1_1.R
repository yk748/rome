######################################################################
# This file is used to generate the two subplots of Figure 1 in Section 2.1.1.
# It includes two examples of estimating the univariate Huberized median 
# with an elastic penalty.
######################################################################
rm(list=ls())

library(latex2exp)

# --------------------------------------------------------------- #
# Utility functions:
fun_L_Hub <- function(v, delta){
  if (abs(v) <= delta){
    return(v^2/2)
  }else{
    return(delta*abs(v) - delta^2/2)
  }
}

fun_dL_Hub <- function(v, delta){
  if (abs(v) <= delta){
    return(v)
  }else{
    return(delta*sign(v))
  }
}

fval <- function(cv, loss, lambda, alpha){
  return(loss + lambda*(alpha*abs(cv) + (1-alpha)/2*cv^2))
}

fprime <- function(cv, influence, lambda, alpha){
  return(influence + ifelse(cv >= 0,  alpha*lambda, 0)
         + ifelse(cv <= 0, -alpha*lambda, 0) + (1-alpha)*lambda*cv)
}


#################################################################################
# Setup for Figure 1 (a):
set.seed(2026)
n <- 50
x <- sort(unique((round(runif(n,-3,8),3)))) #DGP 1
print(x)
n <- length(x)

## Set up hyperparameters (elastic net)
lambda <- 3; alpha <- 0.5

## Set up hyperparameters (Huber)
delta <- 0.5

# --------------------------------------------------------------- #
# Use built-in optimizers:
fun_pen_Huber <- function(b, x, delta, lambda, alpha){
  sum(mapply(i=1:length(x),function(i) fun_L_Hub(b-x[i],delta))) + lambda*(alpha*abs(b)+(1-alpha)/2*b^2)
}
opt_optim <- optimize(f=fun_pen_Huber,c(min(x)-delta,max(x)+delta),x=x,delta=delta,lambda=lambda,alpha=alpha)
opt_optim$deriv <- fprime(cv=opt_optim$minimum,
                          sum(mapply(i=1:length(x),function(i) 
                            fun_dL_Hub(opt_optim$minimum-x[i],delta))),lambda,alpha)
round(c(opt_optim$minimum,opt_optim$objective,opt_optim$deriv),3)

# --------------------------------------------------------------- #
## Brute force: 
cseq <- sort(unique(c(seq(x[1]-delta, x[n]+delta, by = 0.01))))
dat_orc_Hub <- array(NA, c(length(cseq), 3))
colnames(dat_orc_Hub) <- c("b", "f(b)", "f'(b)")
dat_orc_Hub[,1] <- cseq

for (k in 1:length(cseq)){
  loss_Hub <- inf_Hub <- 0
  for (i in 1:n){
    # Huber:
    loss_Hub <- loss_Hub + fun_L_Hub(cseq[k] - x[i], delta)
    inf_Hub <- inf_Hub + fun_dL_Hub(cseq[k] - x[i], delta)
  }
  dat_orc_Hub[k,2:3] <- c(fval(cseq[k],loss_Hub,lambda,alpha),fprime(cseq[k],inf_Hub,lambda,alpha))
}
j_Hub <- which.min(dat_orc_Hub[,2])
best_c_Hub <- dat_orc_Hub[j_Hub,1]
best_val_Hub <- dat_orc_Hub[j_Hub,2]
best_dval_Hub <- dat_orc_Hub[j_Hub,3]
round(c(best_c_Hub,best_val_Hub,best_dval_Hub),3)

# --------------------------------------------------------------- #
## Exact coordinate descent: 
dat_Hub <- array(NA, c(2*n, 6))
colnames(dat_Hub) <- c("x_i", "i", "x pm delta", "A", "Slope", "f'(b)")
dat_Hub[,1] <- c(x,x) # x_i.
dat_Hub[,2] <- c(rep(1:n, 2)) # indices i of x_i.
dat_Hub[,3] <- c(x-delta, x+delta) # kinks; x_i pm delta.
dat_Hub[,4] <- c(rep(c(1,-1), each = n)) # increment when passing the kinks with respect to c.
dat_Hub <- dat_Hub[order(dat_Hub[,3]),] # rearrange data in the order of kinks.
dat_Hub[,5] <- cumsum(dat_Hub[,4]) # This tell us about the slope.

S_Hub <- -n*delta
dat_Hub[1,6] <- S_Hub + ifelse(dat_Hub[1,3] >= 0, alpha*lambda, 0) + 
  ifelse(dat_Hub[1,3] <= 0, -alpha*lambda,0) + 
  (1-alpha)*lambda*dat_Hub[1,3]
for (i in 1:(nrow(dat_Hub)-1)){
  S_Hub <- S_Hub + dat_Hub[i,5]*(dat_Hub[i+1,3]-dat_Hub[i,3])
  
  dat_Hub[i+1,6] <- S_Hub + ifelse(dat_Hub[(i+1),3] >= 0, alpha*lambda, 0) + 
    ifelse(dat_Hub[(i+1),3] <= 0, -alpha*lambda,0) + 
    (1-alpha)*lambda*dat_Hub[(i+1),3]
  if (i>1){
    if ((dat_Hub[(i-1),6] < 0) & (dat_Hub[i,6] >= 0)){
      idx <- i
      x1 <- dat_Hub[i-1,3]; x2 <- dat_Hub[i,3]
      y0 <- dat_Hub[i-1,6]; y1 <- dat_Hub[i,6]
      cd_c_Hub <- x1 - y0*(x2-x1)/(y1-y0)
      names(cd_c_Hub) <- "b"
        
      cd_val_Hub <- fval(cd_c_Hub, sum(sapply(x, function(xi) fun_L_Hub(cd_c_Hub - xi, delta))), lambda, alpha)
      names(cd_val_Hub) <- "f(b)"
      cd_dval_Hub <- fprime(cd_c_Hub, sum(sapply(x, function(xi) fun_dL_Hub(cd_c_Hub - xi, delta))), lambda, alpha)
      names(cd_dval_Hub) <- "f'(b)"
    }
  }
}
round(c(cd_c_Hub,cd_val_Hub,cd_dval_Hub),3)

# --------------------------------------------------------------- #
# Figure 1 (a):
pdf("rome_univariate1.pdf", width = 8, height = 6)
plot(dat_orc_Hub[,1], dat_orc_Hub[,3], type="l",col="black",
     ylab=TeX("$F'(b)$"),
     xlab=TeX("$b$"),
     main=TeX("$X_i\\sim Unif(-3,8), F'(b) = \\sum_{i=1}^{n}\\rho_{\\delta}'(x_{i}-b) + \\lambda(\\alpha sign(b) + (1-\\alpha)b)$"))
axis(side = 1, at = seq(-3, 8, by = 1))
points(opt_optim$minimum,opt_optim$deriv, pch=18,cex=1.5,col="blue")
points(dat_Hub[,3], dat_Hub[,6], pch=5, cex=0.5, col="red")
points(cd_c_Hub,cd_dval_Hub, pch=18, cex=1.2,col="red")
abline(h=0, col="black")
legend("bottomright", legend=c("Optim","Brute force (grids)","Coordinate descent (kinks)","Coordinate descent (optimal)"),
       col=c("blue","black","red","red"), lty=c(NA,1,NA,NA), 
       pch=c(18,NA,5,18), cex=1, bty = "n")
dev.off()



#################################################################################
# Setup for Figure 1 (b):
set.seed(2024)
n <- 100
x <- sort(unique(round(rt(n+1,df=4),3))) # DGP2
print(x)
n <- length(x)

## Set up hyperparameters (elastic net)
lambda <- 5; alpha <- 1

## Set up hyperparameters (Huber)
delta <- 0.5


# --------------------------------------------------------------- #
# Use built-in optimizers:
fun_pen_Huber <- function(b, x, delta, lambda, alpha){
  sum(mapply(i=1:length(x),function(i) fun_L_Hub(b-x[i],delta))) + lambda*(alpha*abs(b)+(1-alpha)/2*b^2)
}
opt_optim <- optimize(f=fun_pen_Huber,c(min(x)-delta,max(x)+delta),x=x,delta=delta,lambda=lambda,alpha=alpha)
opt_optim$deriv <- fprime(cv=opt_optim$minimum,
                          sum(mapply(i=1:length(x),function(i) 
                            fun_dL_Hub(opt_optim$minimum-x[i],delta))),lambda,alpha)
round(c(opt_optim$minimum,opt_optim$objective,opt_optim$deriv),3)

# --------------------------------------------------------------- #
## Brute force: 
cseq <- sort(unique(c(seq(x[1]-delta, x[n]+delta, by = 0.01))))
dat_orc_Hub <- array(NA, c(length(cseq), 3))
colnames(dat_orc_Hub) <- c("b", "f(b)", "f'(b)")
dat_orc_Hub[,1] <- cseq

for (k in 1:length(cseq)){
  loss_Hub <- inf_Hub <- 0
  for (i in 1:n){
    # Huber:
    loss_Hub <- loss_Hub + fun_L_Hub(cseq[k] - x[i], delta)
    inf_Hub <- inf_Hub + fun_dL_Hub(cseq[k] - x[i], delta)
  }
  dat_orc_Hub[k,2:3] <- c(fval(cseq[k],loss_Hub,lambda,alpha),fprime(cseq[k],inf_Hub,lambda,alpha))
}
j_Hub <- which.min(dat_orc_Hub[,2])
best_c_Hub <- dat_orc_Hub[j_Hub,1]
best_val_Hub <- dat_orc_Hub[j_Hub,2]
best_dval_Hub <- dat_orc_Hub[j_Hub,3]
round(c(best_c_Hub,best_val_Hub,best_dval_Hub),3)

# --------------------------------------------------------------- #
## Exact coordinate descent: 
dat_Hub <- array(NA, c(2*n, 6))
colnames(dat_Hub) <- c("x_i", "i", "x pm delta", "A", "Slope", "f'(b)")
dat_Hub[,1] <- c(x,x) # x_i.
dat_Hub[,2] <- c(rep(1:n, 2)) # indices i of x_i.
dat_Hub[,3] <- c(x-delta, x+delta) # kinks; x_i pm delta.
dat_Hub[,4] <- c(rep(c(1,-1), each = n)) # increment when passing the kinks with respect to c.
dat_Hub <- dat_Hub[order(dat_Hub[,3]),] # rearrange data in the order of kinks.
dat_Hub[,5] <- cumsum(dat_Hub[,4]) # This tell us about the slope.

S_Hub <- -n*delta
dat_Hub[1,6] <- S_Hub + ifelse(dat_Hub[1,3] >= 0, alpha*lambda, 0) + 
  ifelse(dat_Hub[1,3] <= 0, -alpha*lambda,0) + 
  (1-alpha)*lambda*dat_Hub[1,3]
for (i in 1:(nrow(dat_Hub)-1)){
  S_Hub <- S_Hub + dat_Hub[i,5]*(dat_Hub[i+1,3]-dat_Hub[i,3])
  
  dat_Hub[i+1,6] <- S_Hub + ifelse(dat_Hub[(i+1),3] >= 0, alpha*lambda, 0) + 
    ifelse(dat_Hub[(i+1),3] <= 0, -alpha*lambda,0) + 
    (1-alpha)*lambda*dat_Hub[(i+1),3]
  if (i>1){
    if ((dat_Hub[(i-1),6] < 0) & (dat_Hub[i,6] >= 0)){
      idx <- i
      x1 <- dat_Hub[i-1,3]; x2 <- dat_Hub[i,3]
      y0 <- dat_Hub[i-1,6]; y1 <- dat_Hub[i,6]
      cd_c_Hub <- x1 - y0*(x2-x1)/(y1-y0)
      names(cd_c_Hub) <- "b"
      
      cd_val_Hub <- fval(cd_c_Hub, sum(sapply(x, function(xi) fun_L_Hub(cd_c_Hub - xi, delta))), lambda, alpha)
      names(cd_val_Hub) <- "f(b)"
      cd_dval_Hub <- fprime(cd_c_Hub, sum(sapply(x, function(xi) fun_dL_Hub(cd_c_Hub - xi, delta))), lambda, alpha)
      names(cd_dval_Hub) <- "f'(b)"
    }
  }
}
round(c(cd_c_Hub,cd_val_Hub,cd_dval_Hub),3)


# --------------------------------------------------------------- #
# Figure 1 (b):
pdf("rome_univariate2.pdf", width = 8, height = 6)
plot(dat_orc_Hub[,1], dat_orc_Hub[,3], type="l",col="black",
     ylab=TeX("$F'(b)$"),
     xlab=TeX("$b$"),
     main=TeX("$X_i\\sim t(4), F'(b) = \\sum_{i=1}^{n}\\rho_{\\delta}'(x_{i}-b) + \\lambda sign(b)$"))
axis(side = 1, at = seq(-10, 10, by = 1))
points(opt_optim$minimum,opt_optim$deriv, pch=18,cex=1.5,col="blue")
points(dat_Hub[,3], dat_Hub[,6], pch=5, cex=0.5, col="red")
points(cd_c_Hub,cd_dval_Hub, pch=18, cex=1.2,col="red")
abline(h=0, col="black")
legend("bottomright", legend=c("Optim","Brute force (grids)","Coordinate descent (kinks)","Coordinate descent (optimal)"),
       col=c("blue","black","red","red"), lty=c(NA,1,NA,NA), 
       pch=c(18,NA,5,18), cex=1, bty = "n")
dev.off()


