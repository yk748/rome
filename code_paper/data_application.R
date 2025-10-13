########################################################################
# Data application
########################################################################
rm(list=ls())

library("rome")
library("superheat")
library("mvtnorm")
library("Matrix")

source("benchmark_huber.R")
########################################################################
load("./glass.RData")
x <- x[,15:500]
y <- y$PbO
n <- dim(x)[1]
p <- dim(x)[2]

set.seed(123)
matplot(x[,sample(1:p,30,replace=FALSE)],type='l',ylab="Values",xlab="Vessels")
superheat(as.matrix(x) %*% t(as.matrix(x)))

scaled_x <- scale(x,center=TRUE,scale=TRUE)
scaled_y <- scale(y,center=TRUE,scale=TRUE)

########################################################################
# Cross-validation results:
cv_rome_05_dev <- cv.rome(x=scaled_x,y=scaled_y,delta=0.5,FUN="rome.adaptive")
pdf("cv_deviance.pdf", width = 13, height = 11)
plot(cv_rome_05_dev)
dev.off()

cv_rome_05_mse <- cv.rome(x=scaled_x,y=scaled_y,delta=0.5,type.measure='mse',FUN="rome.adaptive")
pdf("cv_mse.pdf", width = 13, height = 11)
plot(cv_rome_05_mse)
dev.off()

cv_rome_05_mae <- cv.rome(x=scaled_x,y=scaled_y,delta=0.5,type.measure='mae',FUN="rome.adaptive")
pdf("cv_mae.pdf", width = 13, height = 11)
plot(cv_rome_05_mae)
dev.off()


########################################################################
# Full data:
cv_rome_05 <- cv.rome(x=scaled_x,y=scaled_y,delta=0.5,FUN="rome.adaptive")
cv_rome_10 <- cv.rome(x=scaled_x,y=scaled_y,delta=1.0,FUN="rome.adaptive")
cv_rome_15 <- cv.rome(x=scaled_x,y=scaled_y,delta=1.5,FUN="rome.adaptive")

fit_rome_05 <- rome.adaptive(y=scaled_y,x=scaled_x,screen="adaptive",delta=0.5,
                           lambda=cv_rome_05$lambda_1se,
                           adapt=FALSE)
fit_rome_10 <- rome.adaptive(y=scaled_y,x=scaled_x,screen="adaptive",delta=1.0,
                             lambda=cv_rome_10$lambda_1se,
                             adapt=FALSE)
fit_rome_15 <- rome.adaptive(y=scaled_y,x=scaled_x,screen="adaptive",delta=1.5,
                             lambda=cv_rome_15$lambda_1se,
                             adapt=FALSE)

sum(fit_rome_05$beta!=0)
sum(fit_rome_10$beta!=0)
sum(fit_rome_15$beta!=0)

########################################################################
# Random partition:
iter_max <- 50
delta_vec <- c(0.5,1.0,1.5)

set.seed(123)
for (j in 1:3){
  delta <- delta_vec[j]
  
  nz <- PE <- rep(NA,iter_max)
  for (iter in 1:iter_max){
    
    idx <- sample(1:n,size=(2/3*n),replace=FALSE)
    x_train <- scaled_x[idx,]
    y_train <- scaled_y[idx,]
    
    x_test <- scaled_x[-idx,]
    y_test <- scaled_y[-idx,]
    
    cv_rome <- cv.rome(x=x_train,y=y_train,delta=delta,FUN="rome.adaptive")
    fit_rome <- rome.adaptive(y=scaled_y,x=scaled_x,screen="adaptive",delta=delta,
                              lambda=cv_rome$lambda_1se,
                              adapt=FALSE)
    
    nz[iter] <- sum(fit_rome$beta!=0)
    PE[iter] <- sum(huber_loss(y_test - x_test %*% fit_rome$beta,delta))/60
    cat(iter,"th iteration done.\n")
  }
  output <- list()
  output$nz <- nz
  output$PE <- PE
  save(output,file=paste0("Data_appl_delta",delta,".RData"))
}

load("./Data_appl_delta0.5.RData")
round(mean(output$nz),3)
round(sd(output$nz),3)

round(mean(output$PE),3)
round(sd(output$PE),3)


load("./Data_appl_delta1.RData")
round(mean(output$nz),3)
round(sd(output$nz),3)

round(mean(output$PE),3)
round(sd(output$PE),3)

load("./Data_appl_delta1.5.RData")
round(mean(output$nz),3)
round(sd(output$nz),3)

round(mean(output$PE),3)
round(sd(output$PE),3)
