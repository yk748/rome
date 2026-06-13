######################################################################
# This file is used to generate Figures 5 and Table 2 in Section 4.5, 
# as well as Figure A.2 in Section A.1.

# It contains a real-data application. Figure 5 displays realizations 
# of selected variables alongside the sample covariance matrix. 
# Table 2 presents the prediction errors using various thresholds (\delta) 
# and compares them with the results of Lasso estimation using glmnet. 
# Figure A.2 shows the cross-validation results across three different 
# metrics.
######################################################################
rm(list=ls())

library("rome")
library("superheat")
library("mvtnorm")
library("Matrix")
library("glmnet") 

########################################################################
# Data preparation:
load("./glass.RData")
x <- x[,15:500]
y <- y$PbO
n <- dim(x)[1]
p <- dim(x)[2]

# ---------------------------------------------------------- #
# Left of Figure 5:
set.seed(123)
pdf("sample_path.pdf", width = 8, height = 6)
matplot(x[,sample(1:p,30,replace=FALSE)],type='l',ylab="Values",xlab="Vessels")
dev.off()

# ---------------------------------------------------------- #
# Right of Figure 5:
scaled_x <- scale(x)
scaled_y <- scale(y)

gram_matrix <- t(as.matrix(scaled_x)) %*% as.matrix(scaled_x) / n
var_ticks <- 1:p
z_limits <- range(gram_matrix)
color_palette <- hcl.colors(256, "Viridis")

pdf("gram_mat.pdf", width = 8, height = 6)
layout(matrix(c(1,2), nrow = 1), widths = c(8.5, 1.5))
par(mar = c(5,5,4,1)) 
image(x = var_ticks, 
      y = var_ticks, 
      z = t(gram_matrix[p:1, ]),
      col = color_palette,
      zlim = z_limits,
      axes = FALSE, 
      xlab = "Variables", 
      ylab = "Variables",
      main = "Covariance Matrix")
axis_interval <- seq(1,p,by=20)
axis(1, at = axis_interval) 
axis(2, at = axis_interval) 
box()   

par(mar = c(5, 1, 4, 4)) 
scale_sequence <- seq(z_limits[1], z_limits[2], length.out = 256)
image(y = scale_sequence, 
      z = matrix(scale_sequence, nrow = 1), 
      col = color_palette, 
      zlim = z_limits,
      xaxt = "n", yaxt = "n", 
      xlab = "", ylab = "")
axis(4, las = 1) 
box()
par(xpd = TRUE) 
text(x = 1.0, y = z_limits[2],
     labels = "Scale", 
     adj = c(0.5, -1.2),
     cex = 0.8, 
     font = 2)
par(xpd = FALSE)
layout(1)
dev.off()


########################################################################
# Cross-validation results:
set.seed(123)
alpha <- 0.3

# ---------------------------------------------------------- #
# Left of Figure A.2:
cv_rome_05_dev <- cv.rome(x=scaled_x,y=scaled_y,
                          alpha=alpha,delta=0.5,FUN="rome.adaptive")
pdf("cv_deviance.pdf", width = 8, height = 6)
plot(cv_rome_05_dev)
dev.off()

# ---------------------------------------------------------- #
# Middle of Figure A.2:
cv_rome_05_mse <- cv.rome(x=scaled_x,y=scaled_y,
                          alpha=alpha,delta=0.5,type.measure='mse',FUN="rome.adaptive")
pdf("cv_mse.pdf", width = 8, height = 6)
plot(cv_rome_05_mse)
dev.off()

# ---------------------------------------------------------- #
# Right of Figure A.2:
cv_rome_05_mae <- cv.rome(x=scaled_x,y=scaled_y,
                          alpha=alpha,delta=0.5,type.measure='mae',FUN="rome.adaptive")
pdf("cv_mae.pdf", width = 8, height = 6)
plot(cv_rome_05_mae)
dev.off()


########################################################################
# Full data:
set.seed(123)

cv_rome_05 <- cv.rome(x=scaled_x,y=scaled_y,alpha=alpha,delta=0.5,FUN="rome.adaptive")
cv_rome_10 <- cv.rome(x=scaled_x,y=scaled_y,alpha=alpha,delta=1.0,FUN="rome.adaptive")
cv_rome_15 <- cv.rome(x=scaled_x,y=scaled_y,alpha=alpha,delta=1.5,FUN="rome.adaptive")

fit_rome_05 <- rome.adaptive(y=scaled_y,x=scaled_x,screen="adaptive",delta=0.5,
                             lambda=cv_rome_05$lambda_1se,alpha=alpha,
                             adapt=FALSE)
fit_rome_10 <- rome.adaptive(y=scaled_y,x=scaled_x,screen="adaptive",delta=1.0,
                             lambda=cv_rome_10$lambda_1se,alpha=alpha,
                             adapt=FALSE)
fit_rome_15 <- rome.adaptive(y=scaled_y,x=scaled_x,screen="adaptive",delta=1.5,
                             lambda=cv_rome_15$lambda_1se,alpha=alpha,
                             adapt=FALSE)

sum(fit_rome_05$beta!=0)
sum(fit_rome_10$beta!=0)
sum(fit_rome_15$beta!=0)


# Benchmark: glmnet
cv_glmnet_full <- cv.glmnet(x=scaled_x, y=scaled_y, alpha=alpha, 
                            family="gaussian", intercept=FALSE, standardize=FALSE)
fit_glmnet_full <- glmnet(x=scaled_x, y=scaled_y, alpha=alpha, family="gaussian", 
                          lambda=cv_glmnet_full$lambda.1se, intercept=FALSE, standardize=FALSE)
sum(fit_glmnet_full$beta!=0)

########################################################################
# Random partition:
iter_max <- 50
delta_vec <- c(0.5,1.0,1.5)

set.seed(123)
for (j in 1:4){
  if (j !=4){
    delta <- delta_vec[j]
    
    nz <- PE <- rep(NA,iter_max)
    for (iter in 1:iter_max){
      
      idx <- sample(1:n,size=(2/3*n),replace=FALSE)
      x_train <- scaled_x[idx,]
      y_train <- scaled_y[idx,]
      
      x_test <- scaled_x[-idx,]
      y_test <- scaled_y[-idx,]
      
      cv_rome <- cv.rome(x=x_train,y=y_train,alpha=alpha,delta=delta,
                         FUN="rome.adaptive")
      fit_rome <- rome.adaptive(y=scaled_y,x=scaled_x,screen="adaptive",delta=delta,
                                lambda=cv_rome$lambda_1se,alpha=alpha,
                                adapt=FALSE)
      
      nz[iter] <- sum(fit_rome$beta!=0)
      residuals <- y_test - x_test %*% fit_rome$beta
      PE[iter] <- sum(residuals^2)/60
      cat(iter,"th iteration done.\n")
    }
    output <- list()
    output$nz <- nz
    output$PE <- PE
    save(output,file=paste0("Data_appl_delta",delta,".RData"))
    
  }else{
    nz_glmnet <- rep(NA, iter_max)
    PE_glmnet <- rep(NA, iter_max)
    
    for (iter in 1:iter_max) {
      
      idx <- sample(1:n, size = (2/3 * n), replace = FALSE)
      x_train <- scaled_x[idx, ]
      y_train <- scaled_y[idx, ]
      
      x_test <- scaled_x[-idx, ]
      y_test <- scaled_y[-idx, ]
      
      cv_glmnet <- cv.glmnet(x = x_train, y = y_train, alpha = 0.2, 
                             family = "gaussian", intercept=FALSE, standardize=FALSE)
      fit_glmnet <- glmnet(x = x_train, y = y_train, alpha = 0.2, 
                           family = "gaussian", intercept=FALSE, standardize=FALSE, 
                           lambda = cv_glmnet$lambda.1se)
      
      nz_glmnet[iter] <- sum(fit_glmnet$beta != 0)
      residuals_test <- y_test - x_test %*% fit_glmnet$beta
      PE_glmnet[iter] <- sum(residuals_test^2)/60
      cat(iter,"th iteration done.\n")
    }
    
    output_glmnet <- list()
    output_glmnet$nz <- nz_glmnet
    output_glmnet$PE <- PE_glmnet
    save(output_glmnet, file = "Data_appl_glmnet.RData")
  }
}

# ---------------------------------------------------------- #
# Result for Table 2:
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


load("./Data_appl_glmnet.RData")
round(mean(output_glmnet$nz),3)
round(sd(output_glmnet$nz),3)

round(mean(output_glmnet$PE),3)
round(sd(output_glmnet$PE),3)
