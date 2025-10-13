#################################################################################
# Second Numerical Study - KKT conditions
#################################################################################
# Main loops:
get_sim <- function(n,p,type,iter_max){
  
  library("doParallel")
  library("dqrng")
  library("doRNG")
  
  cl <- makeCluster(50)
  registerDoParallel(cl)
  
  list_dummy <- list()
  set.seed(123)
  run_sim <- foreach(iter = 1:iter_max,
                     .errorhandling = 'stop',
                     .packages = c("Matrix","mvtnorm","rome")) %dopar% {
                       
                       source("benchmark_huber.R")
                       
                       if (type == 5){ # AR t4
                         rho <- 0.4
                         
                         sigma = matrix(0,p,p)
                         for (j in 1:p){
                           for (k in 1:p){
                             sigma[j,k] <-  rho^(abs(j-k))
                           }
                         }
                         x <- rmvt(n,sigma=sigma,df=4)
                         e <- rnorm(n,0,1)
                         b <- mapply(j=1:p,function(j)(-1)^j*exp(-(j-1)/10))
                         y <- x %*% b + e
                         
                       }else if (type == 6){ # Block AR (N/N)
                         rho1 <- 0.2
                         rho2 <- 0.8
                         
                         sigma1 <- matrix(NA,p/2,p/2)
                         for (j in 1:p/2){
                           for (k in 1:p/2){
                             sigma1[j,k] <-  rho1^(abs(j-k))
                           }
                         }
                         sigma2 <- matrix(NA,p/2,p/2)
                         for (j in 1:p/2){
                           for (k in 1:p/2){
                             sigma2[j,k] <-  rho2^(abs(j-k))
                           }
                         }
                         x <- cbind(rmvnorm(n,mean=rep(0,p/2),sigma=sigma1),
                                    rmvnorm(n,mean=rep(0,p/2),sigma=sigma2))
                         b <- mapply(j=1:p,function(j)(-1)^j*exp(-(j-1)/10))
                         e <- rnorm(n,0,1)
                         y <- x %*% b + e
                         
                       }
                       
                       gamma <- 0.5
                       lambda_max <- max(mapply(j=1:p,function(j)abs(sum(mapply(i=1:n,function(i)
                         huber_grad(y[i],thresh=gamma) * x[i,j]))/n)))/gamma
                       lambdas <- mapply(j=0:99,function(j)lambda_max*exp((log(0.05)/(100 - 1)))^j)
                       
                       
                       # ---------------------------------------------------------- #
                       # rome with KKT check:
                       start_t_KKT <- proc.time()
                       rome_KKT <- rome.adaptive(y=y,x=x,screen="none",
                                                 lambda=lambdas,delta=gamma,
                                                 adapt=FALSE)
                       end_t_KKT <- proc.time()
                       end_t_KKT - start_t_KKT
                       
                       # ---------------------------------------------------------- #
                       # rome without KKT check:
                       start_t_noKKT <- proc.time()
                       rome_noKKT <- rome.adaptive(y=y,x=x,screen="none",KKT=FALSE,
                                                 lambda=lambdas,delta=gamma,
                                                 adapt=FALSE)
                       end_t_noKKT <- proc.time()
                       end_t_noKKT - start_t_noKKT
                       
                       # -------------------------------------- #
                       # Output
                       output <- list()
                       output$lambdas <- lambdas
                       
                       output$time <- c(c(end_t_KKT - start_t_KKT)[3],
                                        c(end_t_noKKT - start_t_noKKT)[3])
                       
                       beta <- list()
                       beta$true <- b
                       beta$KKT <- rome_KKT$beta
                       beta$noKKT <- rome_noKKT$beta
                       output$beta <- beta
                       
                       rmse <- list()
                       rmse$KKT <- mapply(i=1:100,function(i)
                         norm(rome_KKT$beta[,i] - b,"2")/norm(b,"2"))
                       rmse$noKKT <- mapply(i=1:100,function(i)
                         norm(rome_noKKT$beta[,i] - b,"2")/norm(b,"2"))
                       output$RMSE <- rmse
                       
                       list_dummy[[iter]] <- output
                     }  
  stopCluster(cl)
  save(run_sim,file=paste0("KKT_comparison_type",type,"_p",p,"_n",n,".RData"))
}

#####################################################################################
# Run simulations:
iter_max <- 100
p_vec <- c(50,100,200,500,1000)
n_vec <- c(100,500)
type_vec <- c(5,6)

for (i in 1:5){
  p <- p_vec[i]
  
  for (j in 1:2){
    n <- n_vec[j]
    
    for (k in 1:2){
      type <- type_vec[k]
      
      cat("The current iteration is p:",p,"n:",n,"type is",type,"\n")
      get_sim (n,p,type,iter_max)
      cat("Done.\n")
    }
    
  }
}

#####################################################################################
# Display results:
iter_max <- 100
p_vec <- c(50,100,200,500,1000)
n_vec <- c(100,500)
type_vec <- c(5,6)

dat_time <- data.frame()
for (i in 1:5){
  p <- p_vec[i]
  
  for (j in 1:2){
    n <- n_vec[j]
    
    for (k in 1:2){
      type <- type_vec[k]
      
      cat("The current iteration is p:",p,"n:",n,"type is",type,"\n")
      load(paste0("KKT_comparison_type",type,"_p",p,"_n",n,".RData"))

      time <- round(rowMeans(mapply(iter=1:iter_max,function(iter)run_sim[[iter]]$time)),3)
      names(time) <- NULL
      dat_time <- rbind(dat_time,
                        data.frame(KKT = time[1],
                                   NoKKT = time[2],
                                   P = p,
                                   n = n,
                                   Type = ifelse(type==5,"AR","Block AR")))
      
    }
    
  }
}

dat_time[order(dat_time$Type),]

