#################################################################################
# First Numerical Study - Comparisons with benchmarks
#################################################################################
get_sim <- function(n,p,type,iter_max){
  
  library("doParallel")
  library("dqrng")
  library("doRNG")
  
  dyn.load("GD_huber.so")
  # dyn.load("GD_huber.dll") # If running on WINDOWS
  
  cl <- makeCluster(50)
  registerDoParallel(cl)

  list_dummy <- list()
  set.seed(123)
  run_sim <- foreach(iter = 1:iter_max,
                     .errorhandling = 'stop',
                     .packages = c("Matrix","mvtnorm","glmnet","hqreg","rome","ILAMM")) %dopar% {
                       
                     huber_grad <- function(v, thresh){
                       if (abs(v) <= thresh){
                         return(v)
                       }else{
                         return(thresh*sign(v))
                       }
                     }
                     
                     # ---------------------------------------------------------- #
                     # Data generation:
                     if (type == 1){ # compound N
                       rho <- 0.8;
                       sigma = matrix(rho,p,p); diag(sigma) <- 1
                       x <- rmvnorm(n,mean=rep(0,p),sigma=sigma)
                       e <- rnorm(n,0,1)
                       b <- c(2,0,1.5,0,0.8,0,0,1,0,1.75,0,0,0.75,0,0,0.3,rep(0,p-16))
                       y <- x %*% b + e
                       
                     }else if (type == 2){ # AR t2
                       rho <- 0.8; v <- 2
                       sigma = matrix(0,p,p)
                       for (j in 1:p){
                         for (k in 1:p){
                           sigma[j,k] <-  rho^(abs(j-k))
                         }
                       }
                       x <- rmvt(n,sigma=sigma,df=v)
                       e <- rnorm(n,0,1)
                       b <- c(2,0,1.5,0,0.8,0,0,1,0,1.75,0,0,0.75,0,0,0.3,rep(0,p-16))
                       y <- x %*% b + e
                       
                     }else if (type == 3){ # Contaminated AR
                       rho <- 0.8;  v <- 1
                       sigma <- matrix(NA,(p-1),(p-1))
                       for (j in 1:(p-1)){
                         for (k in 1:(p-1)){
                           sigma[j,k] <-  rho^(abs(j-k))
                         }
                       }
                       x <- cbind(rmvnorm(n,mean=rep(0,p-1),sigma=sigma),
                                  rt(n,df=v))
                       e <- rnorm(n,0,1)
                       b <- c(2,0,1.5,0,0.8,0,0,1,0,1.75,0,0,0.75,0,0,0.3,rep(0,p-16))
                       y <- x %*% b + e
                       
                     }else if (type == 4){ # Block AR
                       rho1 <- 0.2; rho2 <- 0.8; v <- 1
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
                       x <- cbind(rmvt(n,sigma=sigma1),
                                  rmvnorm(n,mean=rep(0,p/2),sigma=sigma2))
                       b <- c(2,0,1.5,0,0.8,0,0,1,0,1.75,0,0,0.75,0,0,0.3,rep(0,p-16))
                       e <- rnorm(n,0,1)
                       y <- x %*% b + e
                     }
                     
                     gamma <- 0.5
                     lambda_max <- max(mapply(j=1:p,function(j)abs(sum(mapply(i=1:n,function(i)
                       huber_grad(y[i],thresh=gamma) * x[i,j]))/n)))/gamma
                     lambdas <- mapply(j=0:99,function(j)lambda_max*exp((log(0.05)/(100 - 1)))^j)
                     
                     # -------------------------------------- #
                     # rome
                     start_t_rome <- proc.time()
                     fit_rome <- rome.adaptive(y=y,x=x,screen="adaptive",
                                                lambda=lambdas,delta=gamma,
                                                adapt=FALSE)
                     end_t_rome <- proc.time()
      
                     # -------------------------------------- #
                     # hqreg
                     start_t_hqreg <- proc.time()
                     fit_hqreg <- hqreg_raw(X=x,y=y,method="huber",
                                             lambda=lambdas/gamma,gamma=gamma,
                                             intercept=FALSE,screen="ASR")
                     end_t_hqreg <- proc.time()
                     
                     # ---------------------------------------------------------- #
                     # ILAMM:
                     start_t_ILAMM <- proc.time()
                     beta_ILAMM <- matrix(NA,p,length(lambdas))
                     for (l in 1:length(lambdas)){
                       fit <- ncvxHuberReg(X=x, Y=y, lambda=lambdas[l],
                                           tau=gamma, intercept=FALSE, penalty="Lasso")
                       beta_ILAMM[,l] <- fit$beta[-1]
                     }
                     end_t_ILAMM <- proc.time()
                     
                     # ---------------------------------------------------------- #
                     # composite gradient descent:
                     start_t_CGD <- proc.time()
                     beta_CGD <- matrix(NA,p,length(lambdas))
                     for (l in 1:length(lambdas)){
                       beta_CGD[,l] <- GD_huber_l1(X=x,y=y,delta=gamma,lambda=lambdas[l])
                     }
                     end_t_CGD <- proc.time()
                     
                     # -------------------------------------- #
                     # Output
                     output <- list()
                     output$time <- c(c(end_t_rome - start_t_rome)[3],
                                    c(end_t_hqreg - start_t_hqreg)[3],
                                    c(end_t_ILAMM - start_t_ILAMM)[3],
                                    c(end_t_CGD - start_t_CGD)[3])
                     
                     beta <- list()
                     beta$true <- b
                     beta$rome <- fit_rome$beta
                     beta$hqreg <- fit_hqreg$beta
                     beta$ILAMM <- beta_ILAMM
                     beta$CGD <- beta_CGD
                     output$beta <- beta
                     
                     rmse <- list()
                     rmse$rome <- mapply(i=1:100,function(i)
                       norm(fit_rome$beta[,i] - b,"2")/norm(b,"2"))
                     rmse$hqreg <- mapply(i=1:100,function(i)
                       norm(fit_hqreg$beta[,i] - b,"2")/norm(b,"2"))
                     rmse$ILAMM <- mapply(i=1:100,function(i)
                       norm(beta_ILAMM[,i] - b,"2")/norm(b,"2"))
                     rmse$CGD <- mapply(i=1:100,function(i)
                       norm(beta_CGD[,i] - b,"2")/norm(b,"2"))
                     output$RMSE <- rmse
                     list_dummy[[iter]] <- output
                   }  
  stopCluster(cl)
  save(run_sim,file=paste0("path_comparison_type_",type,"_p",p,"_n",n,".RData"))
}

#####################################################################################
# Run simulations:
iter_max <- 100
n_vec <- c(100,500,1000)
p_vec <- c(100,500,1000,5000)
type_vec <- c(1,2,3,4)
type_vec <- c(1,5,4,6)

for (i in 1:4){
  p <- p_vec[i]
  
  for (j in 1:3){
    n <- n_vec[j]
    
    for (k in 1:4){
     type <- type_vec[k] 
     
     cat("The current iteration is p:",p,"n:",n,"type is",type,"\n")
     get_sim (n,p,type,iter_max)
     cat("Done.\n")
    }
  }
}


#####################################################################################
# Draw plots:
library(ggplot2)
library(dplyr)
library(tidyr)

p_vec <- c(100,500,1000,5000)
n_vec <- c(100,500,1000)
type_vec <- c(1,4,5,6)

all_data <- list()
time_data <- list()
for (i in 1:4){
  p <- p_vec[i]
  
  for (j in 1:3){
    n <- n_vec[j]
    
    for (k in 1:4){
      type <- type_vec[k]

      cat("The current iteration is p:",p,"n:",n,"type is",type,"\n")
      
      load(paste0("path_comparison_type_",type,"_p",p,"_n",n,".RData"))
      
      # ----------------------------------------------- #
      df_time <- data.frame(rome=mapply(l=1:100,function(l)run_main[[l]]$time[1]),
                            hqreg=mapply(l=1:100,function(l)run_main[[l]]$time[2]),
                            ILAMM=mapply(l=1:100,function(l)run_ILAMM[[l]]$time[1]),
                            GD=mapply(l=1:100,function(l)run_GD[[l]]$time[1]),
                            n = rep(n,100),
                            p = rep(p,100),
                            Type = rep(type,100))
      time_data[[length(time_data) + 1]] <- df_time

      # ----------------------------------------------- #
      rmse_rome <- log(mapply(l=1:100,function(l)run_main[[l]]$RMSE$rome))
      mean_rome <- rowMeans(rmse_rome)
      upper_rome <- mapply(l=1:100,function(l)quantile(rmse_rome[l,],0.95))
      lower_rome <- mapply(l=1:100,function(l)quantile(rmse_rome[l,],0.05))

      rmse_hqreg <- log(mapply(l=1:100,function(l)run_main[[l]]$RMSE$hqreg))
      mean_hqreg <- rowMeans(rmse_hqreg)
      upper_hqreg <- mapply(l=1:100,function(l)quantile(rmse_hqreg[l,],0.95))
      lower_hqreg <- mapply(l=1:100,function(l)quantile(rmse_hqreg[l,],0.05))
      
      rmse_ILAMM <- log(mapply(l=1:100,function(l)run_ILAMM[[l]]$RMSE$ILAMM))
      mean_ILAMM <- rowMeans(rmse_ILAMM)
      upper_ILAMM <- mapply(l=1:100,function(l)quantile(rmse_ILAMM[l,],0.95))
      lower_ILAMM <- mapply(l=1:100,function(l)quantile(rmse_ILAMM[l,],0.05))
      
      rmse_GD <- log(mapply(l=1:100,function(l)run_GD[[l]]$RMSE$CGD))
      mean_GD <- rowMeans(rmse_GD)
      upper_GD <- mapply(l=1:100,function(l)quantile(rmse_GD[l,],0.95))
      lower_GD <- mapply(l=1:100,function(l)quantile(rmse_GD[l,],0.05))

      df <- data.frame(
        Lambda = rep(1:100,4),
        Mean = c(mean_rome,mean_hqreg,mean_ILAMM,mean_GD),
        Q95 = c(upper_rome,upper_hqreg,upper_ILAMM,upper_GD),
        Q05 = c(lower_rome,lower_hqreg,lower_ILAMM,lower_GD),
        n = rep(n,400),
        p = rep(p,400),
        Type = rep(type,400),
        Method = c(rep("rome",100),rep("hqreg",100),rep("ILAMM",100),rep("GD",100))
      )
      all_data[[length(all_data) + 1]] <- df
    }
  }
}

# ----------------------------------------------- #
# Time comparison:
plot_time <- do.call(rbind, time_data)
df_long <- pivot_longer(plot_time, cols = c(rome, hqreg, ILAMM, GD),
                        names_to = "Method", values_to = "Time")
df_long$Method <- factor(df_long$Method, levels = c("rome","hqreg","ILAMM","GD"))
df_long$pn_label <- factor(
  paste0("p=", df_long$p, ", n=", df_long$n),
  levels = c(
    "p=100, n=100", "p=100, n=500", "p=100, n=1000",
    "p=500, n=100", "p=500, n=500", "p=500, n=1000",
    "p=1000, n=100", "p=1000, n=500", "p=1000, n=1000",
    "p=5000, n=100", "p=5000, n=500", "p=5000, n=1000"))
df_long$Type <- factor(df_long$Type,
                         levels = c(1, 5, 4, 6),
                         labels = c("Compound Normal","AR","Contaminated AR","Block AR"))


df_no_GD <- subset(df_long, Method != "GD")
plot_time1 <- ggplot(df_no_GD, aes(x = Method, y = Time, fill = Method)) +
                    geom_boxplot(outlier.size = 0.5, width = 0.6, alpha = 0.5) +
                    facet_grid(Type ~ pn_label, scales = "free_y") +
                    scale_fill_manual(values = c("rome" = "red", "hqreg" = "blue", "ILAMM" = "green")) +
                    labs( title = "Computation Time for all Lambdas (excluding GD)",
                          y = "Time (CPU seconds)",
                          x = "" ) +
                    theme_minimal() +
                    theme(
                      strip.text = element_text(size = 8),
                      legend.position = "bottom",
                      axis.ticks.x = element_blank(),
                      axis.text.x = element_blank() )

pdf("comparison_time_noGD.pdf", width = 13, height = 11)
plot_time1
dev.off()


df_GD <- subset(df_long, Method == "GD")
plot_time2 <- ggplot(df_GD, aes(x = Method, y = Time, fill = Method)) +
                    geom_boxplot(outlier.size = 0.5, width = 0.6, alpha = 0.5) +
                    facet_grid(Type ~ pn_label, scales = "free_y") +
                    scale_fill_manual(values = c("GD" = "orange")) +
                    labs( title = "Computation Time for GD",
                      y = "Time (CPU seconds)",
                      x = "" ) +
                    theme_minimal() +
                    theme( strip.text = element_text(size = 8),
                      legend.position = "bottom",
                      axis.ticks.x = element_blank(),
                      axis.text.x = element_blank() )

pdf("comparison_time_GD.pdf", width = 13, height = 11)
plot_time2
dev.off()


# ----------------------------------------------- #
# RMSE comparison:
plot_data <- do.call(rbind, all_data)
plot_data$Method <- factor(plot_data$Method, levels = c("rome","hqreg","ILAMM","GD"))
plot_data$pn_label <- factor(
  paste0("p=", plot_data$p, ", n=", plot_data$n),
  levels = c(
    "p=100, n=100", "p=100, n=500", "p=100, n=1000",
    "p=500, n=100", "p=500, n=500", "p=500, n=1000",
    "p=1000, n=100", "p=1000, n=500", "p=1000, n=1000",
    "p=5000, n=100", "p=5000, n=500", "p=5000, n=1000"))
plot_data$Type <- factor(plot_data$Type,
                         levels = c(1, 5, 4, 6),
                         labels = c("Compound Normal","AR","Contaminated AR","Block AR"))


plot_rmse <- ggplot(plot_data, aes(x = Lambda, y = Mean, color = Method, fill = Method)) +
                    geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.2, color = NA, show.legend = FALSE) +
                    geom_line(size = 0.8) +
                    facet_grid(Type ~ pn_label, scales = "free_y") +  
                    labs(title = "RMSE over Lambda Sequences", y = "log(RMSE)", x="Lambda sequence") +
                    scale_color_manual(name="Method",
                                       values=c("rome"="red","hqreg"="blue","ILAMM"="green","GD"="orange")) +
                    scale_fill_manual(name="Method",
                                      values=c("rome"="red","hqreg"="blue","ILAMM"="green","GD"="orange")) +
                    theme_minimal() +
                    theme(
                      strip.text = element_text(size = 8),
                      legend.position = "bottom",
                      axis.ticks.x = element_blank(),
                      axis.text.x = element_blank()
                    )

pdf("comparison_rmse.pdf", width = 13, height = 11)
plot_rmse
dev.off()



