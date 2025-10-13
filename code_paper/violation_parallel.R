#################################################################################
# Second Numerical Study - screening rule
#################################################################################
# Main loops:
get_sim <- function(n,p,type,rho_type,iter_max){
  
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
                     
                     # ---------------------------------------------------------- #
                     # Data generation:
                     if (type == 5){ # AR t4
                       if (rho_type == 1){
                         rho <- 0.0
                       }else if (rho_type == 2){
                         rho <- 0.8
                       }
                       
                       sigma = matrix(0,p,p)
                       for (j in 1:p){
                         for (k in 1:p){
                           sigma[j,k] <-  rho^(abs(j-k))
                         }
                       }
                       x <- rmvt(n,sigma=sigma,df=4)
                       e <- rnorm(n,0,1)
                       ones <- rep(1,floor(p*0.1))
                       b <- sample(c(ones,rep(0,p-floor(p*0.1))),p,replace=FALSE)
                       y <- x %*% b + e
                       
                     }else if (type == 6){ # Block AR (N/N)
                       if (rho_type == 1){
                         rho1 <- 0.0
                       }else if (rho_type == 2){
                         rho1 <- 0.4
                       }
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
                         
                         ones <- rep(1,floor(p*0.1))
                         b <- sample(c(ones,rep(0,p-floor(p*0.1))),p,replace=FALSE)
                         e <- rnorm(n,0,1)
                         y <- x %*% b + e
                         
                       }
                       
                       gamma <- 0.5
                       lambda_max <- max(mapply(j=1:p,function(j)abs(sum(mapply(i=1:n,function(i)
                         huber_grad(y[i],thresh=gamma) * x[i,j]))/n)))/gamma
                       lambdas <- mapply(j=0:99,function(j)lambda_max*exp((log(0.05)/(100 - 1)))^j)
                       
                       # ---------------------------------------------------------- #
                       # rome with ASR:
                       start_t_ASR <- proc.time()
                       rome_ASR <- rome.adaptive(y=y,x=x,screen="adaptive",
                                                 lambda=lambdas,delta=gamma,
                                                 adapt=FALSE)
                       end_t_ASR <- proc.time()
                       
                       
                       # ---------------------------------------------------------- #
                       # rome with SSR:
                       start_t_SSR <- proc.time()
                       rome_SSR <- rome.adaptive(y=y,x=x,screen="strong",
                                                 lambda=lambdas,delta=gamma,
                                                 adapt=FALSE)
                       end_t_SSR <- proc.time()
                       
                       # ---------------------------------------------------------- #
                       # violation
                       c_path_ASR <- matrix(NA,100,p)
                       c_path_SSR <- matrix(NA,100,p)
                       for (l in 1:100){
                         r_ASR_l <- y - x %*% rome_ASR$beta[,l]
                         c_path_ASR[l,] <- mapply(j=1:p,function(j)-sum(mapply(i=1:100,function(i)
                           x[i,j]*huber_grad(r_ASR_l[i],thresh=gamma)))/n)
                         
                         r_SSR_l <- y - x %*% rome_SSR$beta[,l]
                         c_path_SSR[l,] <- mapply(j=1:p,function(j)-sum(mapply(i=1:100,function(i)
                           x[i,j]*huber_grad(r_SSR_l[i],thresh=gamma)))/n)
                       }
                       
                       # Rules setup:
                       M_l_ASR <- ASR_rule <- rep(NA,100)
                       
                       M_l_ASR[1] <- M_l_hqreg[1] <- 1
                       ASR_rule[1] <- 2*lambda_max - lambdas[2]
                       for (l in 2:100){
                         # ASR
                         M_l_ASR[l] <- max(mapply(j=1:p,function(j)abs(c_path_ASR[(l-1),j] 
                                                                       - c_path_ASR[l,j])))/(lambdas[(l-1)]-lambdas[l])
                         ASR_rule[l] <- lambdas[l] + M_l_ASR[l]*(lambdas[(l-1)]-lambdas[l])
                       }
                       # SSR
                       SSR_rule <- c(lambda_max,lambdas[-100] + diff(lambdas))
                       
                       viol_ASR <- viol_SSR <- matrix(NA,100,p)
                       for (j in 1:p){
                         viol_ASR[,j] <- 1*(abs(c_path_ASR[,j]) >= ASR_rule)
                         viol_SSR[,j] <- 1*(abs(c_path_SSR[,j]) >= SSR_rule)
                       }
                       
                       # -------------------------------------- #
                       # Output
                       output <- list()
                       output$lambdas <- lambdas
                       
                       output$time <- c(c(end_t_ASR - start_t_ASR)[3],
                                        c(end_t_SSR - start_t_SSR)[3])
                       
                       beta <- list()
                       beta$true <- b
                       beta$ASR <- rome_ASR$beta
                       beta$SSR <- rome_SSR$beta
                       output$beta <- beta
                       
                       rmse <- list()
                       rmse$ASR <- mapply(i=1:100,function(i)
                         norm(rome_ASR$beta[,i] - b,"2")/norm(b,"2"))
                       rmse$SSR <- mapply(i=1:100,function(i)
                         norm(rome_SSR$beta[,i] - b,"2")/norm(b,"2"))
                       output$RMSE <- rmse
                       
                       violation <- list()
                       violation$ASR <- rowSums(viol_ASR)
                       violation$SSR <- rowSums(viol_SSR)
                       output$violation <- violation
                       
                       list_dummy[[iter]] <- output
                     }  
 stopCluster(cl)
 save(run_sim,file=paste0("violation_comparison_rho_type",rho_type,"_type",type,"_p",p,"_n",n,".RData"))
}
  
#####################################################################################
# Run simulations:
iter_max <- 100
n <- 100
p_vec <- c(100,500,1000)
type_vec <- c(5,6)
rho_vec <- c(1,2)

for (i in 1:3){
  p <- p_vec[i]
  
  for (j in 1:2){
    type <- type_vec[j]
    
    for (k in 1:2){
      rho_type <- rho_vec[k]
      
      cat("The current iteration is p:",p,"n:",n,"type is",type,"rho type is",rho_type,"\n")
      get_sim (n,p,type,rho_type,iter_max)
      cat("Done.\n")
      
    }
  }
}

#####################################################################################
# Draw plots:
library(ggplot2)
library(dplyr)
library(tidyr)

iter_max <- 100
n <- 100
p_vec <- c(100,500,1000)
type_vec <- c(5,6)
rho_vec <- c(1,2,3)

df_result <- data.frame()
for (i in 1:2){
  type <- type_vec[i]
  
  for (j in 1:2){
    rho_type <- rho_vec[j]
    
    for (k in 1:3){
      p <- p_vec[k]
      
      cat("The current iteration is p:",p,"n:",n,"type is",type,"rho type is",rho_type,"\n")
      load(paste0("violation_comparison_rho_type",rho_type,"_type",type,"_p",p,"_n",n,".RData"))
      
      df_result <- rbind(df_result, 
                         data.frame(
                           ASR = rowMeans(mapply(iter=1:iter_max,function(iter)
                             run_sim[[iter]]$violation$ASR))/p,
                           SSR = rowMeans(mapply(iter=1:iter_max,function(iter)
                             run_sim[[iter]]$violation$SSR))/p,
                           lambda = seq(1,100),
                           p = rep(p,100),
                           Corr_type = rep(ifelse(rho_type==1,"Type 1","Type 2"),100), 
                           Type = rep(ifelse(type==7,"AR","Block AR"),100)) )
    }
    
  }
}
df_long <- pivot_longer(data = df_result,
                        cols = c("ASR", "SSR"),
                        names_to = "Screening_rule",
                        values_to = "Value")
  

plot_violation <- ggplot(data = df_long, aes(x = lambda, y = Value, color = factor(p), linetype = Screening_rule)) +
  geom_line(size = 0.8) +
  facet_grid(Corr_type ~ Type, scales = "free_y") +
  scale_color_manual(values = c("100" = "red", "500" = "blue", "1000" = "green")) +
  scale_y_continuous( limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.1) ) +
  labs( title = "Proportion of Violations over Lambda Sequences",
        x = "Lambda sequence",
        y = "Proportion of violations",
        color = "Number of predictors",
        linetype = "Screening rules" ) +
  scale_x_continuous( breaks = seq(0, 100, by = 10)) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.ticks.x = element_blank(),  # remove ticks
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(size = 0.2), # finer vertical grid
    panel.grid.major.y = element_line(size = 0.2) )

pdf("violations.pdf", width = 13, height = 9)
plot_violation
dev.off()

