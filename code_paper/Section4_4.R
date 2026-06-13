#################################################################################
# This file is used to generate Figure 4 in Section 4.4.
# This code compares the adaptive sequential strong rule (ASR) and the 
# sequential strong rule (SSR) applied to the proposed algorithm, evaluating 
# the number of coordinates that violate the selection rules along the lambda path.
#################################################################################
# Main loops:
get_sim <- function(n,p,type,rho_type,iter_max){
  
  library("doParallel")
  library("dqrng")
  library("doRNG")
  
  available_cores <- detectCores()
  num_cores <- min(available_cores - 1, 25)
  num_cores <- max(num_cores, 1)
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  list_dummy <- list()
  run_sim <- foreach(iter = 1:iter_max,
                     .errorhandling = 'stop',
                     .packages = c("Matrix","mvtnorm","rome"),
                     .options.RNG = 123) %dopar% {
                     
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
                       
                       huber_grad_vec <- function(v, thresh) {
                         ifelse(abs(v) <= thresh, v, thresh * sign(v))
                       }
                       gamma <- 0.5
                       lambda0 <- max(abs(t(huber_grad_vec(y, gamma)) %*% x)/n)
                       lstep   <- log(0.05) / 100
                       lambdas <- lambda0 * exp(lstep * 1:100)
                       
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
                       
                       # -------------------------------------------------------- #
                       # c_j(lambda_l) for all j and l 
                       compute_c_path <- function(beta_path) {
                         L <- ncol(beta_path)
                         c_matrix <- matrix(0, nrow = L, ncol = p)
                         for (l in 1:L) {
                           r <- as.numeric(y - x %*% beta_path[, l])
                           c_matrix[l, ] <- -colSums(huber_grad_vec(r, gamma) * x) / n  
                         }
                         return(c_matrix)
                       }
                       c_path_ASR <- compute_c_path(rome_ASR$beta)
                       c_path_SSR <- compute_c_path(rome_SSR$beta)
                       
                       # -------------------------------------------------------- #
                       # M_ell and screening rule thresholds
                       M_l_ASR    <- rep(NA, 100)
                       M_l_ASR[1] <- 1
                       ASR_rule   <- rep(NA, 100)
                       ASR_rule[1] <- 2*lambdas[1] - lambda0  
                       SSR_rule    <- rep(NA, 100)
                       SSR_rule[1] <- 2*lambdas[1] - lambda0
                       
                       for (l in 2:100) {
                         lam_diff    <- lambdas[l-1] - lambdas[l] 
                         M_l_ASR[l]  <- max(abs(c_path_ASR[l-1, ] - c_path_ASR[l, ])) / lam_diff
                         ASR_rule[l] <- lambdas[l] + M_l_ASR[l-1] * (lambdas[l] - lambdas[l-1])
                         SSR_rule[l] <- 2*lambdas[l] - lambdas[l-1]
                       }
                       
                       # -------------------------------------------------------- #
                       # Violations: predictor j violates at step l if
                       # |c_j(lambda_l)| >= rule_l  (was screened out but should be active)
                       kkt_viol_ASR <- matrix(FALSE, 100, p)
                       kkt_viol_SSR <- matrix(FALSE, 100, p)
                       
                       for (l in 2:100) {
                         # Step 1: Identify which predictors were screened OUT by the rules based on the PREVIOUS lambda step
                         screened_out_ASR <- abs(c_path_ASR[l-1, ]) < ASR_rule[l]
                         screened_out_SSR <- abs(c_path_SSR[l-1, ]) < SSR_rule[l]
                         
                         # Step 2: A violation occurs if a predictor was screened out, 
                         # but its gradient at the CURRENT step violates the KKT condition (meaning it should be active)
                         kkt_viol_ASR[l, ] <- screened_out_ASR & (abs(c_path_ASR[l, ]) > lambdas[l])
                         kkt_viol_SSR[l, ] <- screened_out_SSR & (abs(c_path_SSR[l, ]) > lambdas[l])
                       }
                       
                       # Sum the violations per lambda step
                       viol_ASR <- rowSums(kkt_viol_ASR)
                       viol_SSR <- rowSums(kkt_viol_SSR)
                       
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
                       violation$ASR <- viol_ASR
                       violation$SSR <- viol_SSR
                       output$violation <- violation
                       
                       list_dummy[[iter]] <- output
                     }  
 stopCluster(cl)
 save(run_sim,file=paste0("violation_comparison_rho_type",rho_type,"_type",type,"_p",p,"_n",n,".RData"))
}
  
#################################################################################
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

#################################################################################
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
                           Type = rep(ifelse(type==5,"AR","Block AR"),100)) )
    }

  }
}
df_long <- pivot_longer(data = df_result,
                        cols = c("ASR", "SSR"),
                        names_to = "Screening_rule",
                        values_to = "Value")


# ----------------------------------------------- #
# Figure 4:
plot_violation <- ggplot(data = df_long, aes(x = lambda, y = Value, color = factor(p), linetype = Screening_rule)) +
  geom_line(linewidth = 0.8) +
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
    panel.grid.major.x = element_line(size = 0.2),
    panel.grid.major.y = element_line(size = 0.2) )

pdf("violations.pdf", width = 9, height = 5)
plot_violation
dev.off()

