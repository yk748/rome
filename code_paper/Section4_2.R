#################################################################################
# This file is used to generate Figure A.1, described in Section 4.3.
# This code compares the proposed algorithm (CD) with gradient descent (GD) and 
# gradient coordinate descent (GCD) to numerically support the convergence 
# analysis in Section 2.3.
# It compares the objective values obtained by each estimator across a given 
# number of iterations.
#################################################################################
rm(list = ls())

library(latex2exp)
library(mvtnorm)

system("R CMD SHLIB convergence_huber.c")
if (.Platform$OS.type == "windows") {
  dyn.load("convergence_huber.dll")
} else {
  dyn.load("convergence_huber.so")
}

# ----------------------------------------------- #
# Utility functions:
huber_grad_vec <- function(v, delta) {
  ifelse(abs(v) <= delta, v, delta*sign(v))
}

gen_data <- function(n, p, type) {
  b <- (-1)^(1:p) * exp(-(0:(p - 1)) / 10)
  
  if (type == 5) {
    rho   <- 0.4
    Sigma <- outer(1:p, 1:p, function(j, k) rho^abs(j - k))
    X     <- rmvt(n, sigma = Sigma, df = 4)
  } else if (type == 6) {
    rho1 <- 0.2; rho2 <- 0.8; half <- p / 2
    S1   <- outer(1:half, 1:half, function(j, k) rho1^abs(j - k))
    S2   <- outer(1:half, 1:half, function(j, k) rho2^abs(j - k))
    X    <- cbind(rmvnorm(n, mean = rep(0, half), sigma = S1),
                  rmvnorm(n, mean = rep(0, half), sigma = S2))
  }
  y <- X %*% b + rnorm(n, 0, 1)
  list(X = X, y = y, b = b)
}

# ----------------------------------------------- #
# R wrappers:
run_gd <- function(X, y, delta, lambda, max_iter = 50){
  .Call("run_gd",  X, as.numeric(y),
        as.numeric(delta), as.numeric(lambda), as.integer(max_iter)) 
}

run_cgd <- function(X, y, delta, lambda, max_iter = 50){
  .Call("run_cgd", X, as.numeric(y),
        as.numeric(delta), as.numeric(lambda), as.integer(max_iter)) 
}

run_cd <- function(X, y, delta, lambda, max_iter = 50){
  .Call("run_cd",  X, as.numeric(y),
        as.numeric(delta), as.numeric(lambda), as.integer(max_iter)) 
}

# ----------------------------------------------- #
# Main loop:
run_sim <- function(max_iter = 50, delta = 0.5) {
  
  set.seed(123)
  
  settings <- list(
    list(p = 100,  n = 100, label = "p=100, n=100"),
    list(p = 500, n = 100, label = "p=500, n=100"),
    list(p = 1000, n = 500, label = "p=1000, n=500")
  )
  types    <- c(5, 6)
  type_lab <- c("5" = "AR", "6" = "Block AR")
  
  results <- list()
  for (type in types) {
    for (s in settings) {
      n <- s$n; p <- s$p; lab <- s$label
      cat("Running:", type_lab[as.character(type)],",", lab, "\n")
      
      dat <- gen_data(n, p, type)
      X <- dat$X; y <- dat$y
      lambda0 <- max(abs(t(huber_grad_vec(y, delta)) %*% X)/n)
      lambda <- 1e-3*lambda0
      cat("  lambda =", round(lambda, 6), "\n")
      
      cat("  GD...\n");  
      obj_gd  <- run_gd( X, y, delta, lambda, max_iter)
      cat("  CGD...\n"); 
      obj_cgd <- run_cgd(X, y, delta, lambda, max_iter)
      cat("  CD...\n");  
      obj_cd  <- run_cd( X, y, delta, lambda, max_iter)
      
      key <- paste0("type", type, "_", gsub("[= ,]", "_", lab))
      results[[key]] <- list(type    = type,
                             setting = lab,
                             lambda  = lambda,
                             GD      = obj_gd,
                             CGD     = obj_cgd,
                             CD      = obj_cd)
      cat("  Done.\n\n")
    }
  }
  results
}

# ----------------------------------------------- #
# Plotter:
plot_convergence <- function(results, max_iter = 50) {
  n_plots <- length(results)
  par(mfrow = c(ceiling(n_plots / 3), 3), mar = c(4, 4, 3, 1))
  
  for (key in names(results)) {
    res     <- results[[key]]
    vals    <- list(res$GD, res$CGD, res$CD)
    ylim    <- range(unlist(vals), na.rm = TRUE)
    ylim[1] <- max(ylim[1], 1e-10)
    
    plot(seq_len(max_iter), vals[[1]],
         type = "l", col = "blue", lwd = 2, lty = 3, log = "y",
         xlab = TeX("Iterations $(k)$"),
         ylab = "Objective values",
         ylim = ylim,
         main = paste0(res$setting, ", ",
                       ifelse(res$type == 5, "AR", "Block AR")))
    lines(seq_len(max_iter), vals[[2]], col = "green", lwd = 2, lty = 2)
    lines(seq_len(max_iter), vals[[3]], col = "red",  lwd = 2, lty = 1)
    legend("topright",
           legend = c("GD", "CGD", "CD"), bg = "transparent",
           col    = c("blue", "green", "red"),
           lwd = 2, lty = 3:1, bty = "n", cex = 0.8)
  }
}

#################################################################################
# Run:
results <- run_sim(max_iter = 50, delta = 0.5)

# ---------------------------------------------------------- #
# Figure A.1:
pdf("convergence.pdf", width = 11, height = 6)
plot_convergence(results, max_iter = 50)
dev.off()
