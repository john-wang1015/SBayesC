#!/usr/bin/env Rscript
library(MCMCpack)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript SBayesC.R <LD matrix file> <GWAS summary file> <output beta file> <output nnz file> <random_seed> <n_size> [n_iter]")
}

# Assign input arguments
ldm_file <- args[1]
gwas_file <- args[2]
output_beta_file <- args[3]
output_nnz_file <- args[4]
random_seed <- as.integer(args[5])
n_size <- as.numeric(args[6])
n_iter <- ifelse(length(args) > 6, as.integer(args[7]), 10000)

# Set the random seed for reproducibility
set.seed(random_seed)
cat("Using random seed:", random_seed, "\n")
cat("Using n_size:", n_size, "\n")

# Load data
gwas_data <- read.table(gwas_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
b <- gwas_data$b
se <- gwas_data$se

ldm <- read.table(ldm_file, sep = "\t", header = FALSE)
R <- as.matrix(ldm)

# Define SBayesR function
sbayesr <- function(b, se, n, R, niter = 10000, gamma = c(0, 1), startPi = c(0.9, 0.01), startH2 = 0.5){
  m <- nrow(R)          
  ndist <- length(startPi)  
  pi <- startPi         
  h2 <- startH2        
  scale <- sqrt(1 / (n * se^2)) 
  bhat <- b * scale     
  vary <- 1             
  varg <- h2           
  vare <- vary         
  sigmaSq <- varg / (m * sum(gamma * pi))    
  nub <- 4                  
  nue <- 4                  
  scaleb <- (nub - 2) / nub * sigmaSq  
  scalee <- (nue - 2) / nue * vare     
  beta <- rep(0, m)        
  beta_mcmc <- matrix(0, niter, m) 
  bhatcorr <- bhat
  probDelta <- numeric(ndist)
  logDelta <- numeric(2)
  keptIter <- NULL
  
  for (iter in 1:niter){
    logPi <- log(pi)
    logPiComp <- log(1 - pi)
    invSigmaSq <- 1 / (gamma * sigmaSq)
    logSigmaSq <- log(gamma * sigmaSq)
    nsnpDist <- rep(0, ndist)  
    ssq <- 0  
    nnz <- 0  
    
    for (i in 1:m){
      oldSample <- beta[i]
      rhs <- (bhatcorr[i] + oldSample) / (vare / n)
      invLhs <- 1.0 / (1 / (vare / n) + invSigmaSq)
      uhat <- invLhs * rhs
      
      logDelta <- 0.5 * (log(invLhs) - logSigmaSq + uhat * rhs) + logPi
      logDelta[1] <- logPi[1]
      for (k in 1:ndist) {
        probDelta[k] <- 1.0 / sum(exp(logDelta - logDelta[k]))
      }
      
      delta <- sample(1:ndist, 1, prob = probDelta)
      nsnpDist[delta] <- nsnpDist[delta] + 1
      
      if (delta > 1) {
        beta[i] <- rnorm(1, uhat[delta], sqrt(invLhs[delta]))
        bhatcorr <- bhatcorr + R[, i] * (oldSample - beta[i])
        ssq <- ssq + beta[i]^2 / gamma[delta]
        nnz <- nnz + 1
      } else {
        if (oldSample) bhatcorr <- bhatcorr + R[, i] * oldSample
        beta[i] <- 0
      }
    }	
    beta_mcmc[iter, ] <- beta / scale  
    
    pi <- rdirichlet(1, nsnpDist + 1)
    
    sigmaSq <- (ssq + nub * scaleb) / rchisq(1, nnz + nub)
    
    bRb <- crossprod(beta, (bhat - bhatcorr))
    varg <- bRb
    h2 <- varg / vary
    
    keptIter <- rbind(keptIter, c(pi, nnz, sigmaSq, h2))
    
    if (!(iter %% 100)) {
      cat(sprintf("\n iter %4s, nnz = %4s, sigmaSq = %6.3f, h2 = %6.3f\n", iter, nnz, sigmaSq, h2))
    }
  }
  
  colnames(keptIter) <- c(paste0("Pi", 1:length(pi)), "Nnz", "SigmaSq", "h2")
  postMean <- colMeans(keptIter)
  postSD <- apply(keptIter, 2, sd)
  cat("\nPosterior mean:\n")
  print(postMean)
  cat("\nPosterior standard deviation:\n")
  print(postSD)
  return(list(par = keptIter, beta = beta_mcmc))
}

# Run SBayesR analysis
res <- sbayesr(b, se, n_size, R, niter = n_iter)

# Save results
write.table(res$beta, file = output_beta_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(res$par, file = output_nnz_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("\nAnalysis complete. Results saved to:\n", output_beta_file, "\n", output_nnz_file, "\n")
