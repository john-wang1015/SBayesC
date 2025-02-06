library(ggplot2)
library(MCMCpack)

# noscale function
sbayesr = function(b, se, n, R, niter = 1000, gamma = c(0, 0.01, 0.1, 1), startPi = c(0.9, 0.06, 0.03, 0.01), startH2 = 0.5){
  m     = nrow(R)          # number of SNPs
  ndist = length(startPi)  # number of mixture distributions
  pi    = startPi          # starting value for pi
  h2    = startH2          # starting value for heritability
  
  bhat  = b                # Remove scale factor since data is already standardized
  vary  = 1                # Phenotypic variance = 1
  varg  = h2               # Genetic variance
  vare  = vary - varg      # Environmental variance should be 1 - h²
  
  sigmaSq = varg / (m * sum(gamma * pi))  # Ensure correct SNP variance scaling
  nub = 4  # Prior degrees of freedom for SNP effect variance
  nue = 4  # Prior degrees of freedom for residual variance
  scaleb = sigmaSq  # No need to transform further
  scalee = vare     # Keep environmental variance as computed
  
  beta = array(0, m)        # Vector of SNP effects
  beta_mcmc = matrix(0, niter, m) # MCMC samples of SNP effects
  bhatcorr = bhat
  probDelta = vector("numeric", ndist)
  logDelta = array(0, 2)
  keptIter = NULL
  
  for (iter in 1:niter){
    logPi = log(pi)
    logPiComp = log(1 - pi)
    invSigmaSq = 1 / (gamma * sigmaSq)
    logSigmaSq = log(gamma * sigmaSq)
    nsnpDist = rep(0, ndist)  # Number of SNPs in each effect distribution
    ssq = 0  # Sum of squares of beta
    nnz = 0  # Number of nonzero beta
    
    for (i in 1:m){
      oldSample = beta[i]
      
      # Fixing RHS and invLhs to remove unnecessary scaling
      rhs = bhatcorr[i] + oldSample
      invLhs = 1.0 / (1.0 + invSigmaSq)
      uhat = invLhs * rhs
      
      # Sampling mixture distribution membership
      logDelta = 0.5 * (log(invLhs) - logSigmaSq + uhat * rhs) + logPi
      logDelta[1] = logPi[1]
      for (k in 1:ndist) {
        probDelta[k] = 1.0 / sum(exp(logDelta - logDelta[k]))
      }
      
      delta = sample(1:ndist, 1, prob = probDelta)
      nsnpDist[delta] = nsnpDist[delta] + 1
      
      if (delta > 1) {
        beta[i] = rnorm(1, uhat[delta], sqrt(invLhs[delta]))
        bhatcorr = bhatcorr + R[, i] * (oldSample - beta[i])
        ssq = ssq + beta[i]^2 / gamma[delta]
        nnz = nnz + 1
      } else {
        if (oldSample) bhatcorr = bhatcorr + R[, i] * oldSample
        beta[i] = 0
      }
    }  
    
    beta_mcmc[iter, ] = beta  # No need to rescale beta
    
    # Sampling pi
    pi = rdirichlet(1, nsnpDist + 1)
    
    # Sampling SNP effect variance
    sigmaSq = (ssq + nub * scaleb) / rchisq(1, nnz + nub)
    
    # Compute genetic variance and heritability
    bRb = crossprod(beta, (bhat - bhatcorr))
    varg = bRb
    h2 = varg / (varg + vare)  # Correct heritability formula
    
    keptIter <- rbind(keptIter, c(pi, nnz, sigmaSq, h2))
    
    if (!(iter %% 100)) {
      cat(sprintf("\n iter %4s, nnz = %4s, sigmaSq = %6.3f, h2 = %6.3f\n", iter, nnz, sigmaSq, h2))
    }
  }
  
  colnames(keptIter) <- c(paste0("Pi", 1:length(pi)), "Nnz", "SigmaSq", "h2")
  postMean = apply(keptIter, 2, mean)
  postSD = apply(keptIter, 2, sd)
  cat("\nPosterior mean:\n")
  print(postMean)
  cat("\nPosterior standard deviation:\n")
  print(postSD)
  
  return(list(par = keptIter, beta = beta_mcmc))
}


