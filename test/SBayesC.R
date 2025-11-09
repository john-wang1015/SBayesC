###############################################################
## Summary-data-based BayesR (SBayesR)
## Author: Jian Zeng
## Date: 16 June 2022
## Input data: 
##   b: marginal SNP effects (from GWAS using 0/1/2 genotypes)
##   se: standard errors
##   n: GWAS sample size
##   R: LD correlation matrix
###############################################################

library(MCMCpack)

sbayesr = function(b, se, n, R, niter = 5000, gamma = c(0, 1), startPi = c(0.9, 0.1), startH2 = 0.2){
  m     = nrow(R)          # number of SNPs
  ndist = length(startPi)  # number of mixture distributions
  pi    = startPi          # starting value for pi
  h2    = startH2          # starting value for heritability
  scale = sqrt(1/(n*se^2)) # scale factors for marginal SNP effects, i.e., sqrt(2pq/vary)
  bhat  = b*scale          # scaled marginal SNP effects (in units of genotype sd / phenotype sd)
  vary  = 1                # phenotypic variance = 1 due to the scaling
  varg  = h2               # genetic variance
  vare  = vary             # assuming each SNP effect is vanishingly small
  sigmaSq = varg/(m*sum(gamma*pi))    # common factor of SNP effect variance
  nub = 4                  # prior degrees of freedom for SNP effect variance
  nue = 4                  # prior degrees of freedom for residual variance
  scaleb = (nub-2)/nub*sigmaSq  # prior scale parameter for SNP effect variance
  scalee = (nue-2)/nue*vare     # prior scale parameter for residual variance
  beta = array(0,m)        # vector of SNP effects
  beta_mcmc = matrix(0,niter,m) # MCMC samples of SNP effects
  bhatcorr = bhat
  probDelta = vector("numeric", ndist)
  logDelta = array(0,2)
  keptIter = NULL
  
  for (iter in 1:niter){
    # sampling SNP effects
    logPi = log(pi)
    logPiComp = log(1-pi)
    invSigmaSq = 1/(gamma*c(sigmaSq))
    logSigmaSq = log(gamma*c(sigmaSq))
    nsnpDist = rep(0, ndist)  # number of SNPs in each effect distribution
    ssq = 0  # sum of squares of beta
    nnz = 0  # number of nonzero beta
    for (i in 1:m){
      oldSample = beta[i]
      rhs = (bhatcorr[i] + oldSample)/(vare/n)
      invLhs = 1.0/(1/c(vare/n) + invSigmaSq)
      uhat = invLhs*c(rhs)
      
      # sampling mixture distribution membership
      logDelta = 0.5*(log(invLhs) - logSigmaSq + uhat*c(rhs)) + logPi
      logDelta[1] = logPi[1];
      for (k in 1:ndist) {
        probDelta[k] = 1.0/sum(exp(logDelta - logDelta[k]))
      }
      
      delta = sample(1:ndist, 1, prob = probDelta)
      nsnpDist[delta] = nsnpDist[delta] + 1
      
      if (delta > 1) {
        beta[i] = rnorm(1, uhat[delta], sqrt(invLhs[delta]))
        bhatcorr = bhatcorr + R[,i]*(oldSample - beta[i])
        ssq = ssq + beta[i]^2 / gamma[delta]
        nnz = nnz + 1
      } else {
        if (oldSample) bhatcorr = bhatcorr + R[,i]*oldSample
        beta[i] = 0
      }
    }	
    beta_mcmc[iter,] = beta / scale  #  store beta in the original scale
    
    # sampling pi
    pi = rdirichlet(1, nsnpDist + 1)
    
    
    # sampling SNP effect variance
    sigmaSq = (ssq + nub*scaleb)/rchisq(1, nnz+nub)
    
    # compute genetic variance and heritability
    beta_unscaled = beta / scale
    varg = crossprod(beta_unscaled, (bhat - bhatcorr)/scale)
    
    h2  = varg/vary
    
    keptIter <- rbind(keptIter,c(pi, nnz, sigmaSq, h2))
    
    if (!(iter%%100)) {
      cat (sprintf("\n iter %4s, nnz = %4s, sigmaSq = %6.3f, h2 = %6.3f\n", iter, nnz, sigmaSq, h2))
    }
  }
  
  colnames(keptIter) <- c(paste0("Pi", 1:length(pi)),"Nnz","SigmaSq","h2")
  postMean = apply(keptIter, 2, mean)
  postSD = apply(keptIter, 2, sd)
  cat("\nPosterior mean:\n")
  print(postMean)
  cat("\nPosterior standard deviation:\n")
  print(postSD)
  return(list(par=keptIter, beta=beta_mcmc))
}


set.seed(1001)

## load genotype data (the first 1000 common SNPs on chromosome 22 in 1KPG data)
#load("1000G_eur_chr22_1ksnp.RData")
#nind = nrow(X)
#nsnp = ncol(X)

## simulate a trait
#h2 = 0.5
#ncv = 100
#cv = sample(1:nsnp, ncv)
#beta.cv = rnorm(ncv)
#g <- X[, cv] %*% beta.cv
#y <- g + rnorm(nind, 0, sqrt(var(g)*(1/h2 - 1)))

## split into GWAS data set and testing data set
#n_gwas = 300
#X_gwas = X[1:n_gwas,]
#X_test = X[(n_gwas+1):nind,]
#y_gwas = y[1:n_gwas]
#y_test = y[(n_gwas+1):nind]

## conduct GWAS analysis
#b  = apply(X_gwas, 2, function(x){summary(lm(y_gwas~x))$coef[2,1]})
#se = apply(X_gwas, 2, function(x){summary(lm(y_gwas~x))$coef[2,2]})

## compute LD correlation matrix
#R = cor(X)

## load SBayesR script
#source("sbayesr.R")
gwas_data <- read.table("GWASss_data4.ma", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
b <- gwas_data$b
se <- gwas_data$se

ldm <- read.table("ldm_data4.ma", sep = "\t", header = FALSE)

R <- as.matrix(ldm)

## run SBayesR analysis
res = sbayesr(b, se, 1000, R)

## you can check the posterior distribution of SNP effects (obtained from MCMC samples)
plot(density(res$beta[,1]))

## compute the posterior mean values as the SNP effect estimates
betaMean = colMeans(res$beta)
plot(betaMean)

## less noisy compared to the GWAS effect estimates
plot(b)

## check the effect size plot
plot(b, betaMean)
abline(h=0, a=0, b=1)

# Save the posterior mean values as a .ma file
write.table(betaMean, 
            file = "betaMean_results_data2.ma", 
            sep = "\t",           # Use tab-delimited format
            row.names = FALSE,    # Exclude row names
            col.names = FALSE,    # Exclude column names
            quote = FALSE)        # Exclude quotes around values

# Save the full MCMC samples of beta as a .ma file (optional)
write.table(res$beta, 
            file = "beta_mcmc_samples_data2.ma", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

## compute PGS for testing individuals
#pgs = X_test %*% betaMean

## evaluate the performance of PGS prediction
#summary(lm(y_test ~ pgs))$r.squared


