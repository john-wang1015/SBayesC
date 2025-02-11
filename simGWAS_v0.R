library(MASS)
simGWAS <- function(N=1e5, # Sample size
                    hsq=0.2, # heritability
                    m, # number of causal variants
                    M, # total number of SNPs
                    R=R){ # LD matrix - MxM matrix
  ## sampple causal effects
  ## note that m should be smaller than M
  idEffect <- sample(1:M,m)
  
  ## y = X * beta + e
  ## GWAS => (1/N)X' multiplication
  beta_x  <- rnorm(m,mean=0,sd=sqrt(hsq/m)) # each SNP explains the same proportion of heritability: hsq/m 
  mu      <- crossprod(R[idEffect,],beta_x)
  gwas_b  <- mvrnorm(n=1,m=mu,Sigma = (1-hsq)*R/N) ## estimated effect sizes
  gwas_se <- sqrt((1-gwas_b*gwas_b)/(N-1)) ## standard error
  results <- cbind.data.frame(b=gwas_b,se=gwas_se,true_b=beta_x,pos=idEffect)
  return(results)
} 

## Example of LD matrix
M <- 3000
r <- 0.9 # auto-regressive correlation : corr(X_i,X_j) = r^|i-j|
ldm <- outer(1:M,1:M,FUN = function(i,j) r**abs(i-j))

tt <- system.time( GWASss <- simGWAS(N=1e5,hsq=0.2,m=3,M=3000,R=ldm) )

write.table(GWASss, file = "GWASss_data2.ma", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#save(ldm, file = "ldm_data.RData")

write.table(ldm, file = "ldm_data2.ma", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ldm, file = "ldm_data2.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
