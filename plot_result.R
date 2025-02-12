#library(ggplot2)
par(mfrow = c(1, 2))

######## c++ code plot
read_matrix_bin <- function(filename) {
  con <- file(filename, "rb")
  rows <- readBin(con, "integer", 1, size = 4)
  cols <- readBin(con, "integer", 1, size = 4)
  data <- readBin(con, "numeric", rows * cols, size = 4)
  close(con)
  matrix(data, nrow = rows, ncol = cols, byrow = FALSE)  
}

beta_mcmc <- read_matrix_bin("ldm_data2_result.bin")
beta_means <- colMeans(beta_mcmc) 
gwas_data <- read.table("GWASss_data2.ma", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
b <- gwas_data['b']

index = 1:3000

true_beta = rep(0, 3000)
true_beta[gwas_data$pos[1:3]] = gwas_data$true_b[1:3]
df <- data.frame(b = b, beta_means = beta_means)  

plot(true_beta, t(beta_means), 
     col = "blue", pch = 16, cex = 1.2,  
     xlab = "b true", 
     ylab = "b", 
     main = "c++ code",
     ylim = c(-0.2,0.2))
abline(0, 1, col = "red", lty = 2, lwd = 2) 
# 
# plot(index, t(beta_means),
#      col = "blue", pch = 16, cex = 1.2,
#      xlab = "Index",
#      ylab = "b",
#      main = "c++ code",
#      ylim = c(-0.5,0.5))
# 
# points(gwas_data$pos[1:3], gwas_data$true_b[1:3],
#        col = "red", pch = 16, cex = 1.2)
# 
# legend("topright", legend = c("Estimated b", "True b"),
#        col = c("blue", "red"), pch = 16)

#########jz code results
beta_means_r <- read.table("betaMean_results_data2.ma", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
gwas_data_r <- read.table("GWASss_data2.ma", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
b_r <- gwas_data_r['b']

true_beta_r = rep(0, 3000)
true_beta_r[gwas_data_r$pos[1:3]] = gwas_data_r$true_b[1:3]

df <- data.frame(b = b_r, beta_means = beta_means_r)  # Create a data frame

plot(true_beta_r, t(beta_means_r), 
     col = "blue", pch = 16, cex = 1.2,  
     xlab = "b true", 
     ylab = "b", 
     main = "JZ R code",
     ylim = c(-0.5,0.5))
abline(0, 1, col = "red", lty = 2, lwd = 2)  # Dashed red line

# plot(index, t(beta_means), 
#      col = "blue", pch = 16, cex = 1.2,  
#      xlab = "Index", 
#      ylab = "b", 
#      main = "JZ R code",
#      ylim = c(-0.2,0.2))
# 
# 
# points(gwas_data$pos, gwas_data$true_b, 
#        col = "red", pch = 16, cex = 1.2)
# 
# legend("topright", legend = c("Estimated b", "True b"), 
#        col = c("blue", "red"), pch = 16)

#############################################################

par(mfrow = c(1, 2))

beta_mcmc <- read_matrix_bin("ldm_data1_result.bin")
beta_means <- colMeans(beta_mcmc) 
gwas_data <- read.table("GWASss_data1.ma", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
b <- gwas_data['b']

index = 1:1000

true_beta = rep(0, 1000)
true_beta[gwas_data$pos] = gwas_data$true_b
df <- data.frame(b = b, beta_means = beta_means)  

plot(true_beta, t(beta_means), 
     col = "blue", pch = 16, cex = 1.2,  
     xlab = "b true", 
     ylab = "b", 
     main = "SBayesC",
     ylim = c(-0.2,0.2))
abline(0, 1, col = "red", lty = 2, lwd = 2) 


beta_mcmc <- read_matrix_bin("ldm_data1_diff_prior_result.bin")
beta_means <- colMeans(beta_mcmc) 
df <- data.frame(b = b, beta_means = beta_means)  

plot(true_beta, t(beta_means), 
     col = "blue", pch = 16, cex = 1.2,  
     xlab = "b true", 
     ylab = "b", 
     main = "SBayesC diff prior",
     ylim = c(-0.2,0.2))
abline(0, 1, col = "red", lty = 2, lwd = 2) 
quartz.save("compare_m1_m2.png", type = "png", width = 15, height = 8)

#############################################################
par(mfrow = c(1, 2))
beta_mcmc <- read_matrix_bin("ldm_data1_result.bin")
beta_means <- colMeans(beta_mcmc) 
gwas_data <- read.table("GWASss_data1.ma", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
b <- gwas_data['b']

index = 1:1000

true_beta = rep(0, 1000)
true_beta[gwas_data$pos] = gwas_data$true_b
df <- data.frame(b = b, beta_means = beta_means) 

plot(index, t(beta_means),
     col = "blue", pch = 16, cex = 1.2,
     xlab = "Index",
     ylab = "b",
     main = "SBayesC",
     ylim = c(-0.2,0.2))

points(gwas_data$pos, gwas_data$true_b,
       col = "red", pch = 16, cex = 1.2)

legend("topright", legend = c("Estimated b", "True b"),
       col = c("blue", "red"), pch = 16)


beta_mcmc <- read_matrix_bin("ldm_data1_diff_prior_result.bin")
beta_means <- colMeans(beta_mcmc) 
df <- data.frame(b = b, beta_means = beta_means)  

plot(index, t(beta_means),
     col = "blue", pch = 16, cex = 1.2,
     xlab = "Index",
     ylab = "b",
     main = "SBayesC diff prior",
     ylim = c(-0.2,0.2))

points(gwas_data$pos, gwas_data$true_b,
       col = "red", pch = 16, cex = 1.2)

legend("topright", legend = c("Estimated b", "True b"),
       col = c("blue", "red"), pch = 16)

quartz.save("compare_m1_m2_index.png", type = "png", width = 15, height = 8)

#############################################################
par(mfrow = c(1, 1))
read_sigmaSq_nnz_from_bin <- function(file_path) {
  con <- file(file_path, "rb")  # Open file in binary mode
  
  # Read sigmaSq size
  sigmaSq_size <- readBin(con, what = "integer", size = 4, endian = "little")
  
  # Read sigmaSq values
  sigmaSq <- readBin(con, what = "numeric", size = 4, n = sigmaSq_size, endian = "little")
  
  # Read nnz size
  nnz_size <- readBin(con, what = "integer", size = 4, endian = "little")
  
  # Read nnz values
  nnz <- readBin(con, what = "numeric", size = 8, n = nnz_size, endian = "little")
  
  close(con)
  
  # Debug: Check if sizes match
  if (length(sigmaSq) != sigmaSq_size) {
    stop("Mismatch: sigmaSq size in file does not match number of values read.")
  }
  if (length(nnz) != nnz_size) {
    stop("Mismatch: nnz size in file does not match number of values read.")
  }
  
  # Ensure sizes match before creating a data frame
  min_size <- min(length(nnz), length(sigmaSq))
  df <- data.frame(nnz = nnz[1:min_size], sigmaSq = sigmaSq[1:min_size])
  
  return(df)
}

# Read files
df1 <- read_sigmaSq_nnz_from_bin("nnz_ssq_result1.bin")
df2 <- read_sigmaSq_nnz_from_bin("nnz_ssq_result1_diff_prior.bin")

# Verify sizes
nnz1 = df1$nnz[2000:10000]
nnz2 = df2$nnz[2000:10000]

# Compute densities
dens1 <- density(nnz1)
dens2 <- density(nnz2)

# Plot the first density curve
plot(dens1, col = "blue", lwd = 4, main = "Comparison of NNZ Distributions",
     xlab = "NNZ Values", ylab = "Density", ylim = range(c(dens1$y, dens2$y)))

# Add the second density curve
lines(dens2, col = "red", lwd = 4)

# Add a vertical line at x = 100 (true NNZ value)
abline(v = 100, col = "black", lwd = 4, lty = 2)  # Dashed black line

# Add legend
legend("topleft", legend = c("NNZ from Model 1", "NNZ from Model 2", "True NNZ = 100"), 
       col = c("blue", "red", "black"), lwd = 4, lty = c(1, 1, 2))
quartz.save("compare_m1_m2_nnz.png", type = "png", width = 10, height = 10)
