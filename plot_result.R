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

beta_mcmc <- read_matrix_bin("ldm_data1_result.bin")
beta_means <- colMeans(beta_mcmc) 
gwas_data <- read.table("GWASss.ma", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
b <- gwas_data['b']

index = 1:1000

true_beta = rep(0, 1000)
true_beta[gwas_data$pos] = gwas_data$true_b
df <- data.frame(b = b, beta_means = beta_means)  

plot(true_beta, t(beta_means), 
     col = "blue", pch = 16, cex = 1.2,  
     xlab = "b true", 
     ylab = "b", 
     main = "c++ code",
     ylim = c(-0.2,0.2))
abline(0, 1, col = "red", lty = 2, lwd = 2) 

plot(index, t(beta_means), 
     col = "blue", pch = 16, cex = 1.2,  
     xlab = "Index", 
     ylab = "b", 
     main = "c++ code",
     ylim = c(-0.2,0.2))

points(gwas_data$pos, gwas_data$true_b, 
       col = "red", pch = 16, cex = 1.2)

legend("topright", legend = c("Estimated b", "True b"), 
       col = c("blue", "red"), pch = 16)

#########jz code results
beta_means <- read.table("betaMean_results.ma", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
gwas_data <- read.table("GWASss.ma", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
b <- gwas_data['b']

true_beta = rep(0, 1000)
true_beta[gwas_data$pos] = gwas_data$true_b

df <- data.frame(b = b, beta_means = beta_means)  # Create a data frame

plot(true_beta, t(beta_means), 
     col = "blue", pch = 16, cex = 1.2,  
     xlab = "Index", 
     ylab = "b", 
     main = "JZ R code",
     ylim = c(-0.2,0.2))
abline(0, 1, col = "red", lty = 2, lwd = 2)  # Dashed red line

plot(index, t(beta_means), 
     col = "blue", pch = 16, cex = 1.2,  
     xlab = "Index", 
     ylab = "b", 
     main = "JZ R code",
     ylim = c(-0.2,0.2))

points(gwas_data$pos, gwas_data$true_b, 
       col = "red", pch = 16, cex = 1.2)

legend("topright", legend = c("Estimated b", "True b"), 
       col = c("blue", "red"), pch = 16)
