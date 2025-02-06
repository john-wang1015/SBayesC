library(ggplot2)


read_matrix_bin <- function(filename) {
  con <- file(filename, "rb")
  rows <- readBin(con, "integer", 1, size = 4)
  cols <- readBin(con, "integer", 1, size = 4)
  data <- readBin(con, "numeric", rows * cols, size = 4)
  close(con)
  matrix(data, nrow = rows, ncol = cols, byrow = FALSE)  # Column-major order
}

beta_mcmc <- read_matrix_bin("ldm_data1_result.bin")
beta_means <- colMeans(beta_mcmc)  # Compute mean of each column

gwas_data <- read.table("GWASss.ma", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
b <- gwas_data['b']

df <- data.frame(b = b, beta_means = beta_means)  # Create a data frame

ggplot(df, aes(x = b, y = beta_means)) +
  geom_point(color = "blue", size = 2) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1.2) +  
  labs(x = "Marginal SNP effect", y = "Estimate SNP effect", 
       title = "Scatter Plot") +
  xlim(-0.2, 0.2) +
  ylim(-0.2, 0.2) +
  theme_minimal()

# if no ggplot
plot(df$b, df$beta_means, 
     col = "blue", pch = 16, cex = 1.2,  # Blue points, solid circles
     xlab = "Marginal SNP effect", 
     ylab = "Estimate SNP effect", 
     xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2),
     main = "Scatter Plot")

abline(0, 1, col = "red", lty = 2, lwd = 2)  # Dashed red line
