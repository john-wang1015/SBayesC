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
beta_means <- colMeans(beta_mcmc) # Compute mean of each column

gwas_data <- read.table("GWASss.ma", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
b <- gwas_data['b']

df <- data.frame(b = b, beta_means = beta_means)  # Create a data frame

ggplot(df, aes(x = b, y = beta_means)) +
  geom_point(color = "blue", size = 2) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1.2) +  
  labs(x = "Marginal SNP effect", y = "Estimate SNP effect", 
       title = "Scatter Plot") +
  xlim(-0.3, 0.3) +
  ylim(-0.3, 0.3) +
  theme_minimal()
