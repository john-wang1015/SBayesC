# Function to read binary matrix
read_matrix_bin <- function(filename) {
  con <- file(filename, "rb")
  rows <- readBin(con, "integer", 1, size = 4)
  cols <- readBin(con, "integer", 1, size = 4)
  data <- readBin(con, "numeric", rows * cols, size = 4)
  close(con)
  matrix(data, nrow = rows, ncol = cols, byrow = FALSE)  
}

# Function to read .ma files (assuming tab-separated values)
read_matrix_ma <- function(filename) {
  as.matrix(read.table(filename, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
}

# Define dataset sizes and corresponding directories
dataset_sizes <- c("1e4", "1e5", "1e6", "1e7", "1e8", "1e9")
dataset_folders <- c("data/small_samples/1e4", "data/small_samples/1e5", "data/small_samples/1e6",
                     "data/large_samples/1e7", "data/large_samples/1e8", "data/large_samples/1e9")

# Save plot as PNG
#png("correlation_boxplots.png", width = 1200, height = 800)

nn = 100

# Set up 2x3 plotting layout
par(mfrow = c(2, 3)) 

# Loop over dataset folders
for (i in seq_along(dataset_sizes)) {
  folder <- dataset_folders[i]
  dataset_size <- dataset_sizes[i]
  
  correlation_values_jz <- numeric(nn)
  correlation_values_sbc <- numeric(nn)
  correlation_values_2 <- numeric(nn)
  
  # Iterate over n = 1 to 100
  for (n in 1:nn) {
    # File paths
    beta_mcmc_file_r <- sprintf("%s/beta_r_results_%d.ma", folder, n)
    beta_mcmc_file <- sprintf("%s/ldm_data%d_result.bin", folder, n)
    gwas_file <- sprintf("%s/GWASss_data%d.ma", folder, n)
    beta_mcmc_2_file <- sprintf("%s/ldm_data%d_nzp_result.bin", folder, n)
    
    # Read data
    beta_mcmc <- read_matrix_bin(beta_mcmc_file)
    beta_means <- colMeans(beta_mcmc)
    
    beta_mcmc_2 <- read_matrix_bin(beta_mcmc_2_file)
    beta_means_2 <- colMeans(beta_mcmc_2)
    
    beta_mcmc_r <- read_matrix_ma(beta_mcmc_file_r)  # Read JZ_R data
    beta_means_r <- colMeans(beta_mcmc_r)            # Compute posterior mean
    
    # Read GWAS data
    gwas_data <- read.table(gwas_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    b <- gwas_data$b
    
    true_beta <- rep(0, 1000)  # Assume 1000 SNPs
    true_beta[gwas_data$pos] <- gwas_data$true_b
    
    # Compute correlations
    correlation_values_jz[n] <- cor(true_beta, beta_means_r, use = "complete.obs")
    correlation_values_sbc[n] <- cor(true_beta, beta_means, use = "complete.obs")
    correlation_values_2[n] <- cor(true_beta, beta_means_2, use = "complete.obs")
  }
  
  # Combine results into a data frame for boxplot
  correlation_data <- rbind(
    data.frame(Method = "JZ_R", Correlation = correlation_values_jz),
    data.frame(Method = "SBayesC", Correlation = correlation_values_sbc),
    data.frame(Method = "Diff Prior", Correlation = correlation_values_2)
  )
  
  # Ensure correct ordering: JZ_R (LHS), SBayesC (Middle), Diff Prior (RHS)
  correlation_data$Method <- factor(correlation_data$Method, levels = c("JZ_R", "SBayesC", "Diff Prior"))
  
  # Create subplot with dataset-specific title
  boxplot(Correlation ~ Method, data = correlation_data,
          main = paste("Dataset:", dataset_size),
          ylab = "Correlation",
          col = c("pink", "lightblue", "lightgreen"))     
}

# Reset plotting layout to default and save the plot
par(mfrow = c(1, 1)) 
#dev.off()
