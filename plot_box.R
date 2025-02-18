# Function to read a binary posterior mean vector (1000 x 1)
read_posterior_mean <- function(file_path) {
  con <- file(file_path, "rb")
  
  # Read number of rows (first 4 bytes, stored as an integer)
  num_rows <- readBin(con, what = "integer", size = 4, n = 1, endian = "little")
  
  # Read the posterior mean values (1000 floats)
  posterior_mean <- readBin(con, what = "numeric", size = 4, n = num_rows, endian = "little")
  
  close(con)
  return(posterior_mean)
}

# Function to read .ma files (assuming tab-separated values)
read_matrix_ma <- function(filename) {
  as.matrix(read.table(filename, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
}

# Define dataset sizes and corresponding directories
dataset_sizes <- c("1e4", "1e5", "1e6", "1e7", "1e8", "1e9")
dataset_folders <- c("data/small_samples/1e4", "data/small_samples/1e5", "data/small_samples/1e6",
                     "data/large_samples/1e7", "data/large_samples/1e8", "data/large_samples/1e9")

nn = 100  # Number of iterations per dataset

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
    # File paths to posterior mean binary files
    beta_mcmc_file_r <- sprintf("%s/beta_r_results_%d.ma", folder, n)
    beta_posterior_file <- sprintf("%s/ldm_data%d_result_posterior_mean.bin", folder, n)
    beta_posterior_2_file <- sprintf("%s/ldm_data%d_nzp_result_posterior_mean.bin", folder, n)
    gwas_file <- sprintf("%s/GWASss_data%d.ma", folder, n)
    
    # Read posterior mean directly
    beta_means_r <- colMeans(read_matrix_ma(beta_mcmc_file_r)) # JZ_R data (from .ma file)
    beta_means <- read_posterior_mean(beta_posterior_file) # SBayesC posterior mean
    beta_means_2 <- read_posterior_mean(beta_posterior_2_file) # Diff Prior posterior mean
    
    # Read GWAS data
    gwas_data <- read.table(gwas_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    true_beta <- rep(0, 1000)  # Assume 1000 SNPs
    true_beta[gwas_data$pos] <- gwas_data$true_b
    
    # Compute correlations
    correlation_values_jz[n] <- cor(true_beta, beta_means_r, use = "complete.obs")
    correlation_values_sbc[n] <- cor(true_beta, beta_means, use = "complete.obs")
    correlation_values_2[n] <- cor(true_beta, beta_means_2, use = "complete.obs")
  }
  
  # Combine results for boxplot
  correlation_data <- rbind(
    data.frame(Method = "JZ_R", Correlation = correlation_values_jz),
    data.frame(Method = "SBayesC", Correlation = correlation_values_sbc),
    data.frame(Method = "Diff Prior", Correlation = correlation_values_2)
  )
  
  # Ensure correct ordering: JZ_R (LHS), SBayesC (Middle), Diff Prior (RHS)
  correlation_data$Method <- factor(correlation_data$Method, levels = c("JZ_R", "SBayesC", "Diff Prior"))
  
  # Create boxplot with dataset-specific title
  boxplot(Correlation ~ Method, data = correlation_data,
          main = paste("Dataset:", dataset_size),
          ylab = "Correlation",
          col = c("pink", "lightblue", "lightgreen"))     
}

# Reset plotting layout to default
par(mfrow = c(1, 1))
