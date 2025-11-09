# Load required libraries
library(ggplot2)

# Function to read binary posterior mean vector (1000 x 1)
read_posterior_bin <- function(file_path) {
  con <- file(file_path, "rb")
  
  # Read number of rows (first 4 bytes, stored as an integer)
  num_rows <- readBin(con, what = "integer", size = 4, n = 1, endian = "little")
  
  # Read the posterior mean values (1000 floats)
  posterior_mean <- readBin(con, what = "numeric", size = 4, n = num_rows, endian = "little")
  
  close(con)
  return(posterior_mean)
}

# Function to read .ma files (tab-separated values)
read_matrix_ma <- function(filename) {
  as.matrix(read.table(filename, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
}

# Define dataset sizes and corresponding directories
dataset_sizes <- c("1e5", "1e6")  # Only 1e5 and 1e6 now
dataset_folders <- file.path("data_mix", dataset_sizes)

# Initialize list to store all results
all_results <- list()

# Loop over dataset folders
for (i in seq_along(dataset_sizes)) {
  folder <- dataset_folders[i]
  dataset_size <- dataset_sizes[i]
  
  message(sprintf("Processing dataset: %s (%d/%d)", dataset_size, i, length(dataset_sizes)))
  
  # Preallocate results for 100 iterations
  dataset_results <- vector("list", 100)
  
  # Iterate sequentially over n = 1 to 100
  for (n in 1:100) {
    # File paths
    beta_mcmc_file_r <- sprintf("%s/beta_r_results_%d_posterior_mean.bin", folder, n)
    beta_posterior_file <- sprintf("%s/ldm_data%d_result_posterior_mean.bin", folder, n)
    beta_posterior_2_file <- sprintf("%s/ldm_data%d_nzp_result_posterior_mean.bin", folder, n)
    beta_posterior_3_file <- sprintf("%s/ldm_data%d_I_result_posterior_mean.bin", folder, n)  # Identity matrix
    gwas_file <- sprintf("%s/GWASss_data3_mix_%d.ma", folder, n)  # Using GWAS summary 3 mix
    
    # Read posterior mean directly from binary files
    beta_means_r <- colMeans(read_matrix_ma(beta_mcmc_file_r))  # JZ_R data (from .ma file)
    beta_means <- read_posterior_bin(beta_posterior_file)  # SBayesC posterior mean
    beta_means_2 <- read_posterior_bin(beta_posterior_2_file)  # Diff Prior posterior mean
    beta_means_3 <- read_posterior_bin(beta_posterior_3_file)  # Identity LD matrix (R=I)
    
    # Read GWAS data
    gwas_data <- read.table(gwas_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    true_beta <- rep(0, 1000)  # Assume 1000 SNPs
    true_beta[gwas_data$pos] <- gwas_data$true_b
    
    # Compute correlations
    correlation_jz <- cor(true_beta, beta_means_r, use = "complete.obs")
    correlation_sbc <- cor(true_beta, beta_means, use = "complete.obs")
    correlation_diff_prior <- cor(true_beta, beta_means_2, use = "complete.obs")
    correlation_identity <- cor(true_beta, beta_means_3, use = "complete.obs")  # R=I case
    
    # Print progress for each iteration
    print(sprintf("Dataset: %s | Iteration: %d/100 completed.", dataset_size, n))
    
    # Store results
    dataset_results[[n]] <- data.frame(
      Dataset = dataset_size,
      Method = c("JZ_R", "SBayesC", "Diff Prior", "R=I"),
      Correlation = c(correlation_jz, correlation_sbc, correlation_diff_prior, correlation_identity)
    )
  }
  
  # Combine results for this dataset
  all_results[[dataset_size]] <- do.call(rbind, dataset_results)
}

# Combine all results into a single data frame
final_data <- do.call(rbind, all_results)

# Ensure factor ordering for ggplot
final_data$Method <- factor(final_data$Method, levels = c("JZ_R", "SBayesC", "Diff Prior", "R=I"))
final_data$Dataset <- factor(final_data$Dataset, levels = dataset_sizes)

# Plot using ggplot2 with white background
p <- ggplot(final_data, aes(x = Method, y = Correlation, fill = Method)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  facet_wrap(~ Dataset, nrow = 1, ncol = 2, scales = "free_y") +  # Adjusted layout
  theme_bw() +  # White background with black grid
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 6, face = "bold"),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA)  # Black border around the plot
  ) +
  labs(title = "Comparison of Correlations for 1e5 and 1e6", y = "Correlation")

# Save as PNG with white background
ggsave("correlation_boxplots_data_mix.png", p, width = 10, height = 5, dpi = 300, bg = "white")

# Display plot
print(p)
