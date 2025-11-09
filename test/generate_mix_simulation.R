library(MASS)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript generate_simulations.R <size1> <size2> ...")
}

simGWAS <- function(N, hsq, m, M, R) {
  idEffect <- sample(1:M, m)
  beta_x  <- rnorm(m, mean = 0, sd = sqrt(hsq / m))
  mu      <- crossprod(R[idEffect, ], beta_x)
  gwas_b  <- mvrnorm(n = 1, mu, Sigma = (1 - hsq) * R / N)
  gwas_se <- sqrt((1 - gwas_b * gwas_b) / (N - 1))
  data.frame(b = gwas_b, se = gwas_se, true_b = beta_x, pos = idEffect)
}

M <- 1000
r1 <- 0.5  # First LD matrix correlation
r2 <- 0.9  # Second LD matrix correlation

# Generate two different LD matrices
ldm1 <- outer(1:M, 1:M, FUN = function(i, j) r1^abs(i - j))
ldm2 <- outer(1:M, 1:M, FUN = function(i, j) r2^abs(i - j))

for (size in args) {
  output_folder <- file.path("data_mix", size)
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  
  for (task_id in 1:100) {
    seed <- 1000 + task_id
    set.seed(seed)
    
    # Generate two GWAS summary statistics using different LD matrices
    GWASss1 <- simGWAS(N = as.numeric(size), hsq = 0.2, m = 100, M = M, R = ldm1)
    GWASss2 <- simGWAS(N = as.numeric(size), hsq = 0.2, m = 100, M = M, R = ldm2)
    
    # Create the third GWAS summary by interleaving rows (odd: GWASss1, even: GWASss2)
    GWASss3 <- GWASss1
    GWASss3[seq(2, nrow(GWASss3), by = 2), ] <- GWASss2[seq(2, nrow(GWASss2), by = 2), ]
    
    # Define file paths
    gwas_file1 <- file.path(output_folder, paste0("GWASss_data1_mix", task_id, ".ma"))
    gwas_file2 <- file.path(output_folder, paste0("GWASss_data2_mix", task_id, ".ma"))
    gwas_file3 <- file.path(output_folder, paste0("GWASss_data3_mix", task_id, ".ma"))
    
    ldm_file1  <- file.path(output_folder, paste0("ldm1_data_mix", task_id, ".ma"))
    ldm_file2  <- file.path(output_folder, paste0("ldm2_data_mix", task_id, ".ma"))
    
    # Save files
    write.table(GWASss1, file = gwas_file1, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    write.table(GWASss2, file = gwas_file2, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    write.table(GWASss3, file = gwas_file3, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    write.table(ldm1, file = ldm_file1, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(ldm2, file = ldm_file2, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    print(paste("Generated dataset:", size, "Task:", task_id, "Seed:", seed))
  }
}
