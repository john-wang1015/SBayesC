library(MASS)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript generate_simulations.R <output_folder>")
}
output_folder <- args[1]

simGWAS <- function(N, hsq, m, M, R) {
  idEffect <- sample(1:M, m)
  beta_x  <- rnorm(m, mean = 0, sd = sqrt(hsq / m))
  mu      <- crossprod(R[idEffect, ], beta_x)
  gwas_b  <- mvrnorm(n = 1, mu, Sigma = (1 - hsq) * R / N)
  gwas_se <- sqrt((1 - gwas_b * gwas_b) / (N - 1))
  data.frame(b = gwas_b, se = gwas_se, true_b = beta_x, pos = idEffect)
}

M <- 2000
r <- 0.9
ldm <- outer(1:M, 1:M, FUN = function(i, j) r^abs(i - j))

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

for (task_id in 1:100) {
  seed <- 1000 + task_id
  set.seed(seed)
  
  GWASss <- simGWAS(N = 1e7, hsq = 0.2, m = 100, M = M, R = ldm)
  
  gwas_file <- file.path(output_folder, paste0("GWASss_data", task_id, ".ma"))
  ldm_file  <- file.path(output_folder, paste0("ldm_data", task_id, ".ma"))
  
  write.table(GWASss, file = gwas_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(ldm, file = ldm_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  print(paste("Generated dataset:", task_id, "Seed:", seed))
}
