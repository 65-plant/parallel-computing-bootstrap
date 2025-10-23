# ==============================================================================
# PARALLEL COMPUTING FOR BOOTSTRAP ANALYSIS
# Master Analysis Script
# ==============================================================================
# Author: Mitra Aftabi
# Supervisor: Dr. Andreas Busjahn
# Institution: CQ Beratung+Bildung GmbH - Berlin
# Date: October 2025
# ==============================================================================

# ==============================================================================
# 1. SETUP
# ==============================================================================

# Set working directory
setwd("C:/Users/aftab/Documents/Weiterbildung-CQ/Internship")

# Load required libraries
library(agridat)      # Agricultural datasets
library(boot)         # Bootstrap methods
library(tidyverse)    # Data manipulation
library(furrr)        # Parallel with purrr syntax
library(future)       # Parallel framework
library(future.apply) # Future lapply functions
library(tictoc)       # Simple timing
library(parallel)     # Base R parallel

# Check available cores
cat("Available CPU cores:", detectCores(), "
")

# ==============================================================================
# 2. LOAD AND EXPLORE DATA
# ==============================================================================

# Load wheat dataset
data(crossa.wheat)

# Dataset summary
cat("
Dataset Summary:
")
cat("Observations:", nrow(crossa.wheat), "
")
cat("Genotypes:", length(unique(crossa.wheat$gen)), "
")
cat("Locations:", length(unique(crossa.wheat$loc)), "
")

# ==============================================================================
# 3. DEFINE ANALYSIS FUNCTION
# ==============================================================================

# Bootstrap confidence interval for one genotype
analyze_one_genotype <- function(genotype_name, data, n_boot = 5000) {
  
  # Filter data for this genotype
  geno_data <- data %>% filter(gen == genotype_name)
  
  # Bootstrap function (mean yield)
  boot_fun <- function(data, indices) {
    d <- data[indices, ]
    return(mean(d$yield))
  }
  
  # Run bootstrap
  boot_result <- boot(
    data = geno_data,
    statistic = boot_fun,
    R = n_boot
  )
  
  # Extract confidence interval
  ci <- boot.ci(boot_result, type = "perc")
  
  # Return results
  return(list(
    genotype = genotype_name,
    n_obs = nrow(geno_data),
    original_mean = boot_result$t0,
    boot_se = sd(boot_result$t),
    ci_lower = ci$percent[4],
    ci_upper = ci$percent[5]
  ))
}

# ==============================================================================
# 4. PREPARE ANALYSIS
# ==============================================================================

# Select genotypes to analyze
all_genotypes <- unique(crossa.wheat$gen)
genotypes_to_analyze <- all_genotypes[1:5]  # Adjust number as needed

cat("
Analyzing", length(genotypes_to_analyze), "genotypes
")
cat("Bootstrap replicates per genotype: 5000

")

# ==============================================================================
# 5. METHOD 1: SERIAL (BASELINE)
# ==============================================================================

cat("=== METHOD 1: SERIAL BASELINE ===
")

tic("Serial")
results_serial <- lapply(genotypes_to_analyze, function(g) {
  analyze_one_genotype(g, crossa.wheat, n_boot = 5000)
})
time_serial <- toc()

serial_time <- as.numeric(time_serial$toc - time_serial$tic)
cat("Time:", round(serial_time, 2), "seconds

")

# ==============================================================================
# 6. METHOD 2: furrr
# ==============================================================================

cat("=== METHOD 2: furrr (TIDYVERSE STYLE) ===
")

# Set up parallel processing
plan(multisession, workers = 6)

tic("furrr")
results_furrr <- future_map(
  genotypes_to_analyze,
  ~analyze_one_genotype(.x, crossa.wheat, n_boot = 5000),
  .options = furrr_options(seed = 123)
)
time_furrr <- toc()

furrr_time <- as.numeric(time_furrr$toc - time_furrr$tic)
speedup_furrr <- serial_time / furrr_time
cat("Time:", round(furrr_time, 2), "seconds
")
cat("Speedup:", round(speedup_furrr, 2), "x

")

# ==============================================================================
# 7. METHOD 3: future
# ==============================================================================

cat("=== METHOD 3: future (BEST PERFORMANCE) ===
")

# Use same parallel plan
tic("future")
results_future <- future_lapply(
  genotypes_to_analyze,
  function(g) analyze_one_genotype(g, crossa.wheat, n_boot = 5000),
  future.seed = 123
)
time_future <- toc()

future_time <- as.numeric(time_future$toc - time_future$tic)
speedup_future <- serial_time / future_time
cat("Time:", round(future_time, 2), "seconds
")
cat("Speedup:", round(speedup_future, 2), "x

")

# ==============================================================================
# 8. METHOD 4: parallel
# ==============================================================================

cat("=== METHOD 4: parallel (BASE R) ===
")

# Create cluster
cl <- makeCluster(6)

# Export data and functions
clusterExport(cl, c("crossa.wheat", "analyze_one_genotype"))

# Load packages on workers
clusterEvalQ(cl, {
  library(boot)
  library(dplyr)
})

tic("parallel")
results_parallel <- parLapply(
  cl,
  genotypes_to_analyze,
  function(g) analyze_one_genotype(g, crossa.wheat, n_boot = 5000)
)
time_parallel <- toc()

# Stop cluster
stopCluster(cl)

parallel_time <- as.numeric(time_parallel$toc - time_parallel$tic)
speedup_parallel <- serial_time / parallel_time
cat("Time:", round(parallel_time, 2), "seconds
")
cat("Speedup:", round(speedup_parallel, 2), "x

")

# ==============================================================================
# 9. COMPARE RESULTS
# ==============================================================================

cat("=== PERFORMANCE COMPARISON ===

")

comparison <- data.frame(
  Method = c("Serial", "furrr", "future", "parallel"),
  Time_sec = c(serial_time, furrr_time, future_time, parallel_time),
  Speedup = c(1.00, speedup_furrr, speedup_future, speedup_parallel)
)

print(comparison, row.names = FALSE)

cat("
Fastest method:", 
    comparison$Method[which.min(comparison$Time_sec)], "
")

# ==============================================================================
# 10. SAVE RESULTS
# ==============================================================================

# Create results directory
dir.create("results", showWarnings = FALSE)

# Save all results
saveRDS(results_serial, "results/results_serial_final.rds")
saveRDS(results_furrr, "results/results_furrr_final.rds")
saveRDS(results_future, "results/results_future_final.rds")
saveRDS(results_parallel, "results/results_parallel_final.rds")
saveRDS(comparison, "results/performance_comparison.rds")

cat("
✅ All results saved to results/ directory
")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================

