# ==============================================================================
# PARALLEL BOOTSTRAP WITH future
# Quick Start Script
# ==============================================================================

library(agridat)
library(boot)
library(dplyr)
library(future)
library(future.apply)

# Load data
data(crossa.wheat)

# Analysis function
analyze_one_genotype <- function(genotype_name, data, n_boot = 5000) {
  geno_data <- data %>% filter(gen == genotype_name)
  
  boot_fun <- function(data, indices) {
    return(mean(data[indices, ]$yield))
  }
  
  boot_result <- boot(data = geno_data, statistic = boot_fun, R = n_boot)
  ci <- boot.ci(boot_result, type = "perc")
  
  return(list(
    genotype = genotype_name,
    mean = boot_result$t0,
    ci_lower = ci$percent[4],
    ci_upper = ci$percent[5]
  ))
}

# Set up parallel processing
plan(multisession, workers = 6)

# Run analysis
genotypes <- unique(crossa.wheat$gen)[1:5]

results <- future_lapply(
  genotypes,
  function(g) analyze_one_genotype(g, crossa.wheat, n_boot = 5000),
  future.seed = 123
)

# Convert to dataframe
results_df <- do.call(rbind, lapply(results, as.data.frame))
print(results_df)

cat("✅ Analysis complete using future!
")

