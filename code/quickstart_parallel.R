# ==============================================================================
# PARALLEL BOOTSTRAP WITH parallel (Base R)
# Quick Start Script
# ==============================================================================

library(agridat)
library(boot)
library(dplyr)
library(parallel)

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

# Create cluster
cl <- makeCluster(6)

# Export to cluster
clusterExport(cl, c("crossa.wheat", "analyze_one_genotype"))
clusterEvalQ(cl, {library(boot); library(dplyr)})

# Run analysis
genotypes <- unique(crossa.wheat$gen)[1:5]

results <- parLapply(
  cl,
  genotypes,
  function(g) analyze_one_genotype(g, crossa.wheat, n_boot = 5000)
)

# Stop cluster (IMPORTANT!)
stopCluster(cl)

# Convert to dataframe
results_df <- do.call(rbind, lapply(results, as.data.frame))
print(results_df)

cat("✅ Analysis complete using parallel!
")

