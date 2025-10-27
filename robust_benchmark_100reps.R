# ==============================================================================
# ROBUST BENCHMARKING WITH 100 REPETITIONS
# Comparing quartiles of runtime
# ==============================================================================
# Author: Mitra Aftabi
# Supervisor: Dr. Andreas Busjahn
# Institution: CQ Beratung+Bildung GmbH - Berlin
# Date: October 2025
# ==============================================================================

library(agridat)
library(boot)
library(tidyverse)
library(furrr)
library(future)
library(future.apply)
library(parallel)
library(ggplot2)

# Setup
data(crossa.wheat)
use_cores <- detectCores() - 2
genotypes_to_analyze <- unique(crossa.wheat$gen)[1:5]

cat("=== ROBUST BENCHMARKING SETUP ===\n")
cat("CPU cores available:", detectCores(), "\n")
cat("Using cores:", use_cores, "\n")
cat("Genotypes to analyze:", length(genotypes_to_analyze), "\n")
cat("Bootstrap replicates:", 5000, "\n")
cat("Repetitions per method:", 100, "\n\n")

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

# ==============================================================================
# RUN EACH METHOD 100 TIMES
# ==============================================================================

n_reps <- 100
cat("Running", n_reps, "repetitions for each method...\n")
cat("This will take approximately 30-45 minutes.\n\n")

# Storage for timing results
times_serial <- numeric(n_reps)
times_furrr <- numeric(n_reps)
times_future <- numeric(n_reps)
times_parallel <- numeric(n_reps)

# METHOD 1: Serial (100 times)
cat("=== METHOD 1: SERIAL BASELINE ===\n")
cat("Testing Serial method (", n_reps, " reps)...\n")
for (i in 1:n_reps) {
  if (i %% 10 == 0) cat("  Rep", i, "/", n_reps, "\n")

  start_time <- Sys.time()
  results_serial <- lapply(genotypes_to_analyze, function(g) {
    analyze_one_genotype(g, crossa.wheat, n_boot = 5000)
  })
  end_time <- Sys.time()

  times_serial[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
}
cat("✅ Serial complete!\n\n")

# METHOD 2: furrr (100 times)
cat("=== METHOD 2: furrr ===\n")
cat("Testing furrr method (", n_reps, " reps)...\n")
plan(multisession, workers = use_cores)

for (i in 1:n_reps) {
  if (i %% 10 == 0) cat("  Rep", i, "/", n_reps, "\n")

  start_time <- Sys.time()
  results_furrr <- future_map(
    genotypes_to_analyze,
    ~ analyze_one_genotype(.x, crossa.wheat, n_boot = 5000),
    .options = furrr_options(seed = 123 + i)
  )
  end_time <- Sys.time()

  times_furrr[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
}
cat("✅ furrr complete!\n\n")

# METHOD 3: future (100 times)
cat("=== METHOD 3: future ===\n")
cat("Testing future method (", n_reps, " reps)...\n")

for (i in 1:n_reps) {
  if (i %% 10 == 0) cat("  Rep", i, "/", n_reps, "\n")

  start_time <- Sys.time()
  results_future <- future_lapply(
    genotypes_to_analyze,
    function(g) analyze_one_genotype(g, crossa.wheat, n_boot = 5000),
    future.seed = 123 + i
  )
  end_time <- Sys.time()

  times_future[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
}
cat("✅ future complete!\n\n")

# METHOD 4: parallel (100 times)
cat("=== METHOD 4: parallel ===\n")
cat("Testing parallel method (", n_reps, " reps)...\n")

for (i in 1:n_reps) {
  if (i %% 10 == 0) cat("  Rep", i, "/", n_reps, "\n")

  # Create cluster
  cl <- makeCluster(use_cores)
  clusterExport(cl, c("crossa.wheat", "analyze_one_genotype"))
  clusterEvalQ(cl, {
    library(boot)
    library(dplyr)
  })

  start_time <- Sys.time()
  results_parallel <- parLapply(
    cl,
    genotypes_to_analyze,
    function(g) analyze_one_genotype(g, crossa.wheat, n_boot = 5000)
  )
  end_time <- Sys.time()

  stopCluster(cl)

  times_parallel[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
}
cat("✅ parallel complete!\n\n")

# ==============================================================================
# CALCULATE QUARTILES
# ==============================================================================

cat("=== QUARTILE ANALYSIS ===\n\n")

calculate_quartiles <- function(times, method_name) {
  q <- quantile(times, probs = c(0.25, 0.5, 0.75))

  cat(method_name, ":\n")
  cat("  Q1 (25%):", round(q[1], 3), "seconds\n")
  cat("  Median (50%):", round(q[2], 3), "seconds\n")
  cat("  Q3 (75%):", round(q[3], 3), "seconds\n")
  cat("  Mean:", round(mean(times), 3), "seconds\n")
  cat("  SD:", round(sd(times), 3), "seconds\n")
  cat("  Min:", round(min(times), 3), "seconds\n")
  cat("  Max:", round(max(times), 3), "seconds\n\n")

  return(data.frame(
    Method = method_name,
    Q1 = q[1],
    Median = q[2],
    Q3 = q[3],
    Mean = mean(times),
    SD = sd(times),
    Min = min(times),
    Max = max(times)
  ))
}

# Calculate for all methods
results_serial_q <- calculate_quartiles(times_serial, "Serial")
results_furrr_q <- calculate_quartiles(times_furrr, "furrr")
results_future_q <- calculate_quartiles(times_future, "future")
results_parallel_q <- calculate_quartiles(times_parallel, "parallel")

# Combine results
all_quartiles <- rbind(
  results_serial_q,
  results_furrr_q,
  results_future_q,
  results_parallel_q
)

# ==============================================================================
# CALCULATE SPEEDUP BASED ON MEDIAN
# ==============================================================================

all_quartiles$Speedup_Median <- all_quartiles$Median[1] / all_quartiles$Median

cat("=== SPEEDUP COMPARISON (Based on Median) ===\n\n")
print(all_quartiles[, c("Method", "Median", "Speedup_Median", "Q1", "Q3")],
  row.names = FALSE, digits = 3
)

# ==============================================================================
# CREATE VISUALIZATIONS
# ==============================================================================

cat("\n=== CREATING VISUALIZATIONS ===\n")

# Prepare data for plotting
timing_data <- data.frame(
  Method = rep(c("Serial", "furrr", "future", "parallel"), each = n_reps),
  Time = c(times_serial, times_furrr, times_future, times_parallel)
)

# 1. BOXPLOT
cat("Creating boxplot...\n")
p_box <- ggplot(timing_data, aes(x = Method, y = Time, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(alpha = 0.2, width = 0.2) +
  labs(
    title = "Runtime Distribution Across 100 Repetitions",
    subtitle = paste("5 genotypes, 5000 bootstrap replicates,", use_cores, "cores"),
    x = "Method",
    y = "Time (seconds)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14)
  )

print(p_box)
ggsave("figures/runtime_distribution_boxplot.png",
  plot = p_box, width = 10, height = 6, dpi = 300
)
cat("✅ Boxplot saved!\n")

# 2. VIOLIN PLOT
cat("Creating violin plot...\n")
p_violin <- ggplot(timing_data, aes(x = Method, y = Time, fill = Method)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", size = 3, color = "red") +
  labs(
    title = "Runtime Distribution with Density (100 Repetitions)",
    subtitle = paste("Violin plot shows distribution shape -", use_cores, "cores used"),
    x = "Method",
    y = "Time (seconds)",
    caption = "Red point = median"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14)
  )

print(p_violin)
ggsave("figures/runtime_violin_plot.png",
  plot = p_violin,
  width = 10, height = 6, dpi = 300
)
cat("✅ Violin plot saved!\n")

# 3. SPEEDUP WITH ERROR BARS
cat("Creating speedup chart...\n")
speedup_data <- all_quartiles %>%
  mutate(
    Speedup_Q1 = all_quartiles$Median[1] / Q3, # Conservative
    Speedup_Median = all_quartiles$Median[1] / Median,
    Speedup_Q3 = all_quartiles$Median[1] / Q1 # Optimistic
  )

p_speedup <- ggplot(speedup_data, aes(x = Method, y = Speedup_Median)) +
  geom_col(aes(fill = Method), alpha = 0.7) +
  geom_errorbar(aes(ymin = Speedup_Q1, ymax = Speedup_Q3),
    width = 0.3, linewidth = 1
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  geom_text(aes(label = paste0(round(Speedup_Median, 2), "x")),
    vjust = -0.5, fontface = "bold", size = 4
  ) +
  labs(
    title = "Speedup with Quartile Range (100 Repetitions)",
    subtitle = "Error bars show Q1-Q3 range | Red line = baseline (1x)",
    x = "Method",
    y = "Speedup (×)",
    caption = "Higher is better"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14)
  )

print(p_speedup)
ggsave("figures/speedup_with_quartiles.png",
  plot = p_speedup,
  width = 10, height = 6, dpi = 300
)
cat("✅ Speedup chart saved!\n")

# 4. QUARTILE COMPARISON
cat("Creating quartile comparison...\n")
quartile_long <- all_quartiles %>%
  select(Method, Q1, Median, Q3) %>%
  pivot_longer(
    cols = c(Q1, Median, Q3),
    names_to = "Quartile",
    values_to = "Time"
  )

p_quartiles <- ggplot(quartile_long, aes(x = Method, y = Time, fill = Quartile)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(aes(label = round(Time, 2)),
    position = position_dodge(width = 0.9),
    vjust = -0.5, size = 3
  ) +
  labs(
    title = "Quartile Comparison Across Methods",
    subtitle = "Lower values indicate better performance",
    x = "Method",
    y = "Time (seconds)"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))

print(p_quartiles)
ggsave("figures/quartile_comparison.png",
  plot = p_quartiles,
  width = 10, height = 6, dpi = 300
)
cat("✅ Quartile comparison saved!\n")

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

cat("\n=== SAVING RESULTS ===\n")

# Create results directory if needed
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Save timing data
saveRDS(timing_data, "results/timing_data_100reps.rds")
saveRDS(all_quartiles, "results/quartile_comparison.rds")
saveRDS(speedup_data, "results/speedup_data.rds")

# Save as CSV for easy viewing
write.csv(all_quartiles, "results/quartile_comparison.csv", row.names = FALSE)
write.csv(timing_data, "results/timing_data_100reps.csv", row.names = FALSE)
write.csv(speedup_data, "results/speedup_comparison.csv", row.names = FALSE)

cat("\n✅ ALL ANALYSIS COMPLETE!\n\n")
cat("================================================================================\n")
cat("RESULTS SUMMARY\n")
cat("================================================================================\n\n")
cat("Files created:\n\n")
cat("DATA FILES:\n")
cat("  - results/timing_data_100reps.rds/csv\n")
cat("  - results/quartile_comparison.rds/csv\n")
cat("  - results/speedup_comparison.rds/csv\n\n")
cat("VISUALIZATIONS:\n")
cat("  - figures/runtime_distribution_boxplot.png\n")
cat("  - figures/runtime_violin_plot.png\n")
cat("  - figures/speedup_with_quartiles.png\n")
cat("  - figures/quartile_comparison.png\n\n")
cat("================================================================================\n")
cat(
  "FASTEST METHOD (by median):",
  all_quartiles$Method[which.min(all_quartiles$Median)], "\n"
)
cat(
  "MAXIMUM SPEEDUP:",
  round(max(all_quartiles$Speedup_Median), 2), "x\n"
)
cat("================================================================================\n\n")
cat("Analysis complete! Ready to share with Dr. Busjahn.\n")

# ==============================================================================
# END OF ROBUST BENCHMARKING
# ==============================================================================
