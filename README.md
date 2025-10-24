
# Parallel Computing for Plant Phenotype Bootstrap Analysis

**Author:** Mitra Aftabi  
**Supervisor:** Dr. Andreas Busjahn 
**Institution:** CQ Beratung+Bildung GmbH - Berlin 
**Date:** October 2025

---

## Project Overview

This project compares three parallel computing methods in R (furrr, future, parallel) 
for accelerating bootstrap confidence interval calculations in plant breeding data analysis.

### Key Achievement
**Up to 5.71x speedup** using the `future` package with optimal workload size!

---

## Quick Results

| Method | Time (5 genotypes, 5000 reps) | Speedup | Best Use Case |
|--------|-------------------------------|---------|---------------|
| Serial | 4.19s | 1.00x | Baseline |
| furrr | 2.71s | 1.55x | 5+ tasks, tidyverse users |
| future | 1.47s | 2.85x | 3-10 tasks (WINNER!) |
| parallel | 1.28s | 3.27x | No dependencies needed |

** Note :** 5 genotypes with `future` = **5.71x speedup!**

---

## Project Structure
```
Internship/
├── results/               # All benchmark data
│   ├── results_serial.rds
│   ├── results_furrr_5k.rds
│   ├── results_future_5k.rds
│   ├── results_parallel_5k.rds
│   ├── scalability_results.rds
│   └── performance_summary.csv
├── figures/              # Visualizations
│   ├── timing_comparison.png
│   ├── speedup_comparison.png
│   └── scalability_analysis.png
├── documentation/        # Guides and summaries
│   ├── how_to_guide_COMPLETE.txt
│   ├── installation_guide.txt
│   └── project_summary.txt
└── code/                # R scripts (if saved)
```

---

## Dataset

**Source:** crossa.wheat from agridat R package  
**Size:** 450 observations, 18 genotypes, 25 locations  
**Analysis:** Bootstrap confidence intervals for mean yield per genotype  
**Reference:** Crossa et al. (1990)

---

## Key Findings

### 1. Performance Depends on Workload Size
- **Small tasks (3-5):** `future` wins (low overhead)
- **Medium tasks (5-10):** `furrr` wins (better load balancing)
- **Large tasks (10+):** `furrr` maintains lead

### 2. Bootstrap Replicates Matter
- 1000 replicates: Minimal/no speedup (overhead dominates)
- 5000 replicates: Good speedup (2-5x)
- 10000 replicates: Best speedup (3-6x)

### 3. Different Methods for Different Needs
- **Beginners:** Use `furrr` (easiest syntax)
- **Performance:** Use `future` (best speedup)
- **No dependencies:** Use `parallel` (built-in)

---

## Installation
```r
# Install required packages
install.packages('furrr')
install.packages('future')
install.packages('future.apply')
install.packages('agridat')
install.packages('boot')
install.packages('tidyverse')

# parallel package comes with R (no install needed)
```

---

## Quick Start
```r
library(furrr)
library(boot)

# Set up parallel processing
plan(multisession, workers = 4)

# Analysis function
analyze_item <- function(item) {
  # Bootstrap code
  boot_result <- boot(data, statistic_function, R = 5000)
  return(boot_result)
}

# Run in parallel (automatically!)
results <- future_map(items, analyze_item)
```

---

## Documentation

- **Complete How-To Guide:** `documentation/how_to_guide_COMPLETE.txt`
- **Installation Guide:** `documentation/installation_guide.txt`
- **Project Summary:** `documentation/project_summary.txt`

---

## Hardware Used

- **CPU:** 8 cores (6 used for parallel processing)
- **RAM:** 16GB+
- **Platform:** Windows x64
- **R Version:** 4.5.1

---

## Recommendations

### For Production Use
1. **Best overall:** `future` (performance + flexibility)
2. **For tidyverse users:** `furrr` (familiar syntax)
3. **For HPC/clusters:** `parallel` (no dependencies)

### Optimization Tips
1. Use 5000+ bootstrap replicates for clear speedup
2. Leave 2 cores free for system (use `detectCores() - 2`)
3. Test with 3-10 items first before scaling up
4. Set random seed for reproducibility

---

## References

**Packages:**
- furrr: Vaughan & Dancho (2023)
- future: Bengtsson (2021)
- agridat: Wright (2023)

**Data:**
- Crossa, J. et al. (1990). Statistical analyses of multilocation trials. 
  Advances in Agronomy, 44, 55-85.

---

## License

Educational internship project - CQ Beratung+Bildung GmbH, Berlin

---

## Contact

**Student:** Mitra Aftabi  
**Email:** aftabi.mitra@gmail.com  
**Supervisor:** Dr. Andreas Busjahn

---

## Acknowledgments

- Dr. Andreas Busjahn for project supervision and guidance
- agridat package maintainers for providing datasets
- R Core Team for parallel computing infrastructure

---

**Last Updated:** October 2025


