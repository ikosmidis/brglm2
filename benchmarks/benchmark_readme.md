# BRGLM2 Benchmarking Suite

Comprehensive performance benchmarking framework for comparing original and optimized versions of the brglm2 R package.

## Overview

This suite compares the performance of two versions of brglm2:
- **Original**: Main branch implementation
- **New**: Optimized implementation with performance improvements

The benchmarks test across multiple datasets (small, medium, and large) and different fitting methods (brglmFit and mdyplFit).

## Repository Structure
```
parent-directory/
├── brglm2/                      # Your working branch 
│   ├── benchmarks/
│   │   ├── benchmark_readme.md  # You're here (if this wasnt obvious)
│   │   ├── setup.R              # Configuration and setup
│   │   ├── check_setup.R        # Validation script
│   │   ├── main_benchmark.R     # Core benchmarking tests
│   │   ├── visualize_results.R  # Plot generation
│   │   ├── run_all.R            # Master execution script
│   │   └── results/             # Output directory (auto-created)
│   │       └── YYYYMMDD_HHMMSS/ # Timestamped results
│   │           ├── benchmark_results.txt
│   │           ├── benchmark_data.RData
│   │           └── plots/
│   ├── R/                       # Your source code
│   ├── tests/
│   └── ...
│
├── brglm2-original/             # Clone of main branch (reference)
│   └── ...
│
└── benchmark_libs/              # Installed package versions
    ├── original/
    │   └── brglm2/              # Compiled original version
    │       ├── R/
    │       ├── data/
    │       ├── libs/
    │       ├── Meta/
    │       ├── DESCRIPTION
    │       └── ...
    └── new/
        └── brglm2/              # Compiled optimized version
            ├── R/
            ├── data/
            ├── libs/
            ├── Meta/
            ├── DESCRIPTION
            └── ...
```

The `benchmark_libs` folder contains **installed (compiled) R packages**. This is essential for dual loading for comparison.

## Initial Setup

### 1. Directory Structure Setup

Starting from your parent directory (containing all projects):

```bash
# Navigate to parent directory
cd /path/to/parent-directory

# If not already present, create benchmark_libs
mkdir -p benchmark_libs/original
mkdir -p benchmark_libs/new

# Clone main branch for reference (if not already done)
git clone <repository-url> brglm2-original
cd brglm2-original
git checkout main
cd ..

# Your working branch (if not already present)
git clone <repository-url> brglm2
cd brglm2
git checkout -b benchmark-testing  # or your optimization branch

# Create results directory
mkdir -p benchmarks/results
```

### 2. Install Package Versions to `benchmark_libs`

From the parent directory:

```bash
# Install original version from main branch clone
cd brglm2-original
R CMD INSTALL . --library=../benchmark_libs/original

# Install new version from your working branch
cd ../brglm2
R CMD INSTALL . --library=../benchmark_libs/new
```

**Important**: After making changes to your optimized code, reinstall to `benchmark_libs/new`:

```bash
cd brglm2
R CMD INSTALL . --library=../benchmark_libs/new
```

### 3. Configure Paths

Edit `brglm2/benchmarks/setup.R` to point to the correct locations.

**If running from `brglm2/` directory** (recommended), use relative paths:

```r
# Paths relative to brglm2/ directory
LIB_ORIGINAL <- "../benchmark_libs/original"
LIB_NEW <- "../benchmark_libs/new"
```

**Or use absolute paths**:

```r
# Example absolute paths (adjust for your system)
LIB_ORIGINAL <- "C:/path/to/parent-directory/benchmark_libs/original"
LIB_NEW <- "C:/path/to/parent-directory/benchmark_libs/new"
```

Also update the package path in `brglm2/benchmarks/check_setup.R`:

```r
# Absolute path to your working brglm2 source
pkg_path <- "C:/path/to/parent-directory/brglm2"
```

### 4. Install Required R Packages

```r
install.packages(c("tictoc", "rbenchmark", "microbenchmark", 
                   "ggplot2", "tinytest"))
```

### 5. Verify Setup

From the `brglm2/` directory:

```r
setwd("C:/path/to/parent-directory/brglm2")  # Set working directory
source("benchmarks/check_setup.R")
```

This will verify:
- R version compatibility
- Required packages installed
- Library paths exist and contain brglm2
- Test datasets accessible
- Results directory structure

## Running Benchmarks

### Quick Start

**From the `brglm2/` directory** (your working branch):

```r
# Set working directory (if not already there)
setwd("C:/path/to/parent-directory/brglm2")

# Run complete benchmark suite
source("benchmarks/run_all.R")
```

This executes:
1. Setup validation
2. Test suite execution (both versions)
3. Microbenchmark comparisons (100 evaluations)
4. Large dataset tests
5. MDYPL method comparison
6. Visualization generation

Expected runtime: ~5 minutes (default settings)

### Individual Components

Run specific benchmark sections:

```r
# Setup and configuration
source("benchmarks/setup.R")

# Main benchmarks only
source("benchmarks/main_benchmark.R")

# Generate plots from existing results
source("benchmarks/visualize_results.R")
```

### Workflow After Code Changes

When you modify your optimized code:

```bash
# 1. Reinstall the new version
cd /path/to/parent-directory/brglm2
R CMD INSTALL . --library=../benchmark_libs/new

# 2. Re-run benchmarks
# In R:
setwd("C:/path/to/parent-directory/brglm2")
source("benchmarks/run_all.R")
```

## Benchmark Components

### 1. Test Suite Validation
- Runs complete tinytest suite for both versions
- Verifies numerical correctness
- Reports pass/fail rates and execution time

### 2. Microbenchmark Tests

#### Small Dataset (Lizards)
- Dataset: 409 observations, binomial GLM
- Default: 100 evaluations
- Tests: brglmFit method

#### Medium Dataset (Endometrial)
- Dataset: 79 observations, probit link
- Default: 50 evaluations
- Tests: brglmFit method

### 3. Large Dataset Tests

#### MultipleFeatures Dataset
- Dataset: 2000 observations, 432 features
- High-dimensional classification problem
- Single timing test (computationally intensive)
- Tests: brglmFit and mdyplFit methods

### 4. Visualizations

Generated plots:
- `01_lizards_comparison.png` - Boxplot comparison
- `02_endometrial_comparison.png` - Boxplot comparison
- `03_speedup_summary.png` - Speedup factor bar chart
- `04_absolute_timing.png` - Absolute execution times
- `05_lizards_distribution.png` - Distribution densities

## Increasing Evaluation Counts

To improve consistency and reduce variability:

### In `main_benchmark.R`

**Lizards microbenchmark** (Line ~68):
```r
lizards_bench <- microbenchmark(
  # ... 
  times = 1000,  # Change to increase/decrease test numbers
  unit = "s"
)
```

**Endometrial microbenchmark** (Line ~91):
```r
endo_bench <- microbenchmark(
  # ...
  times = 500,  # Change to increase/decrease test numbers
  unit = "s"
)
```

**MultipleFeatures tests** (Lines ~138-162):
For more robust timing on large datasets, increase replications.

```r
time_mf_orig <- system.time({
  replicate(5, {  # Change from 5 replications
    fit_mf_orig <- glm(full_mf_fm, data = MultipleFeatures, 
                       family = binomial(),
                       method = brglm2_original, 
                       subset = training, 
                       maxit = 200)
  })
})

# Divide elapsed time by number of replications for average
cat("Average time per run:", time_mf_orig["elapsed"] / 5, "s\n")
```

**MDYPL tests** (Lines ~195-221):
Similarly change replications:

```r
time_mdypl_orig <- system.time({
  replicate(5, {  # Change from 5 replications
    fit_mdypl_orig <- glm(full_mf_fm, data = MultipleFeatures,
                          family = binomial(),
                          method = mdypl_original,
                          alpha = alpha_val,
                          subset = training,
                          maxit = 200)
  })
})
```

### Recommended Settings for 10-Minute Runtime

```r
# Microbenchmarks
lizards:     times = 400 
endometrial: times = 200  

# Large dataset tests
MultipleFeatures (brglmFit): replicate(5, ...) 
MDYPL:                       replicate(5, ...) 
```

## Output Files

### Results Directory Structure

Each run creates a timestamped directory:

```
benchmarks/results/YYYYMMDD_HHMMSS/
├── benchmark_results.txt      # Complete text output
├── benchmark_data.RData       # R objects for plotting
└── plots/
    ├── 01_lizards_comparison.png
    ├── 02_endometrial_comparison.png
    ├── 03_speedup_summary.png
    ├── 04_absolute_timing.png
    └── 05_lizards_distribution.png
```

### Key Metrics Reported

- **Speedup factors**: New vs Original execution time ratios
- **Absolute timings**: Median/mean execution times
- **Test results**: Pass/fail counts for validation
- **Numerical accuracy**: Coefficient and deviance comparisons

## Troubleshooting

### Common Issues

**"Original library path does not exist"**
- Update `LIB_ORIGINAL` in `benchmarks/setup.R`
- Ensure brglm2 is installed: `R CMD INSTALL . --library=path/to/original`

**"brglm2 package NOT found"**
- Reinstall packages in both library locations
- Verify installation: `library(brglm2, lib.loc="path/to/lib")`

**"Package directory does NOT exist"**
- Update `pkg_path` in `benchmarks/check_setup.R` and `main_benchmark.R`
- Ensure you're running from project root

**Plots fail to generate**
- Check that `benchmark_data.RData` was created
- Verify ggplot2 is installed
- Review error messages in console output

### Verifying Results

Successful benchmarks should show:
- All tests passing (or matching pass rate between versions)
- Speedup factors > 1.0 (new faster than original)
- Numerical accuracy: coefficient differences < 1e-10

## Notes

- The test suite execution times are included in overall speedup metrics
- Large dataset tests use single evaluations due to computational cost
- All visualizations use consistent color scheme: Red (original), Green (new)
- Benchmark data is saved for regenerating plots without re-running tests

## Citations

- Kosmidis, I., & Firth, D. (2021). Jeffreys-prior penalty, finiteness and shrinkage in binomial-response generalized linear models. *Biometrika*, 108(1), 71-82.
- Sterzinger, P., & Kosmidis, I. (2024). An iteratively reweighted least squares algorithm for the maximum Diaconis-Ylvisaker prior penalized likelihood in binomial-response generalized linear models.

## License

Same as parent brglm2 package.