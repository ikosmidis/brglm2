# BRGLM2 Benchmarking Suite

Benchmarking framework for comparing original and optimised versions of the brglm2 R package.

## Overview

The suite compares two versions of brglm2:
- **Original**: Main branch implementation (reference)
- **New**: Optimised implementation under development

Tests cover small, medium, and large datasets using both `brglmFit` and `mdyplFit` methods.

---

## File Structure

```
parent-directory/
├── brglm2/                          # Your working branch
│   ├── benchmarks/
│   │   ├── benchmark_readme.md      # This file
│   │   ├── setup.R                  # Config, validation, CPU pinning, utilities
│   │   ├── main_benchmark.R         # Benchmark tests (sourced by run_all.R)
│   │   ├── visualize_results.R      # Plot generation (sourced by run_all.R)
│   │   ├── run_all.R                # Master execution script — run this
│   │   └── results/                 # Auto-created output directory
│   │       └── YYYYMMDD_HHMMSS/
│   │           ├── benchmark_results.txt
│   │           ├── benchmark_data.RData
│   │           └── plots/
│   ├── R/
│   ├── tests/
│   └── ...
│
├── brglm2-original/                 # Clone of main branch (reference)
│
└── benchmark_libs/
    ├── original/brglm2/             # Compiled original version
    └── new/brglm2/                  # Compiled optimised version
```

There is no separate `check_setup.R` — validation is now handled inside `setup.R`.

---

## Initial Setup

### 1. Directory Structure

From your parent directory:

```bash
mkdir -p benchmark_libs/original benchmark_libs/new

# Reference clone of main branch
git clone <repository-url> brglm2-original
cd brglm2-original && git checkout main && cd ..

# Your working branch (if not already present)
git clone <repository-url> brglm2
cd brglm2 && git checkout -b your-optimisation-branch
```

### 2. Install Both Package Versions

```bash
# Install original (reference) version
cd brglm2-original
R CMD INSTALL . --library=../benchmark_libs/original
# or
Rscript -e "install.packages('.', repos=NULL, type='source', lib='../benchmark_libs/new')"

# Install new (optimised) version
cd ../brglm2
R CMD INSTALL . --library=../benchmark_libs/new
# or
Rscript -e "install.packages('.', repos=NULL, type='source', lib='../benchmark_libs/new')"
```

After any code change, reinstall the new version before re-running benchmarks:

```bash
cd brglm2
R CMD INSTALL . --library=../benchmark_libs/new
```

### 3. Configure Paths

Open `benchmarks/setup.R` and update the two path variables at the top:

```r
LIB_ORIGINAL <- "../benchmark_libs/original"   # relative to brglm2/
LIB_NEW      <- "../benchmark_libs/new"
```

Absolute paths also work if you prefer:

```r
LIB_ORIGINAL <- "C:/path/to/parent-directory/benchmark_libs/original"
LIB_NEW      <- "C:/path/to/parent-directory/benchmark_libs/new"
```

### 4. Install Required R Packages

```r
install.packages(c("microbenchmark", "bench", "ggplot2", "tinytest"))
```

Note: `tictoc` and `rbenchmark` are no longer used.

---

## Running Benchmarks

### Option A: From the R Terminal (Easiest)
If you are already inside an R session (like RStudio or the R GUI), you can run the suite by manually setting the "Pinned" flag to bypass the auto-relaunch logic. Note that this runs without CPU pinning (slightly higher timing variance).

```r
setwd("C:/path/to/brglm2")
Sys.setenv(BRGLM2_PINNED = "1")
source("benchmarks/run_all.R")
```

### Option B: From PowerShell (Recommended for Accuracy)
To get the benefit of CPU pinning while ensuring your keyboard input (the `y/n` prompt) works correctly, use this two-step approach in your PowerShell terminal:

```powershell
# 1. Set the environment variable so setup.R doesn't try to relaunch
$env:BRGLM2_PINNED="1"

# 2. Run the script directly
Rscript.exe benchmarks\run_all.R
```

### Option C: Manual Pinning (CMD / Batch)
If you are using a standard Windows Command Prompt (not PowerShell), use the `start` command:

```bat
set BRGLM2_PINNED=1
start /AFFINITY 1 /B /WAIT Rscript.exe benchmarks\run_all.R
```
*Note: If you try this in PowerShell, you must wrap it: `cmd.exe /c "..."`*

---


## Benchmark Components

### Sections in `main_benchmark.R`

| Section | Dataset | Method | Tool | Iterations |
|---------|---------|--------|------|------------|
| 1 | — | tinytest suite | — | — |
| 2 | Lizards (n=409) | brglmFit / binomial | `microbenchmark` | 500 |
| 3 | Endometrial (n=79) | brglmFit / probit | `microbenchmark` | 250 |
| 4 | MultipleFeatures (n=2000, p=432) | brglmFit | `bench::mark` | 30 |
| 5 | MultipleFeatures (n=2000, p=432) | mdyplFit | `bench::mark` | 30 |

All sections interleave `original` and `new` evaluations automatically. Sections 4 and 5 previously ran all original reps then all new reps, which was the main cause of inconsistent speedup estimates.

### Why Sections 4 & 5 Use `bench::mark` Instead of `microbenchmark`

`bench::mark` runs a gc() before each iteration and interleaves expressions by default, both of which `microbenchmark` does not do. For long-running fits (seconds each), this produces noticeably more stable medians. It also reports memory allocation and gc counts, which is useful for diagnosing regressions.

### Speedup Reporting

All speedup factors are reported with a 95% confidence interval derived from the Q25/Q75 ratio bounds, not just a single median ratio. A result like `1.8× (CI: 1.6–2.0×)` tells you the improvement is consistent; a wide CI like `1.8× (CI: 0.9–3.2×)` means you need more iterations.

---

## Adjusting Iteration Counts

Edit these values in `main_benchmark.R`:

```r
# Section 2 — Lizards
times = 500          # change for faster/slower runs

# Section 3 — Endometrial
times = 250

# Section 4 — MultipleFeatures brglmFit
iterations = 30      # bench::mark parameter

# Section 5 — MultipleFeatures mdyplFit
iterations = 30
```

Rough runtime guide:

| Config | Approx. runtime |
|--------|----------------|
| Default (as shipped) | ~10–15 min |
| Lizards 100, Endo 50, Large 10 | ~3–5 min |
| Lizards 1000, Endo 500, Large 50 | ~30–45 min |

---

## Outputs

Each run creates a timestamped directory under `benchmarks/results/`:

```
benchmarks/results/YYYYMMDD_HHMMSS/
├── benchmark_results.txt            # Full console output
├── benchmark_data.RData             # R objects (reload for re-plotting)
└── plots/
    ├── 01_lizards_comparison.png
    ├── 02_endometrial_comparison.png
    ├── 03_mf_brglmfit_comparison.png
    ├── 04_mf_mdyplfit_comparison.png
    ├── 05_speedup_summary.png
    └── 06_lizards_distribution.png
```

To regenerate plots from a previous run without re-running benchmarks:

```r
setwd("C:/path/to/brglm2")
source("benchmarks/setup.R")        # sets results_dir to new timestamp — override if needed:
# results_dir <- "benchmarks/results/20250101_120000"
source("benchmarks/visualize_results.R")
```

---

## Troubleshooting

| Error | Fix |
|-------|-----|
| `Original library path does not exist` | Update `LIB_ORIGINAL` in `setup.R` |
| `brglm2 not found in original library` | Run `R CMD INSTALL . --library=../benchmark_libs/original` from `brglm2-original/` |
| `__RELAUNCHED__` error in console | Normal — this is the CPU-pinning re-launch. The real run continues in the new process. |
| Plots fail | Check `benchmark_data.RData` exists; verify `ggplot2` is installed |
| Wide CI / high variance | Increase iteration counts; use the manual CPU-pinning command (Option B above) |

---

## Citations

- Kosmidis, I., & Firth, D. (2021). Jeffreys-prior penalty, finiteness and shrinkage in binomial-response generalized linear models. *Biometrika*, 108(1), 71–82.
- Sterzinger, P., & Kosmidis, I. (2024). An iteratively reweighted least squares algorithm for the maximum Diaconis-Ylvisaker prior penalized likelihood in binomial-response generalized linear models.

## License

Same as the parent brglm2 package.