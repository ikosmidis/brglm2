# =============================================================================
# run_all.R
# Master execution script for the brglm2 benchmark suite.
# Run this from the brglm2/ working directory:
#
#   setwd("C:/path/to/brglm2")
#   source("benchmarks/run_all.R")
#
# Or from a terminal (recommended for CPU pinning):
#   set BRGLM2_PINNED=1
#   start /AFFINITY 1 /B /WAIT Rscript.exe benchmarks\run_all.R
# =============================================================================

cat(strrep("=", 70), "\n")
cat("BRGLM2 COMPREHENSIVE BENCHMARK SUITE\n")
cat(strrep("=", 70), "\n\n")

cat("This will run:\n")
cat("  1. Setup verification\n")
cat("  2. tinytest suite validation\n")
cat("  3. Microbenchmarks: Lizards (500 iters) & Endometrial (250 iters)\n")
cat("  4. Large dataset: MultipleFeatures brglmFit  (30 iters, interleaved)\n")
cat("  5. Large dataset: MultipleFeatures mdyplFit  (30 iters, interleaved)\n")
cat("  6. Visualisation generation\n\n")
cat("Expected runtime: ~10-15 minutes\n\n")

response <- readline(prompt = "Continue? (y/n): ")
if (tolower(trimws(response)) != "y") {
  cat("Benchmark cancelled.\n")
  stop("Cancelled by user.", call. = FALSE)
}

# ---- Phase 1: Setup --------------------------------------------------------
cat("\n--- Phase 1: Setup & Verification ---\n\n")
source("benchmarks/setup.R")

# ---- Phase 2: Benchmarks ---------------------------------------------------
cat("\n--- Phase 2: Running Benchmarks ---\n\n")
suite_start <- Sys.time()
source("benchmarks/main_benchmark.R")
suite_end <- Sys.time()

# ---- Phase 3: Visualisations -----------------------------------------------
cat("\n--- Phase 3: Generating Visualisations ---\n\n")
source("benchmarks/visualize_results.R")

# ---- Summary ---------------------------------------------------------------
cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK SUITE COMPLETE\n")
cat(strrep("=", 70), "\n\n")

cat("Total benchmark runtime:", format(difftime(suite_end, suite_start)), "\n")
cat("Results saved to:       ", results_dir, "\n\n")

# Load the saved speedups and print a clean table
load(file.path(results_dir, "benchmark_data.RData"))

cat("Speedup summary:\n")
cat(sprintf("  %-35s %s\n", "Test", "Speedup"))
cat("  ", strrep("-", 45), "\n")
for (nm in names(speedups)) {
  cat(sprintf("  %-35s %.2fx\n", nm, speedups[[nm]]))
}
cat("  ", strrep("-", 45), "\n")
cat(sprintf("  %-35s %.2fx\n", "Mean", mean(unlist(speedups))))
cat("\n")