# Master script to run all benchmarks
# This script runs the complete benchmarking suite and generates visualizations

cat("======================================================================\n")
cat("BRGLM2 COMPREHENSIVE BENCHMARK SUITE\n")
cat("======================================================================\n\n")

cat("This will run:\n")
cat("  1. Test suite validation\n")
cat("  2. Microbenchmark comparisons (small & medium datasets)\n")
cat("  3. Large dataset tests (MultipleFeatures)\n")
cat("  4. MDYPL method comparison\n")
cat("  5. Visualization generation\n\n")

cat("Expected runtime: ~ 5 minutes depending on system\n\n")

response <- readline(prompt = "Continue? (y/n): ")

if (tolower(response) != "y") {
  cat("Benchmark cancelled.\n")
  quit(save = "no")
}

cat("\n======================================================================\n")
cat("PHASE 1: Running validation and setup\n")
cat("======================================================================\n\n")

# Verify and Setup
source("benchmarks/check_setup.R")
source("benchmarks/setup.R")

cat("\n======================================================================\n")
cat("PHASE 2: Running main benchmark script\n")
cat("======================================================================\n\n")

# Run main benchmark
start_time <- Sys.time()
source("benchmarks/main_benchmark.R")
end_time <- Sys.time()

cat("\n======================================================================\n")
cat("PHASE 3: Generating visualizations\n")
cat("======================================================================\n\n")

# Generate plots
source("benchmarks/visualize_results.R")

cat("\n======================================================================\n")
cat("BENCHMARK SUITE COMPLETE\n")
cat("======================================================================\n\n")

cat("Total runtime:", format(difftime(end_time, start_time)), "\n")
cat("Results directory:", file.path("benchmarks", "results"), "\n\n")

cat("Summary of benchmarks:\n")
results_dir <- file.path("benchmarks", "results")
dirs <- list.dirs(results_dir, recursive = FALSE)
latest_dir <- dirs[which.max(file.info(dirs)$mtime)]
load(file.path(latest_dir, "benchmark_data.RData"))

speedups <- c("Test suite" = orig_results$elapsed / new_results$elapsed, speedups)

for (i in seq_along(speedups)) {
  cat(sprintf("  %-30s: %.2fx faster\n", names(speedups)[i], speedups[i]))
}
cat(sprintf("  %-30s: %.2fx faster\n", "Average", mean(speedups)))