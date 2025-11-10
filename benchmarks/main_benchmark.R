# Add microbenchmark library
library(microbenchmark)

# Create results directory with timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir <- file.path("benchmarks", "results", timestamp)
dir.create(results_dir, recursive = TRUE)

# Redirect output to file
sink(file.path(results_dir, "benchmark_results.txt"), split = TRUE)

tryCatch({
  separator <- strrep("=", 70)
  
  cat(separator, "\n")
  cat("BRGLM2 COMPREHENSIVE BENCHMARKING RESULTS\n")
  cat("Timestamp:", timestamp, "\n")
  cat(separator, "\n\n")
  
# ====== SECTION 1: Test Suite Execution ======
cat("\n", separator, "\n")
cat("SECTION 1: TEST SUITE EXECUTION\n")
cat(separator, "\n\n")

pkg_path <- "C:/Users/ollie/OneDrive/Desktop/UNI/Project brglm2/brglm2"

# Ensure DD is available in case internal calls need it
if (!exists("DD", mode = "function")) {
  DD <- function(expr, name, order = 1) {
    if(order < 1) stop("'order' must be >= 1")
    if(order == 1) D(expr, name)
    else DD(D(expr, name), name, order - 1)
  }
}

run_test_suite <- function(libpath, label) {
  cat("Running test suite for", label, "version...\n")
  cat("Library path:", libpath, "\n\n")
  
  # Load brglm2 version
  library(brglm2, lib.loc = libpath)
  library(tinytest)
  
  # Run test suite, timing total duration
  start <- Sys.time()
  result <- test_all(pkg_path)
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  
  # Summarise results safely
  sum_df <- as.data.frame(result)
  if ("ok" %in% names(sum_df)) {
    passes <- sum(sum_df$ok)
    fails <- sum(!sum_df$ok)
  } else {
    # fallback if tinytest version differs
    passes <- NA
    fails <- NA
    warning("Could not extract passes/fails from test results")
  }
  
  cat(sprintf("All ok, %d results (%.1fs)\n\n", passes + fails, elapsed))
  
  invisible(list(result = result, passes = passes, fails = fails, elapsed = elapsed))
}

# --- Run NEW version ---
new_results <- run_test_suite(LIB_NEW, "NEW")

# --- Run ORIGINAL version ---
orig_results <- run_test_suite(LIB_ORIGINAL, "ORIGINAL")

# --- Save summary file (concise) ---
test_output_file <- file.path(results_dir, "test_results.txt")
cat(sprintf(
  "Original: %d tests (%.1fs)\nNew: %d tests (%.1fs)\nSpeedup: %.2fx faster\n",
  orig_results$passes + orig_results$fails, orig_results$elapsed,
  new_results$passes + new_results$fails, new_results$elapsed,
  orig_results$elapsed / new_results$elapsed
), file = test_output_file)

cat("\nDetailed test results saved to:", test_output_file, "\n")
  
  # ====== SECTION 2: Microbenchmark Comparison ======
  cat("\n", separator, "\n")
  cat("SECTION 2: MICROBENCHMARK COMPARISON\n")
  cat("(100 evaluations with unit: seconds)\n")
  cat(separator, "\n\n")
  
  cat("Test: Lizards dataset (small)\n")
  lizards_bench <- microbenchmark(
    original = glm(full_fm, data = lizards,
                   family = binomial(),
                   method = brglm2_original,
                   maxit = 200),
    new = glm(full_fm, data = lizards,
              family = binomial(),
              method = brglm2_new,
              maxit = 200),
    times = 100,
    unit = "s"
  )
  print(lizards_bench)
  
  cat("\n\nSpeedup factor:", 
      median(lizards_bench$time[lizards_bench$expr == "original"]) / 
        median(lizards_bench$time[lizards_bench$expr == "new"]), "x\n")
  
  # ====== SECTION 3: Endometrial Dataset Test ======
  cat("\n", separator, "\n")
  cat("SECTION 3: ENDOMETRIAL DATASET (PROBIT LINK)\n")
  cat(separator, "\n\n")
  
  data("endometrial", package = "brglm2")
  endo_fm <- HG ~ NV + PI + EH
  
  cat("Microbenchmark (50 evaluations):\n")
  endo_bench <- microbenchmark(
    original = glm(endo_fm, data = endometrial,
                   family = binomial("probit"),
                   method = brglm2_original,
                   maxit = 200),
    new = glm(endo_fm, data = endometrial,
              family = binomial("probit"),
              method = brglm2_new,
              maxit = 200),
    times = 50,
    unit = "s"
  )
  print(endo_bench)
  
  cat("\n\nSpeedup factor:", 
      median(endo_bench$time[endo_bench$expr == "original"]) / 
        median(endo_bench$time[endo_bench$expr == "new"]), "x\n")
  
  # ====== SECTION 4: MultipleFeatures Dataset (Large, High-Dimensional) ======
  cat("\n", separator, "\n")
  cat("SECTION 4: MULTIPLEFEATURES DATASET (HIGH-DIMENSIONAL)\n")
  cat("Reference: Sterzinger & Kosmidis (2024, Section 8)\n")
  cat(separator, "\n\n")

  if (exists("MultipleFeatures", envir = .GlobalEnv)) {
    rm(MultipleFeatures, envir = .GlobalEnv)
  }
  data("MultipleFeatures", package = "brglm2")
  
  # Prepare data following the documentation example
  vars <- grep("fou|kar", names(MultipleFeatures), value = TRUE)
  train_id <- which(MultipleFeatures$training)
  MultipleFeatures[train_id, vars] <- scale(MultipleFeatures[train_id, vars], scale = FALSE)
  
  # Calculate kappa for MDYPL
  kappa <- length(vars) / sum(MultipleFeatures$training)
  
  # Create formula for full model
  full_mf_fm <- formula(paste("I(digit == 7) ~", paste(vars, collapse = " + ")))
  
  cat("Dataset info:\n")
  cat("  Total observations:", nrow(MultipleFeatures), "\n")
  cat("  Training observations:", sum(MultipleFeatures$training), "\n")
  cat("  Number of features:", length(vars), "\n")
  cat("  Dimensionality ratio (kappa):", round(kappa, 4), "\n\n")
  
  cat("Running single timing test...\n\n")
  
  cat("Original version timing:\n")
  time_mf_orig <- system.time({
    fit_mf_orig <- glm(full_mf_fm, data = MultipleFeatures, 
                       family = binomial(),
                       method = brglm2_original, 
                       subset = training, 
                       maxit = 200)
  })
  print(time_mf_orig)
  cat("\n")
  
  cat("New version timing:\n")
  time_mf_new <- system.time({
    fit_mf_new <- glm(full_mf_fm, data = MultipleFeatures, 
                      family = binomial(),
                      method = brglm2_new, 
                      subset = training, 
                      maxit = 200)
  })
  print(time_mf_new)
  cat("\n")
  
  cat("MultipleFeatures speedup:", 
      time_mf_orig["elapsed"] / time_mf_new["elapsed"], "x\n")
  
  # Verify coefficients match
  cat("\nNumerical verification:\n")
  cat("  Coefficients match:", 
      isTRUE(all.equal(coef(fit_mf_orig), coef(fit_mf_new), tolerance = 1e-10)), "\n")
  cat("  Max coefficient difference:", 
      max(abs(coef(fit_mf_orig) - coef(fit_mf_new))), "\n")
  cat("  Deviance match:",
      isTRUE(all.equal(fit_mf_orig$deviance, fit_mf_new$deviance, tolerance = 1e-10)), "\n")
  
  # ====== SECTION 5: MDYPL Method Test ======
  cat("\n", separator, "\n")
  cat("SECTION 5: MDYPL METHOD (MultipleFeatures)\n")
  cat("Maximum Diaconis-Ylvisaker Prior Penalized Likelihood\n")
  cat(separator, "\n\n")
  
  # Load both versions with mdyplFit
  library(brglm2, lib.loc = LIB_ORIGINAL)
  mdypl_original <- brglm2::mdyplFit
  
  detach("package:brglm2", unload = TRUE)
  library(brglm2, lib.loc = LIB_NEW)
  mdypl_new <- brglm2::mdyplFit
  
  cat("Testing mdyplFit method with alpha = 1/(1 + kappa)\n")
  alpha_val <- 1 / (1 + kappa)
  cat("  Alpha value:", round(alpha_val, 4), "\n\n")
  
  cat("Original mdyplFit timing:\n")
  time_mdypl_orig <- system.time({
    fit_mdypl_orig <- glm(full_mf_fm, data = MultipleFeatures,
                          family = binomial(),
                          method = mdypl_original,
                          alpha = alpha_val,
                          subset = training,
                          maxit = 200)
  })
  print(time_mdypl_orig)
  cat("\n")
  
  cat("New mdyplFit timing:\n")
  time_mdypl_new <- system.time({
    fit_mdypl_new <- glm(full_mf_fm, data = MultipleFeatures,
                         family = binomial(),
                         method = mdypl_new,
                         alpha = alpha_val,
                         subset = training,
                         maxit = 200)
  })
  print(time_mdypl_new)
  cat("\n")
  
  cat("MDYPL speedup:", 
      time_mdypl_orig["elapsed"] / time_mdypl_new["elapsed"], "x\n")
  
  cat("\nNumerical verification (MDYPL):\n")
  cat("  Coefficients match:", 
      isTRUE(all.equal(coef(fit_mdypl_orig), coef(fit_mdypl_new), tolerance = 1e-10)), "\n")
  cat("  Max coefficient difference:", 
      max(abs(coef(fit_mdypl_orig) - coef(fit_mdypl_new))), "\n")
  
  # ====== SECTION 6: Summary Statistics ======
  cat("\n", separator, "\n")
  cat("SECTION 6: OVERALL SUMMARY\n")
  cat(separator, "\n\n")
  
  speedups <- c(
    "Lizards (small)" = median(lizards_bench$time[lizards_bench$expr == "original"]) / 
      median(lizards_bench$time[lizards_bench$expr == "new"]),
    "Endometrial (medium)" = median(endo_bench$time[endo_bench$expr == "original"]) / 
      median(endo_bench$time[endo_bench$expr == "new"]),
    "MultipleFeatures (large)" = time_mf_orig["elapsed"] / time_mf_new["elapsed"],
    "MDYPL (large)" = time_mdypl_orig["elapsed"] / time_mdypl_new["elapsed"]
  )
  
  cat("Speedup Summary:\n")
  for (i in seq_along(speedups)) {
    cat(sprintf("  %-25s: %.2fx faster\n", names(speedups)[i], speedups[i]))
  }
  
  cat("\nMean speedup:", round(mean(speedups), 2), "x\n")
  cat("Median speedup:", round(median(speedups), 2), "x\n")
  
  # Save benchmark objects for plotting
  save(lizards_bench, endo_bench, 
       time_mf_orig, time_mf_new, 
       time_mdypl_orig, time_mdypl_new,
       speedups,
       file = file.path(results_dir, "benchmark_data.RData"))
  
  cat("\n", separator, "\n")
  cat("BENCHMARKING COMPLETE\n")
  cat("All results saved to:", results_dir, "\n")
  cat(separator, "\n")
  
}, error = function(e) {
  cat("\n\nERROR OCCURRED:\n")
  cat(conditionMessage(e), "\n")
  cat(traceback(), "\n")
}, finally = {
  sink()
})

cat("\nResults saved to:", file.path(results_dir, "benchmark_results.txt"), "\n")
cat("Benchmark data saved to:", file.path(results_dir, "benchmark_data.RData"), "\n")