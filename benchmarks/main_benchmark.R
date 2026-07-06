# =============================================================================
# main_benchmark.R
# Core benchmarking tests.  Sourced by run_all.R after setup.R.
#
# Key reproducibility improvements over the previous version:
#   - Large dataset tests now interleave original/new (via bench::mark) instead
#     of running all original first then all new, which biased results.
#   - gc() + Sys.sleep() between sections to drain GC pressure and let CPU
#     thermals settle.
#   - Speedup factors reported with 95% confidence intervals, not just medians.
# =============================================================================

SEP <- strrep("=", 70)

sink(file.path(results_dir, "benchmark_results.txt"), split = TRUE)

tryCatch({

  cat(SEP, "\n")
  cat("BRGLM2 COMPREHENSIVE BENCHMARKING RESULTS\n")
  cat("Timestamp:", timestamp, "\n")
  cat(SEP, "\n\n")

  # ===========================================================================
  # SECTION 1: Test Suite
  # ===========================================================================
  cat("\n", SEP, "\n")
  cat("SECTION 1: TEST SUITE (NEW IMPLEMENTATION)\n")
  cat(SEP, "\n\n")

  library(tinytest)
  cat("Running tinytest suite against currently loaded brglm2 (new version)...\n")
  start_ts <- Sys.time()

  sink(NULL)                                        # suppress test noise briefly
  test_result <- tinytest::test_package("brglm2")
  sink(file.path(results_dir, "benchmark_results.txt"), append = TRUE, split = TRUE)

  elapsed_ts <- as.numeric(difftime(Sys.time(), start_ts, units = "secs"))
  sum_df     <- as.data.frame(test_result)
  n_passed   <- sum(sum_df$result == TRUE)
  n_total    <- nrow(sum_df)

  cat(sprintf("Result: %d/%d tests passed  (%.2fs)\n\n", n_passed, n_total, elapsed_ts))

  gc(full = TRUE); Sys.sleep(2)

  # ===========================================================================
  # SECTION 2: Lizards — small dataset (microbenchmark, interleaved)
  # ===========================================================================
  cat("\n", SEP, "\n")
  cat("SECTION 2: LIZARDS DATASET (SMALL, BINOMIAL)\n")
  cat(SEP, "\n\n")

  lizards_bench <- microbenchmark(
    original = glm(full_fm, data = lizards,
                   family = binomial(), method = brglm2_original, maxit = 200),
    new      = glm(full_fm, data = lizards,
                   family = binomial(), method = brglm2_new,      maxit = 200),
    times = 5000,    # interleaved by microbenchmark automatically
    unit  = "s"
  )
  print(lizards_bench)
  cat("\n")
  lizards_speedup <- report_speedup(lizards_bench, "Lizards")

  gc(full = TRUE); Sys.sleep(2)

  # ===========================================================================
  # SECTION 3: Endometrial — medium dataset (microbenchmark, interleaved)
  # ===========================================================================
  cat("\n", SEP, "\n")
  cat("SECTION 3: ENDOMETRIAL DATASET (MEDIUM, PROBIT)\n")
  cat(SEP, "\n\n")

  endo_bench <- microbenchmark(
    original = glm(endo_fm, data = endometrial,
                   family = binomial("probit"), method = brglm2_original, maxit = 200),
    new      = glm(endo_fm, data = endometrial,
                   family = binomial("probit"), method = brglm2_new,      maxit = 200),
    times = 2500,
    unit  = "s"
  )
  print(endo_bench)
  cat("\n")
  endo_speedup <- report_speedup(endo_bench, "Endometrial")

  gc(full = TRUE); Sys.sleep(2)

  # ===========================================================================
  # SECTION 4: MultipleFeatures brglmFit — large dataset (bench::mark)
  #
  # bench::mark interleaves iterations of both expressions automatically and
  # forces a gc() before each iteration, giving a much fairer comparison than
  # two separate replicate() blocks.
  # ===========================================================================
  cat("\n", SEP, "\n")
  cat("SECTION 4: MULTIPLEFEATURES DATASET — brglmFit (HIGH-DIMENSIONAL)\n")
  cat("Reference: Sterzinger & Kosmidis (2024, Section 8)\n")
  cat(SEP, "\n\n")

  cat("Dataset info:\n")
  cat("  Total observations:      ", nrow(MultipleFeatures), "\n")
  cat("  Training observations:   ", sum(MultipleFeatures$training), "\n")
  cat("  Features:                ", length(vars), "\n")
  cat("  Dimensionality ratio:    ", round(kappa, 4), "\n\n")

  # Warm-up fits (not timed) — populate caches so first timed iteration isn't
  # penalised by cold-cache effects.
  invisible(glm(full_mf_fm, data = MultipleFeatures, family = binomial(),
                method = brglm2_original, subset = training, maxit = 200))
  invisible(glm(full_mf_fm, data = MultipleFeatures, family = binomial(),
                method = brglm2_new,      subset = training, maxit = 200))
  gc(full = TRUE); Sys.sleep(1)

  mf_bench <- bench::mark(
    original = glm(full_mf_fm, data = MultipleFeatures, family = binomial(),
                   method = brglm2_original, subset = training, maxit = 200),
    new      = glm(full_mf_fm, data = MultipleFeatures, family = binomial(),
                   method = brglm2_new,      subset = training, maxit = 200),
    iterations = 150,   # bench::mark interleaves; 30 each is enough here
    check      = FALSE
  )
  print(mf_bench[, c("expression", "min", "median", "itr/sec", "mem_alloc", "n_gc")])

  orig_mf_med <- median(mf_bench$time[[1]], na.rm = TRUE)
  new_mf_med  <- median(mf_bench$time[[2]], na.rm = TRUE)
  orig_mf_q25 <- quantile(mf_bench$time[[1]], 0.25, na.rm = TRUE)
  new_mf_q75  <- quantile(mf_bench$time[[2]], 0.75, na.rm = TRUE)
  orig_mf_q75 <- quantile(mf_bench$time[[1]], 0.75, na.rm = TRUE)
  new_mf_q25  <- quantile(mf_bench$time[[2]], 0.25, na.rm = TRUE)

  cat(sprintf(
    "\nMultipleFeatures (brglmFit) speedup: %.2fx  (95%% CI: %.2f, %.2fx)\n",
    as.numeric(orig_mf_med / new_mf_med),
    as.numeric(orig_mf_q25 / new_mf_q75),
    as.numeric(orig_mf_q75 / new_mf_q25)
  ))

  # Numerical verification
  cat("\nNumerical verification:\n")
  fit_mf_orig <- glm(full_mf_fm, data = MultipleFeatures, family = binomial(),
                     method = brglm2_original, subset = training, maxit = 200)
  fit_mf_new  <- glm(full_mf_fm, data = MultipleFeatures, family = binomial(),
                     method = brglm2_new,      subset = training, maxit = 200)
  compare_fits(fit_mf_orig, fit_mf_new)

  gc(full = TRUE); Sys.sleep(2)

  # ===========================================================================
  # SECTION 5: MultipleFeatures mdyplFit — large dataset (bench::mark)
  # ===========================================================================
  cat("\n", SEP, "\n")
  cat("SECTION 5: MULTIPLEFEATURES DATASET — mdyplFit (MDYPL METHOD)\n")
  cat("Maximum Diaconis-Ylvisaker Prior Penalised Likelihood\n")
  cat(SEP, "\n\n")

  cat("  Alpha value:", round(alpha_val, 4), "\n\n")

  # Warm-up
  invisible(glm(full_mf_fm, data = MultipleFeatures, family = binomial(),
                method = mdypl_original, alpha = alpha_val,
                subset = training, maxit = 200))
  invisible(glm(full_mf_fm, data = MultipleFeatures, family = binomial(),
                method = mdypl_new, alpha = alpha_val,
                subset = training, maxit = 200))
  gc(full = TRUE); Sys.sleep(1)

  mdypl_bench <- bench::mark(
    original = glm(full_mf_fm, data = MultipleFeatures, family = binomial(),
                   method = mdypl_original, alpha = alpha_val,
                   subset = training, maxit = 200),
    new      = glm(full_mf_fm, data = MultipleFeatures, family = binomial(),
                   method = mdypl_new,      alpha = alpha_val,
                   subset = training, maxit = 200),
    iterations = 500,
    check      = FALSE
  )
  print(mdypl_bench[, c("expression", "min", "median", "itr/sec", "mem_alloc", "n_gc")])

  orig_md_med <- median(mdypl_bench$time[[1]], na.rm = TRUE)
  new_md_med  <- median(mdypl_bench$time[[2]], na.rm = TRUE)
  orig_md_q25 <- quantile(mdypl_bench$time[[1]], 0.25, na.rm = TRUE)
  new_md_q75  <- quantile(mdypl_bench$time[[2]], 0.75, na.rm = TRUE)
  orig_md_q75 <- quantile(mdypl_bench$time[[1]], 0.75, na.rm = TRUE)
  new_md_q25  <- quantile(mdypl_bench$time[[2]], 0.25, na.rm = TRUE)

  cat(sprintf(
    "\nMultipleFeatures (mdyplFit) speedup: %.2fx  (95%% CI: %.2f, %.2fx)\n",
    as.numeric(orig_md_med / new_md_med),
    as.numeric(orig_md_q25 / new_md_q75),
    as.numeric(orig_md_q75 / new_md_q25)
  ))

  # Numerical verification
  cat("\nNumerical verification:\n")
  fit_mdypl_orig <- glm(full_mf_fm, data = MultipleFeatures, family = binomial(),
                        method = mdypl_original, alpha = alpha_val,
                        subset = training, maxit = 200)
  fit_mdypl_new  <- glm(full_mf_fm, data = MultipleFeatures, family = binomial(),
                        method = mdypl_new,      alpha = alpha_val,
                        subset = training, maxit = 200)
  compare_fits(fit_mdypl_orig, fit_mdypl_new)

  # ===========================================================================
  # Save objects needed by visualize_results.R
  # ===========================================================================
  speedups <- c(
    "Lizards (small)"          = lizards_speedup,
    "Endometrial (medium)"     = endo_speedup,
    "MultipleFeatures brglmFit" = as.numeric(orig_mf_med  / new_mf_med),
    "MultipleFeatures mdyplFit" = as.numeric(orig_md_med  / new_md_med)
  )

  save(lizards_bench, endo_bench, mf_bench, mdypl_bench, speedups,
       fit_mf_orig, fit_mf_new, fit_mdypl_orig, fit_mdypl_new,
       file = file.path(results_dir, "benchmark_data.RData"))

  cat("\n", SEP, "\n")
  cat("BENCHMARKING COMPLETE\n")
  cat("Results saved to:", results_dir, "\n")
  cat(SEP, "\n")

}, error = function(e) {
  cat("\n\nERROR:\n", conditionMessage(e), "\n")
}, finally = {
  sink()
})