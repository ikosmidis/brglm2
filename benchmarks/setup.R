# Configuration for all benchmark scripts

# ====== Library Paths ======
LIB_ORIGINAL <- "C:/Users/ollie/OneDrive/Desktop/UNI/Project brglm2/benchmark_libs/original"
LIB_NEW <- "C:/Users/ollie/OneDrive/Desktop/UNI/Project brglm2/benchmark_libs/new"

# ====== Load Required Packages ======
cat("Loading benchmark packages...\n")

required_packages <- c("tictoc", "rbenchmark", "microbenchmark", "ggplot2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

cat("All packages loaded successfully.\n\n")

# ====== Load Both Versions of brglm2 ======
cat("Loading brglm2 versions...\n")

# Load original version
if (!dir.exists(LIB_ORIGINAL)) {
  stop("Original library path does not exist: ", LIB_ORIGINAL)
}

library(brglm2, lib.loc = LIB_ORIGINAL)
brglm2_original <- brglm2::brglmFit
cat("  Original brglmFit loaded from:", LIB_ORIGINAL, "\n")

# Detach and load new version
detach("package:brglm2", unload = TRUE)

if (!dir.exists(LIB_NEW)) {
  stop("New library path does not exist: ", LIB_NEW)
}

library(brglm2, lib.loc = LIB_NEW)
brglm2_new <- brglm2::brglmFit
cat("  New brglmFit loaded from:", LIB_NEW, "\n\n")

# ====== Verify Version Loading ======
if (!exists("brglm2_original") || !exists("brglm2_new")) {
  stop("Failed to load one or both brglm2 versions")
}

# ====== Setup Test Data ======
cat("Loading test datasets...\n")

# Lizards dataset (small)
data("lizards", package = "brglm2")
full_fm <- cbind(grahami, opalinus) ~ height + diameter + light + time
cat("  Lizards dataset loaded (n =", nrow(lizards), ")\n")

# Endometrial dataset (medium)
data("endometrial", package = "brglm2")
cat("  Endometrial dataset loaded (n =", nrow(endometrial), ")\n")

# MultipleFeatures dataset (large, high-dimensional)
if (exists("MultipleFeatures", envir = .GlobalEnv)) {
  rm(MultipleFeatures, envir = .GlobalEnv)
}
data("MultipleFeatures", package = "brglm2")
cat("  MultipleFeatures dataset loaded (n =", nrow(MultipleFeatures), ")\n")

# ====== Configuration Summary ======
cat("\n")
cat(strrep("=", 70), "\n")
cat("BENCHMARK SETUP COMPLETE\n")
cat(strrep("=", 70), "\n")
cat("Original library:", LIB_ORIGINAL, "\n")
cat("New library:", LIB_NEW, "\n")
cat("Datasets loaded: lizards, endometrial, MultipleFeatures\n")
cat("Packages loaded:", paste(required_packages, collapse = ", "), "\n")
cat(strrep("=", 70), "\n\n")

# ====== Utility Functions ======

# Function to safely run and time code
safe_time <- function(expr, name = "Expression") {
  cat("Timing:", name, "...\n")
  tryCatch({
    timing <- system.time(result <- expr)
    cat("  Completed in", timing["elapsed"], "seconds\n")
    list(result = result, timing = timing, success = TRUE)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    list(result = NULL, timing = NULL, success = FALSE, error = e)
  })
}

# Function to compare two fits for numerical equality
compare_fits <- function(fit1, fit2, tolerance = 1e-10, verbose = TRUE) {
  checks <- list(
    coefficients = all.equal(coef(fit1), coef(fit2), tolerance = tolerance),
    fitted_values = all.equal(fitted(fit1), fitted(fit2), tolerance = tolerance),
    deviance = all.equal(fit1$deviance, fit2$deviance, tolerance = tolerance)
  )
  
  if (verbose) {
    cat("\nNumerical Comparison:\n")
    cat("  Coefficients match:", isTRUE(checks$coefficients), "\n")
    if (!isTRUE(checks$coefficients)) {
      cat("    Max difference:", max(abs(coef(fit1) - coef(fit2))), "\n")
    }
    cat("  Fitted values match:", isTRUE(checks$fitted_values), "\n")
    if (!isTRUE(checks$fitted_values)) {
      cat("    Max difference:", max(abs(fitted(fit1) - fitted(fit2))), "\n")
    }
    cat("  Deviance match:", isTRUE(checks$deviance), "\n")
    if (!isTRUE(checks$deviance)) {
      cat("    Difference:", abs(fit1$deviance - fit2$deviance), "\n")
    }
  }
  
  all(sapply(checks, isTRUE))
}

# Function to format speedup nicely
format_speedup <- function(time_old, time_new) {
  speedup <- time_old / time_new
  improvement_pct <- (1 - 1/speedup) * 100
  
  list(
    speedup = speedup,
    improvement_pct = improvement_pct,
    text = sprintf("%.2fx faster (%.1f%% improvement)", speedup, improvement_pct)
  )
}

cat("Setup complete! Ready to run benchmarks.\n\n")