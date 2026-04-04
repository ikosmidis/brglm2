# =============================================================================
# setup.R
# Configuration, validation, and environment preparation for brglm2 benchmarks.
# Sourced by run_all.R before any benchmarking starts.
# =============================================================================

SEP <- strrep("=", 70)

cat(SEP, "\n")
cat("BRGLM2 BENCHMARK SETUP & VERIFICATION\n")
cat(SEP, "\n\n")

# ====== 1. Library Paths ======
# Edit these two lines to point at your compiled package copies.
# Relative paths work when the working directory is brglm2/ (recommended).
LIB_ORIGINAL <- "../benchmark_libs/original"
LIB_NEW      <- "../benchmark_libs/new"

# ====== 2. Required Packages ======
cat("Checking required packages...\n")

required_packages <- c("microbenchmark", "bench", "ggplot2", "tinytest")

missing_pkgs <- Filter(function(p) !requireNamespace(p, quietly = TRUE),
                       required_packages)

if (length(missing_pkgs) > 0) {
  cat("  Installing missing packages:", paste(missing_pkgs, collapse = ", "), "\n")
  install.packages(missing_pkgs)
}

suppressPackageStartupMessages({
  for (pkg in required_packages) library(pkg, character.only = TRUE)
})

cat("  OK - all packages loaded:", paste(required_packages, collapse = ", "), "\n\n")

# ====== 3. R Version Check ======
cat("Checking R version...\n")
cat(" ", R.version.string, "\n")
if (getRversion() < "4.0.0") warning("R >= 4.0.0 recommended.")
cat("  OK\n\n")

# ====== 4. Library Path Validation ======
cat("Checking library paths...\n")

check_lib <- function(path, label) {
  if (!dir.exists(path)) {
    stop(label, " library path does not exist: ", path,
         "\nPlease update setup.R.")
  }
  if (!dir.exists(file.path(path, "brglm2"))) {
    stop("brglm2 not found in ", label, " library: ", path,
         "\nRun: R CMD INSTALL . --library=", path)
  }
  cat("  OK -", label, "library:", path, "\n")
}

check_lib(LIB_ORIGINAL, "original")
check_lib(LIB_NEW,      "new")
cat("\n")

# ====== 5. Load Both brglm2 Versions ======
cat("Loading brglm2 versions...\n")

library(brglm2, lib.loc = LIB_ORIGINAL)
brglm2_original <- brglm2::brglmFit
mdypl_original  <- brglm2::mdyplFit
cat("  Original loaded from:", LIB_ORIGINAL, "\n")

detach("package:brglm2", unload = TRUE)

library(brglm2, lib.loc = LIB_NEW)
brglm2_new <- brglm2::brglmFit
mdypl_new  <- brglm2::mdyplFit
cat("  New loaded from:", LIB_NEW, "\n\n")

# ====== 6. Load Datasets ======
cat("Loading datasets...\n")

data("lizards",          package = "brglm2")
data("endometrial",      package = "brglm2")
data("MultipleFeatures", package = "brglm2")

full_fm  <- cbind(grahami, opalinus) ~ height + diameter + light + time
endo_fm  <- HG ~ NV + PI + EH

vars     <- grep("fou|kar", names(MultipleFeatures), value = TRUE)
train_id <- which(MultipleFeatures$training)
MultipleFeatures[train_id, vars] <- scale(MultipleFeatures[train_id, vars],
                                           scale = FALSE)
kappa      <- length(vars) / sum(MultipleFeatures$training)
full_mf_fm <- formula(paste("I(digit == 7) ~", paste(vars, collapse = " + ")))
alpha_val  <- 1 / (1 + kappa)

cat("  lizards           n =", nrow(lizards), "\n")
cat("  endometrial       n =", nrow(endometrial), "\n")
cat("  MultipleFeatures  n =", nrow(MultipleFeatures),
    " | features =", length(vars),
    " | kappa =", round(kappa, 4), "\n\n")

# ====== 7. Results Directory ======
timestamp   <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir <- file.path("benchmarks", "results", timestamp)
dir.create(file.path(results_dir, "plots"), recursive = TRUE)
cat("Results will be saved to:", results_dir, "\n\n")

# ====== 8. Windows CPU Affinity (Single Core) ======
# Pinning R to one logical core eliminates OS scheduler migration noise,
# which is the main source of run-to-run timing variance on Windows.
# This uses the Windows built-in `start` command — no extra software needed.
#
# HOW IT WORKS: We re-launch this script under a new Rscript process that is
# pinned to a single CPU core via the /AFFINITY flag of cmd.exe's `start`.
# Affinity 0x1 = core 0 only.  Change to 0x2 for core 1, 0x4 for core 2, etc.
#
# IMPORTANT: This block only fires when the script is run directly
# (i.e. not already pinned).  Set BRGLM2_PINNED=1 in your environment to skip.

pin_to_single_core_windows <- function() {
  if (.Platform$OS.type != "windows") return(invisible(NULL))
  if (identical(Sys.getenv("BRGLM2_PINNED"), "1")) {
    cat("CPU affinity: already pinned to single core.\n\n")
    return(invisible(NULL))
  }

  cat("CPU affinity: re-launching R pinned to core 0 for reproducibility...\n")
  cat("(Set BRGLM2_PINNED=1 in your environment to skip this step.)\n\n")

  # Build the command that will be relaunched
  rscript  <- file.path(R.home("bin"), "Rscript.exe")
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("--file=", args)]
  if (length(file_arg) == 0) {
    script <- ""
  } else {
    script <- normalizePath(sub("--file=", "", file_arg), mustWork = FALSE)
  }
  script   <- sub("--file=", "", script)

  if (!nzchar(script)) {
    cat("  NOTE: Cannot detect script path in interactive mode.\n")
    cat("  To pin manually, run from a terminal:\n\n")
    cat("    set BRGLM2_PINNED=1\n")
    cat("    start /AFFINITY 1 /B /WAIT Rscript.exe benchmarks\\run_all.R\n\n")
    cat("  Continuing unpinned in this session.\n\n")
    return(invisible(NULL))
  }

  cmd <- sprintf(
    'cmd.exe /C "set BRGLM2_PINNED=1 && start /AFFINITY 1 /B /WAIT "%s" "%s""',
    rscript, script
  )
  cat("  Running:", cmd, "\n\n")
  system(cmd)
  stop("__RELAUNCHED__")   # stop the current (unpinned) session cleanly
}

pin_to_single_core_windows()

# ====== 9. Utility Functions ======

# Print speedup with a 95% CI (conservative: Q25/Q75 ratio bounds)
report_speedup <- function(bench_obj, label) {
  orig <- bench_obj$time[bench_obj$expr == "original"]
  new  <- bench_obj$time[bench_obj$expr == "new"]

  med_speedup  <- median(orig) / median(new)
  lo_speedup   <- quantile(orig, 0.25) / quantile(new, 0.75)
  hi_speedup   <- quantile(orig, 0.75) / quantile(new, 0.25)

  cat(sprintf(
    "%s speedup: %.2fx  (95%% CI: %.2f, %.2fx)\n",
    label, med_speedup, lo_speedup, hi_speedup
  ))
  invisible(med_speedup)
}

# Numerical equality check between two fitted GLMs
compare_fits <- function(fit1, fit2, tolerance = 1e-10) {
  checks <- list(
    coefficients  = all.equal(coef(fit1),    coef(fit2),    tolerance = tolerance),
    fitted_values = all.equal(fitted(fit1),  fitted(fit2),  tolerance = tolerance),
    deviance      = all.equal(fit1$deviance, fit2$deviance, tolerance = tolerance)
  )
  cat("  Coefficients match: ",  isTRUE(checks$coefficients),  "\n")
  cat("  Fitted values match:",  isTRUE(checks$fitted_values), "\n")
  cat("  Deviance match:     ",  isTRUE(checks$deviance),      "\n")
  if (!isTRUE(checks$coefficients))
    cat("  Max coef diff:", max(abs(coef(fit1) - coef(fit2))), "\n")
  invisible(all(sapply(checks, isTRUE)))
}

cat(SEP, "\n")
cat("SETUP COMPLETE — ready to benchmark.\n")
cat(SEP, "\n\n")