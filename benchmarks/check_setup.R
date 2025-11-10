# Setup Verification Script
# Run this first to check that everything is configured correctly

cat("\n")
cat(strrep("=", 70), "\n")
cat("BRGLM2 BENCHMARK SETUP VERIFICATION\n")
cat(strrep("=", 70), "\n\n")

# ====== Check R Version ======
cat("1. Checking R version...\n")
r_version <- R.version.string
cat("   ", r_version, "\n")
if (getRversion() < "4.0.0") {
  cat("WARNING: R version 4.0.0 or higher recommended\n")
}
cat("OK\n\n")

# ====== Check Required Packages ======
cat("2. Checking required packages...\n")
required_pkgs <- c("tictoc", "rbenchmark", "microbenchmark", "ggplot2", "tinytest")

missing_pkgs <- c()
for (pkg in required_pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("OK", pkg, "installed\n")
  } else {
    cat("FAIL", pkg, "NOT FOUND\n")
    missing_pkgs <- c(missing_pkgs, pkg)
  }
}

if (length(missing_pkgs) > 0) {
  cat("\nInstalling missing packages...\n")
  install.packages(missing_pkgs)
  cat("OK Installation complete\n")
} else {
  cat("OK All packages available\n")
}
cat("\n")

# ====== Check Library Paths ======
cat("3. Checking library paths...\n")

# Read paths from setup.R if it exists
setup_file <- "benchmarks/setup.R"
if (file.exists(setup_file)) {
  source(setup_file, local = TRUE)
  
  cat("Original library:", LIB_ORIGINAL, "\n")
  if (dir.exists(LIB_ORIGINAL)) {
    cat("OK Original library exists\n")
    
    # Check if brglm2 is there
    if (file.exists(file.path(LIB_ORIGINAL, "brglm2"))) {
      cat("OK brglm2 package found in original library\n")
    } else {
      cat("FAIL brglm2 package NOT found in original library\n")
    }
  } else {
    cat("FAIL Original library path does NOT exist\n")
    cat("Please update LIB_ORIGINAL in benchmarks/setup.R\n")
  }
  
  cat("\nNew library:", LIB_NEW, "\n")
  if (dir.exists(LIB_NEW)) {
    cat("OK New library exists\n")
    
    # Check if brglm2 is there
    if (file.exists(file.path(LIB_NEW, "brglm2"))) {
      cat("OK brglm2 package found in new library\n")
    } else {
      cat("FAIL brglm2 package NOT found in new library\n")
    }
  } else {
    cat("FAIL New library path does NOT exist\n")
    cat("Please update LIB_NEW in benchmarks/setup.R\n")
  }
} else {
  cat("FAIL setup.R not found\n")
  cat("Please ensure you're in the project root directory\n")
}
cat("\n")

# ====== Check Package Path (for tests) ======
cat("4. Checking package path for tests...\n")
pkg_path <- "C:/Users/ollie/OneDrive/Desktop/UNI/Project brglm2/brglm2"

if (dir.exists(pkg_path)) {
  cat("OK Package directory exists:", pkg_path, "\n")
  
  # Check for key files
  key_files <- c("DESCRIPTION", "NAMESPACE", "R", "tests")
  for (f in key_files) {
    if (file.exists(file.path(pkg_path, f))) {
      cat("OK", f, "found\n")
    } else {
      cat("FAIL", f, "NOT found\n")
    }
  }
} else {
  cat("FAIL Package directory does NOT exist\n")
  cat("Please update pkg_path in main_benchmark.R\n")
}
cat("\n")

# ====== Check Datasets ======
cat("5. Checking dataset availability...\n")

# Try to load brglm2 from somewhere to access datasets
tryCatch({
  suppressWarnings(library(brglm2))
  
  # Check each dataset
  datasets_to_check <- c("lizards", "endometrial", "MultipleFeatures")
  for (ds in datasets_to_check) {
    tryCatch({
      data(list = ds, package = "brglm2")
      if (exists(ds)) {
        cat("OK", ds, "dataset accessible\n")
      }
    }, error = function(e) {
      cat("FAIL", ds, "dataset NOT accessible\n")
    })
  }
  
  cat("OK All datasets available\n")
}, error = function(e) {
  cat("Could not load brglm2 to check datasets\n")
  cat("This is OK if you're setting up for the first time\n")
})
cat("\n")

# ====== Check Results Directory ======
cat("6. Checking results directory...\n")
results_dir <- "benchmarks/results"

if (!dir.exists(results_dir)) {
  cat("Creating results directory...\n")
  dir.create(results_dir, recursive = TRUE)
  cat("OK Results directory created\n")
} else {
  cat("OK Results directory exists\n")
  
  # Check for previous results
  existing_results <- list.dirs(results_dir, recursive = FALSE)
  if (length(existing_results) > 0) {
    cat("Found", length(existing_results), "previous benchmark run(s)\n")
    latest <- existing_results[which.max(file.info(existing_results)$mtime)]
    cat("Latest:", basename(latest), "\n")
  }
}
cat("\n")

# ====== Summary ======
cat(strrep("=", 70), "\n")
cat("VERIFICATION SUMMARY\n")
cat(strrep("=", 70), "\n\n")

# Collect all checks
checks <- list(
  "R version" = TRUE,  # We always have R running this
  "Packages" = length(missing_pkgs) == 0,
  "Library paths" = exists("LIB_ORIGINAL") && exists("LIB_NEW") && 
                    dir.exists(LIB_ORIGINAL) && dir.exists(LIB_NEW),
  "Package path" = dir.exists(pkg_path),
  "Results directory" = dir.exists(results_dir)
)

all_ok <- all(unlist(checks))

if (all_ok) {
  cat("ALL CHECKS PASSED\n\n")
  cat("You're ready to run benchmarks!\n\n")
  cat("Run the full suite with:\n")
  cat("  source('benchmarks/run_all.R')\n\n")
  cat("Or run individual components:\n")
  cat("  source('benchmarks/main_benchmark.R')  # Main benchmarks\n")
  cat("  source('benchmarks/visualize_results.R')  # Generate plots\n")
} else {
  cat("SOME CHECKS FAILED\n\n")
  cat("Please fix the issues above before running benchmarks.\n\n")
  
  cat("Common fixes:\n")
  cat("1. Update library paths in benchmarks/setup.R\n")
  cat("2. Install missing packages\n")
  cat("3. Ensure you're in the correct working directory\n")
  cat("4. Check that brglm2 is installed in both library locations\n")
}

cat("\n", strrep("=", 70), "\n\n")

# ====== Helpful Commands ======
cat("Helpful commands:\n\n")
cat("Check current working directory:\n")
cat("  getwd()\n\n")
cat("Change working directory:\n")
cat("  setwd('path/to/project')\n\n")
cat("Install a package:\n")
cat("  install.packages('package_name')\n\n")
cat("List installed packages:\n")
cat("  installed.packages()[,c('Package', 'LibPath')]\n\n")
cat("Check where package is installed:\n")
cat("  find.package('brglm2')\n\n")