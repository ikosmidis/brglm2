# Load required packages
library(brglm2)
library(enrichwith)  # Needed for enriching family objects
library(numDeriv)    # Needed for numerical derivatives
library(profvis)
library(microbenchmark)

# Source the helper file to get compute_fit
# Try different paths depending on where script is run from
if (file.exists("R/brglmFit-helper.R")) {
  source("R/brglmFit-helper.R")
} else if (file.exists("../R/brglmFit-helper.R")) {
  source("../R/brglmFit-helper.R")
} else {
  stop("Cannot find brglmFit-helper.R. Please check your working directory.")
}

# Verify compute_fit is available
if (!exists("compute_fit")) {
  stop("compute_fit function not found after sourcing. Check the helper file.")
}

cat("compute_fit function loaded successfully\n")

# Load and prepare the dataset (same as your existing script)
data("MultipleFeatures", package = "brglm2")

vars <- grep("fou|kar", names(MultipleFeatures), value = TRUE)
training <- which(MultipleFeatures$training)
MultipleFeatures[training, vars] <- scale(MultipleFeatures[training, vars], 
                                           scale = FALSE)

# Create the formula and prepare data
full_mf_fm <- formula(paste("I(digit == 7) ~", paste(vars, collapse = " + ")))

# Prepare inputs for compute_fit (mimicking what brglmFit does)
train_data <- MultipleFeatures[training, ]
y <- as.numeric(train_data$digit == 7)
x <- model.matrix(full_mf_fm, data = train_data)
weights <- rep(1, length(y))
offset <- rep(0, length(y))

# Use enrichment approach (matching brglmFit.R line 426)
family <- enrichwith::enrich(binomial(), with = c("d1afun", "d2afun", "d3afun", "d1variance"))

# Add d2mu.deta if not already present
if (is.null(family$d2mu.deta)) {
  mu.eta <- family$mu.eta
  family$d2mu.deta <- function(eta) {
    numDeriv::grad(mu.eta, eta)
  }
}

# Get dimensions
nobs <- length(y)
nvars <- ncol(x)
keep <- weights > 0

# Get starting values (use a quick ML fit)
ml_fit <- glm(full_mf_fm, 
              data = train_data,
              family = binomial())

# Starting parameter vector
pars <- c(coef(ml_fit), 1)  # betas + dispersion

# ===== FOCUSED PROFILING OF compute_fit =====

cat("\n=== Profiling compute_fit function ===\n")

# First, test that compute_fit works outside profvis
cat("\nTesting compute_fit directly...\n")
test_result <- try({
  compute_fit(
    pars = pars,
    y = y,
    x = x,
    weights = weights,
    offset = offset,
    family = family,
    fixed_totals = NULL,
    row_totals = NULL,
    no_dispersion = FALSE,
    nobs = nobs,
    nvars = nvars,
    keep = keep,
    need_qr = TRUE,
    need_hatvalues = TRUE
  )
}, silent = FALSE)

if (inherits(test_result, "try-error")) {
  cat("\nERROR: compute_fit failed. Error details:\n")
  print(test_result)
  stop("Cannot proceed with profiling - compute_fit is not working")
} else {
  cat("Success! compute_fit executed without errors.\n")
  cat("Result class:", class(test_result), "\n")
}

# 1. Detailed line-by-line profiling with profvis
cat("\n1. Visual profiling with profvis (check browser/viewer)...\n")
p <- profvis({
  for (i in 1:100) {  # Run multiple times to get good sampling
    fit_result <- compute_fit(
      pars = pars,
      y = y,
      x = x,
      weights = weights,
      offset = offset,
      family = family,
      fixed_totals = NULL,
      row_totals = NULL,
      no_dispersion = FALSE,
      nobs = nobs,
      nvars = nvars,
      keep = keep,
      need_qr = TRUE,
      need_hatvalues = TRUE
    )
  }
})
print(p)

htmlwidgets::saveWidget(p, "fit_profile(03.12.25).html")

# 2. Benchmark different scenarios
cat("\n2. Benchmarking compute_fit under different conditions...\n")

# Test with QR and hatvalues
bench_full <- microbenchmark(
  with_qr_and_hat = compute_fit(
    pars = pars, y = y, x = x, weights = weights, offset = offset,
    family = family, fixed_totals = NULL, row_totals = NULL,
    no_dispersion = FALSE, nobs = nobs, nvars = nvars, keep = keep,
    need_qr = TRUE, need_hatvalues = TRUE
  ),
  with_qr_no_hat = compute_fit(
    pars = pars, y = y, x = x, weights = weights, offset = offset,
    family = family, fixed_totals = NULL, row_totals = NULL,
    no_dispersion = FALSE, nobs = nobs, nvars = nvars, keep = keep,
    need_qr = TRUE, need_hatvalues = FALSE
  ),
  no_qr = compute_fit(
    pars = pars, y = y, x = x, weights = weights, offset = offset,
    family = family, fixed_totals = NULL, row_totals = NULL,
    no_dispersion = FALSE, nobs = nobs, nvars = nvars, keep = keep,
    need_qr = FALSE, need_hatvalues = FALSE
  ),
  times = 50
)
print(bench_full)

# 3. Detailed Rprof for function-level timing
cat("\n3. Function-level profiling with Rprof...\n")
Rprof("compute_fit_profile.out", interval = 0.001, line.profiling = TRUE)
for (i in 1:200) {
  fit_result <- compute_fit(
    pars = pars, y = y, x = x, weights = weights, offset = offset,
    family = family, fixed_totals = NULL, row_totals = NULL,
    no_dispersion = FALSE, nobs = nobs, nvars = nvars, keep = keep,
    need_qr = TRUE, need_hatvalues = TRUE
  )
}
Rprof(NULL)

prof_summary <- summaryRprof("compute_fit_profile.out", lines = "both")
cat("\nTop functions by total time:\n")
print(head(prof_summary$by.total, 20))
cat("\nTop functions by self time:\n")
print(head(prof_summary$by.self, 20))

if (!is.null(prof_summary$by.line)) {
  cat("\nTop lines by total time:\n")
  print(head(prof_summary$by.line$line.out, 30))
}

# 4. Memory profiling
cat("\n4. Memory profiling...\n")
Rprof("compute_fit_memory.out", memory.profiling = TRUE)
for (i in 1:100) {
  fit_result <- compute_fit(
    pars = pars, y = y, x = x, weights = weights, offset = offset,
    family = family, fixed_totals = NULL, row_totals = NULL,
    no_dispersion = FALSE, nobs = nobs, nvars = nvars, keep = keep,
    need_qr = TRUE, need_hatvalues = TRUE
  )
}
Rprof(NULL)

mem_summary <- summaryRprof("compute_fit_memory.out", memory = "both")
cat("\nMemory usage by function:\n")
print(head(mem_summary$by.total, 15))

cat("\n=== Profiling complete ===\n")
cat("Profile output files created:\n")
cat("  - compute_fit_profile.out (detailed timing)\n")
cat("  - compute_fit_memory.out (memory usage)\n")