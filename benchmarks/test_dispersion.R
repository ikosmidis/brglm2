## ============================================================
## test_dispersion.R
##
## Diagnostic script to compare dispersion estimation between:
##   - brglmFit.R        (modified / trust-region version)
##   - brglmFit_original.R  (original step-halving version)
##   - glm()                (baseline MLE)
##
## Run this script from within your brglm2 package source tree,
## i.e. after sourcing/loading the package so that all helpers
## (compute_fit, estimate_dispersion, AS_*_adjustment, etc.) are
## available, then source this file.
##
## Usage (from R, with package loaded via devtools):
##   devtools::load_all(".")          # load helpers from current package
##   source("test_dispersion.R")
## ============================================================

library(brglm2)   # for data sets and brglmControl

## ---- helpers -------------------------------------------------------

fmt <- function(x, d = 6) formatC(x, digits = d, format = "f")

compare_dispersion <- function(label, fit_orig, fit_new, fit_glm = NULL) {
  cat("\n", strrep("=", 60), "\n")
  cat("DATASET / SCENARIO:", label, "\n")
  cat(strrep("-", 60), "\n")

  d_orig <- fit_orig$dispersion
  d_new  <- fit_new$dispersion

  cat(sprintf("  Original  dispersion : %s\n", fmt(d_orig)))
  cat(sprintf("  Modified  dispersion : %s\n", fmt(d_new)))
  if (!is.null(fit_glm)) {
    d_glm <- summary(fit_glm)$dispersion
    cat(sprintf("  glm()     dispersion : %s\n", fmt(d_glm)))
    cat(sprintf("  |Modified - glm|     : %s\n", fmt(abs(d_new - d_glm))))
    cat(sprintf("  |Original - glm|     : %s\n", fmt(abs(d_orig - d_glm))))
  }
  cat(sprintf("  |Modified - Original|: %s\n", fmt(abs(d_new - d_orig))))

  ## Beta comparison (should be close for ML type)
  b_orig <- coef(fit_orig)
  b_new  <- coef(fit_new)
  max_beta_diff <- max(abs(b_new - b_orig), na.rm = TRUE)
  cat(sprintf("  max |beta_mod - beta_orig| : %s\n", fmt(max_beta_diff)))

  ## Convergence flags
  cat(sprintf("  Original converged: %s  |  Modified converged: %s\n",
              fit_orig$converged, fit_new$converged))

  ## transformed_dispersion
  cat(sprintf("  Original  transformed_dispersion : %s\n",
              fmt(fit_orig$transformed_dispersion)))
  cat(sprintf("  Modified  transformed_dispersion : %s\n",
              fmt(fit_new$transformed_dispersion)))

  invisible(list(d_orig = d_orig, d_new = d_new))
}


## ---- source both fitters -------------------------------------------
## Temporarily rename so both can coexist in the same session.


cat("Both fitters loaded.\n")

## ---- wrappers for glm() method argument ----------------------------

fit_orig <- function(formula, family, data, type = "AS_mean", ...) {
  glm(formula, family = family, data = data,
      method = "brglmFit_original", type = type, ...)
}

fit_new <- function(formula, family, data, type = "AS_mean", ...) {
  glm(formula, family = family, data = data,
      method = "brglmFit", type = type, ...)
}


## ====================================================================
## SCENARIO 1 – Gaussian / identity: dispersion should equal
##              the residual variance (MLE or bias-reduced variant).
## ====================================================================
data("anorexia", package = "MASS")

cat("\n\n>>> SCENARIO 1: Gaussian (anorexia data)\n")
for (tp in c("ML", "AS_mean", "AS_median", "AS_mixed", "correction")) {
  fo <- Postwt ~ Prewt + Treat + offset(Prewt)
  g  <- glm(fo, family = gaussian, data = anorexia)
  o  <- fit_orig(fo, gaussian, anorexia, type = tp)
  n  <- fit_new (fo, gaussian, anorexia, type = tp)
  compare_dispersion(paste("Gaussian /", tp), o, n, if (tp == "ML") g else NULL)
}


## ====================================================================
## SCENARIO 2 – Gamma: continuous positive response, free dispersion.
##              This is one of the most sensitive cases.
## ====================================================================
data("coalition", package = "brglm2")

cat("\n\n>>> SCENARIO 2: Gamma (coalition data)\n")
for (tp in c("ML", "AS_mean", "AS_median", "AS_mixed", "correction")) {
  fo <- duration ~ fract + numst2
  g  <- glm(fo, family = Gamma, data = coalition)
  o  <- fit_orig(fo, Gamma, coalition, type = tp)
  n  <- fit_new (fo, Gamma, coalition, type = tp)
  compare_dispersion(paste("Gamma /", tp), o, n, if (tp == "ML") g else NULL)
}


## ====================================================================
## SCENARIO 3 – Inverse Gaussian
## ====================================================================
cat("\n\n>>> SCENARIO 3: Inverse Gaussian\n")
if (requireNamespace("statmod", quietly = TRUE)) {
  set.seed(42)
  n   <- 80
  dat_ig <- data.frame(
    x1 = rnorm(n), x2 = rbinom(n, 1, 0.5),
    y  = statmod::rinvgauss(n, mean = 2, dispersion = 0.5)
  )
  for (tp in c("ML", "AS_mean", "AS_median", "AS_mixed")) {
    fo <- y ~ x1 + x2
    g  <- glm(fo, family = inverse.gaussian, data = dat_ig)
    o  <- try(fit_orig(fo, inverse.gaussian, dat_ig, type = tp), silent = TRUE)
    nf <- try(fit_new (fo, inverse.gaussian, dat_ig, type = tp), silent = TRUE)
    if (!inherits(o, "try-error") && !inherits(nf, "try-error"))
      compare_dispersion(paste("InvGauss /", tp), o, nf, if (tp == "ML") g else NULL)
    else cat("  [SKIPPED – error in fit for type", tp, "]\n")
  }
} else {
  cat("[Skipping Scenario 3 – install statmod to enable]\n")
}


## ====================================================================
## SCENARIO 4 – Binomial: dispersion FIXED at 1, should stay 1.
## ====================================================================
data("endometrial", package = "brglm2")

## TODO check ML fitting and whether brglmFit_original is still correct under sep and this type

cat("\n\n>>> SCENARIO 4: Binomial (dispersion should == 1)\n")
for (tp in c("ML", "AS_mean", "AS_median", "AS_mixed")) {
  fo <- HG ~ NV + PI + EH
  g  <- glm(fo, family = binomial("probit"), data = endometrial)
  o  <- fit_orig(fo, binomial("probit"), endometrial, type = tp)
  nf <- fit_new (fo, binomial("probit"), endometrial, type = tp)
  compare_dispersion(paste("Binomial /", tp), o, nf, g)
}


## ====================================================================
## SCENARIO 5 – Poisson: dispersion FIXED at 1, should stay 1.
## ====================================================================
set.seed(7)
dat_pois <- data.frame(
  x = rnorm(50),
  y = rpois(50, lambda = exp(0.5 + 0.3 * rnorm(50)))
)
cat("\n\n>>> SCENARIO 5: Poisson (dispersion should == 1)\n")
for (tp in c("ML", "AS_mean", "AS_median")) {
  fo <- y ~ x
  g  <- glm(fo, family = poisson, data = dat_pois)
  o  <- fit_orig(fo, poisson, dat_pois, type = tp)
  nf <- fit_new (fo, poisson, dat_pois, type = tp)
  compare_dispersion(paste("Poisson /", tp), o, nf, g)
}


## ====================================================================
## SCENARIO 6 – Deep-dive: trace the dispersion step internals.
##   Runs one iteration manually with both fitters and prints
##   grad_zeta, inverse_info_zeta, and the resulting step.
##   Uses Gamma/coalition as the test case.
## ====================================================================
cat("\n\n", strrep("=", 60), "\n")
cat("SCENARIO 6: Internal dispersion step trace (Gamma / AS_mean)\n")
cat(strrep("=", 60), "\n")

cat("\n--- Original (trace=TRUE, maxit=3) ---\n")
tryCatch(
  glm(duration ~ fract + numst2, family = Gamma, data = coalition,
      method = "brglmFit_original", type = "AS_mean",
      control = list(maxit = 3, trace = TRUE)),
  error = function(e) cat("ERROR:", conditionMessage(e), "\n")
)

cat("\n--- Modified (trace=TRUE, maxit=3) ---\n")
tryCatch(
  glm(duration ~ fract + numst2, family = Gamma, data = coalition,
      method = "brglmFit", type = "AS_mean",
      control = list(maxit = 3, trace = TRUE)),
  error = function(e) cat("ERROR:", conditionMessage(e), "\n")
)


## ====================================================================
## SCENARIO 7 – Verify d1zeta scaling: check that the post-loop
##   transformation back to the requested scale is correct for
##   is_ML / is_AS_median / is_AS_mixed cases.
##   Extracts transformed_dispersion and info_transformed_dispersion.
## ====================================================================
cat("\n\n", strrep("=", 60), "\n")
cat("SCENARIO 7: Transformed dispersion scale check (Gamma)\n")
cat(strrep("=", 60), "\n")

for (tp in c("ML", "AS_mixed", "AS_median")) {
  o  <- fit_orig(duration ~ fract + numst2, Gamma, coalition, type = tp)
  nf <- fit_new (duration ~ fract + numst2, Gamma, coalition, type = tp)
  cat(sprintf("\n  type = %-12s  transformation = %s\n", tp,
              o$transformation))
  cat(sprintf("    Orig  transf_disp = %s   info_transf_disp = %s\n",
              fmt(o$transformed_dispersion), fmt(o$info_transformed_dispersion)))
  cat(sprintf("    New   transf_disp = %s   info_transf_disp = %s\n",
              fmt(nf$transformed_dispersion), fmt(nf$info_transformed_dispersion)))
  cat(sprintf("    diff  transf_disp = %s   diff info         = %s\n",
              fmt(abs(nf$transformed_dispersion - o$transformed_dispersion)),
              fmt(abs(nf$info_transformed_dispersion - o$info_transformed_dispersion))))
}


## ====================================================================
## SCENARIO 8 – vcov dispersion column: the final information matrix
##   entry used by vcov(..., model = "dispersion") must be consistent.
## ====================================================================
cat("\n\n", strrep("=", 60), "\n")
cat("SCENARIO 8: vcov() dispersion block comparison\n")
cat(strrep("=", 60), "\n")

for (tp in c("AS_mean", "AS_mixed", "AS_median")) {
  o  <- fit_orig(duration ~ fract + numst2, Gamma, coalition, type = tp)
  nf <- fit_new (duration ~ fract + numst2, Gamma, coalition, type = tp)
  vo <- tryCatch(vcov(o,  model = "dispersion"), error = function(e) NA)
  vn <- tryCatch(vcov(nf, model = "dispersion"), error = function(e) NA)
  cat(sprintf("  type = %s\n", tp))
  cat(sprintf("    vcov orig  = %s\n", fmt(vo)))
  cat(sprintf("    vcov new   = %s\n", fmt(vn)))
  cat(sprintf("    diff       = %s\n", fmt(abs(vn - vo))))
}

cat("\n\n>>> Diagnostic script complete.\n")