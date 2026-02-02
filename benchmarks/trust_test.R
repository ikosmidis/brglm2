# Minimal Test - Quick validation of trust region implementation

library(brglm2)

# Simple test: lizards dataset
data("lizards", package = "brglm2")

cat("Running simple test on lizards dataset...\n\n")

# Quasi-Fisher (reference implementation)
fit_qf <- glm(cbind(grahami, opalinus) ~ height + diameter + light + time, 
              family = binomial(logit), 
              data = lizards,
              method = "brglmFit")

# Trust Region (new implementation)
fit_tr <- glm(cbind(grahami, opalinus) ~ height + diameter + light + time, 
              family = binomial(logit), 
              data = lizards,
              method = "brglmFit_trustregion")

# Compare coefficients
cat("Coefficients comparison:\n")
comparison <- data.frame(
    QuasiFisher = coef(fit_qf),
    TrustRegion = coef(fit_tr),
    Difference = coef(fit_tr) - coef(fit_qf)
)
print(comparison)

cat("\n")
cat("Max absolute difference:", max(abs(coef(fit_tr) - coef(fit_qf))), "\n")
cat("Quasi-Fisher iterations:", fit_qf$iter, "\n")
cat("Trust Region iterations:", fit_tr$iter, "\n")
cat("Both converged:", fit_qf$converged & fit_tr$converged, "\n")

# Success criterion
if (max(abs(coef(fit_tr) - coef(fit_qf))) < 1e-4) {
    cat("\nTEST PASSED: Coefficients agree to 4 decimal places\n")
} else {
    cat("\nTEST FAILED: Coefficients differ by more than tolerance\n")
}