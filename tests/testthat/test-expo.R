context("exponentiating regression coefficients")

## source(system.file("inst", "brglm0/brglm0.R", package = "brglm2"))
data("lizards", package = "brglm2")

lfit <- glm(cbind(grahami, opalinus) ~ height + diameter +
                       light + time, family = binomial(), data=lizards,
                   method = "brglmFit", type = "ML", epsilon = 1e-10, maxit = 1000)


or_est <- exp(coef(lfit))
or_ses <- sqrt(diag(diag(or_est) %*% vcov(lfit) %*% diag(or_est)))

ml_from_ml <- expo(update(lfit, type = "ML"), type = "ML")
ml_from_meaBR <- expo(update(lfit, type = "AS_mean"), type = "ML")
ml_from_medBR <- expo(update(lfit, type = "AS_median"), type = "ML")
ml_from_corr <- expo(update(lfit, type = "correction"), type = "ML")
ml_from_jeffreys <- expo(update(lfit, type = "MPL_Jeffreys"), type = "ML")

keep <- c("coef", "se", "ci")

expect_identical(ml_from_ml[keep], ml_from_meaBR[keep])
expect_identical(ml_from_ml[keep], ml_from_medBR[keep])
expect_identical(ml_from_ml[keep], ml_from_corr[keep])
expect_identical(ml_from_ml[keep], ml_from_jeffreys[keep])
expect_equal(coef(ml_from_ml), or_est)
expect_equal(ml_from_ml$se, or_ses, check.attributes = FALSE)
expect_equal(ml_from_ml$ci, exp(confint(lfit)))










library(parallel)
nsimu <- 1000
set.seed(123)
Y <- simulate(lfit, nsimu)
true_psi <- exp(coef(lfit))
X <- model.matrix(lfit)

methods <- c("ML", "correction+", "correction*", "Lylesetal2012", "AS_median")
results <- as.list(numeric(length(methods)))
names(results) <- methods
for (met in methods) {
    results[[met]] <- mclapply(1:nsimu, function(k) {
        temp_data <- lizards
        temp_data[c("grahami", "opalinus")] <- Y[[k]]
        mod <- update(lfit, data = temp_data)
        expo(mod, met)
    }, mc.cores = 8)
}

get_bias <- function(res, truth) rowMeans(sapply(res, function(x) x$coef - truth))
get_pu <- function(res, truth) rowMeans(sapply(res, function(x) x$coef < truth))
get_mse <- function(res, truth) rowMeans(sapply(res, function(x) (x$coef - truth)^2))
get_avgv <- function(res, truth) rowMeans(sapply(res, function(x) x$se^2))
get_pu <- function(res, truth) rowMeans(sapply(res, function(x) x$coef < truth))
get_coverage <- function(res, truth) {
    rowMeans(sapply(res, function(x) {
        (x$ci[, 1] < truth) & (x$ci[, 2] > truth)
    }))
}

round(sapply(results, get_bias, truth = true_psi), 2)
round(sapply(results, get_mse, truth = true_psi), 2)
round(sapply(results, get_pu, truth = true_psi), 2)
round(sapply(results, get_coverage, truth = true_psi), 2)
vars <- sapply(results, get_mse, truth = true_psi) - sapply(results, get_bias, truth = true_psi)^2
avgvars <- sapply(results, get_avgv, truth = true_psi)
