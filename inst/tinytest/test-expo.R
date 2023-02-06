## source(system.file("inst", "brglm0/brglm0.R", package = "brglm2"))
data("lizards", package = "brglm2")

lfit <- glm(cbind(grahami, opalinus) ~ height + diameter +
                       light + time, family = binomial(), data=lizards,
                   method = "brglmFit", type = "ML", epsilon = 1e-10, maxit = 1000)

info <- get_information_function(lfit)

start_methods <- c("AS_mixed", "AS_mean", "AS_median", "correction", "MPL_Jeffreys", "ML")
keep <- c("coef", "se", "ci")

## ML
or_est_ml <- exp(coef(lfit))
or_ses_ml <- sqrt(diag(diag(or_est_ml) %*% vcov(lfit) %*% diag(or_est_ml)))
ml <- expo(update(lfit, type = "ML"), type = "ML")
for (met in start_methods) {
    obj <- expo(update(lfit, type = met), type = "ML")
    expect_identical(ml[keep], obj[keep])
}
expect_equal(coef(ml), or_est_ml)
expect_equal(ml$se, or_ses_ml, check.attributes = FALSE)
expect_equal(ml$ci, exp(confint(lfit)))

## Lyles et al (2012)
lfit_mix <- update(lfit, type = "AS_mixed")
or_est_bc1 <- exp(coef(lfit_mix)) * exp(-diag(vcov(lfit_mix))/2)
or_ses_bc1 <- sqrt(diag(diag(or_est_bc1) %*% solve(info(log(or_est_bc1))) %*% diag(or_est_bc1)))
for (met in start_methods) {
    bc1 <- expo(update(lfit_mix, type = met), type = "Lylesetal2012")
    expect_equal(coef(bc1), or_est_bc1)
    expect_equal(bc1$se, or_ses_bc1, check.attributes = FALSE)
}

## correction*
or_est_bc2 <- exp(coef(lfit_mix)) / (1 + diag(vcov(lfit_mix))/2)
or_ses_bc2 <- sqrt(diag(diag(or_est_bc2) %*% solve(info(log(or_est_bc2))) %*% diag(or_est_bc2)))
for (met in start_methods) {
    bc2 <- expo(update(lfit_mix, type = met), type = "correction*")
    expect_equal(coef(bc2), or_est_bc2)
    expect_equal(bc2$se, or_ses_bc2, check.attributes = FALSE)
}

## correction+
or_est_bc3 <- exp(coef(lfit_mix)) * (1 - diag(vcov(lfit_mix))/2)
or_ses_bc3 <- sqrt(diag(diag(or_est_bc3) %*% solve(info(log(or_est_bc3))) %*% diag(or_est_bc3)))
for (met in start_methods) {
    bc3 <- expo(update(lfit, type = met), type = "correction+")
    expect_equal(coef(bc3), or_est_bc3)
    expect_equal(bc3$se, or_ses_bc3, check.attributes = FALSE)
}

## AS median
lfit_med <- update(lfit, type = "AS_median")
or_est_med <- exp(coef(lfit_med))
or_ses_med <- sqrt(diag(diag(or_est_med) %*% solve(info(log(or_est_med))) %*% diag(or_est_med)))
for (met in start_methods) {
    med <- expo(update(lfit_med, type = met), type = "AS_median")
    expect_equal(coef(med), or_est_med)
    expect_equal(med$se, or_ses_med, check.attributes = FALSE)
}

## starting from a glm object
expo_methods <- c("ML", "correction*", "correction+", "Lylesetal2012", "AS_median")
lfit_glm <- glm(cbind(grahami, opalinus) ~ height + diameter +
                    light + time, family = binomial(), data=lizards)
for (met in expo_methods) {
    out_brglmFit <- expo(lfit_med, type = met)[keep]
    out_glm <- expo(lfit_glm, type = met)[keep]
    expect_equal(out_brglmFit, out_glm, tolerance = 1e-06)
}

expect_stdout(print(expo(lfit_glm)), "Odds ratios")

## ## Interpretation
## set.seed(123)
## dat <- data.frame(y = rexp(10), x = rnorm(10))
## mod <- glm(y ~ x, family = Gamma("log"), data = dat)
## expo_mod <- expo(mod, type = "ML")
## expect_stdout(print(expo_mod), "Multiplicative effects to the mean")

## set.seed(111)
## dat <- data.frame(y = exp(rnorm(10)), x = rnorm(10))
## mod <- glm(y ~ x, family = inverse.gaussian("log"), data = dat)
## expo_mod <- expo(mod, type = "ML")
## expect_stdout(print(expo_mod), "Multiplicative effects to the mean")

## set.seed(111)
## dat <- data.frame(x = rexp(10), y = rpois(10, 10))
## mod <- glm(y ~ x, family = poisson("log"), data = dat)
## expo_mod <- expo(mod, type = "ML")
## expect_stdout(print(expo_mod), "Multiplicative effects to the mean")

## set.seed(111)
## dat <- data.frame(x = rexp(100, 0.5))
## dat$y <- rbinom(100, 1, exp(-1 -0.2 * dat$x))
## mod <- glm(y ~ x, family = binomial("log"), data = dat)
## expo_mod <- expo(mod, type = "ML")
## expect_stdout(print(expo_mod), "Relative risks")


## library(parallel)
## nsimu <- 1000
## set.seed(123)
## Y <- simulate(lfit, nsimu)
## true_psi <- exp(coef(lfit))
## X <- model.matrix(lfit)

## methods <- c("ML", "correction+", "correction*", "Lylesetal2012", "AS_median")
## results <- as.list(numeric(length(methods)))
## names(results) <- methods
## for (met in methods) {
##     results[[met]] <- mclapply(1:nsimu, function(k) {
##         temp_data <- lizards
##         temp_data[c("grahami", "opalinus")] <- Y[[k]]
##         mod <- update(lfit, data = temp_data)
##         expo(mod, met)
##     }, mc.cores = 8)
## }

## get_bias <- function(res, truth) rowMeans(sapply(res, function(x) x$coef - truth))
## get_pu <- function(res, truth) rowMeans(sapply(res, function(x) x$coef < truth))
## get_mse <- function(res, truth) rowMeans(sapply(res, function(x) (x$coef - truth)^2))
## get_avgv <- function(res, truth) rowMeans(sapply(res, function(x) x$se^2))
## get_pu <- function(res, truth) rowMeans(sapply(res, function(x) x$coef < truth))
## get_coverage <- function(res, truth) {
##     rowMeans(sapply(res, function(x) {
##         (x$ci[, 1] < truth) & (x$ci[, 2] > truth)
##     }))
## }

## round(sapply(results, get_bias, truth = true_psi), 2)
## round(sapply(results, get_mse, truth = true_psi), 2)
## round(sapply(results, get_pu, truth = true_psi), 2)
## round(sapply(results, get_coverage, truth = true_psi), 2)
## vars <- sapply(results, get_mse, truth = true_psi) - sapply(results, get_bias, truth = true_psi)^2
## avgvars <- sapply(results, get_avgv, truth = true_psi)



