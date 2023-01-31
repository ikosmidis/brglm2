context("exponentiating regression coefficients")

## source(system.file("inst", "brglm0/brglm0.R", package = "brglm2"))
data("lizards", package = "brglm2")

lizardsML <- glm(cbind(grahami, opalinus) ~ height + diameter +
                     light + time, family = binomial(), data=lizards,
                 method = "brglmFit", type = "ML", epsilon = 1e-10, maxit = 1000)

lizardsBR <- glm(cbind(grahami, opalinus) ~ height + diameter +
                     light + time, family = binomial(), data=lizards,
                 method = "brglmFit", type = "AS_mean", epsilon = 1e-10, maxit = 1000)


or_est <- exp(coef(lizardsML))
or_ses <- sqrt(diag(diag(or_est) %*% vcov(lizardsML) %*% diag(or_est)))


ml_from_ml <- expo(lizardsML, type = "ML")
ml_from_br <- expo(lizardsBR, type = "ML")

med_from_br <- expo(lizardsBR, type = "AS_median")



library(parallel)
nsimu <- 10000
set.seed(123)
Y <- simulate(lizardsML, nsimu)
true_psi <- exp(coef(lizardsML))
X <- model.matrix(lizardsML)

methods <- c("ML", "correction+", "correction*", "Lylesetal2012", "AS_median")
results <- as.list(numeric(length(methods)))
names(results) <- methods
for (met in methods) {
    results[[met]] <- lapply(1:nsimu, function(k) {
        temp_data <- lizards
        temp_data[c("grahami", "opalinus")] <- Y[[k]]
        mod <- glm(cbind(grahami, opalinus) ~ height + diameter +
                       light + time, family = binomial(), data = temp_data,
                   method = "brglmFit", type = "ML", epsilon = 1e-10, maxit = 1000)
        expo(mod, met)
    }, mc.cores = 10)
}

get_bias <- function(res, truth) rowMeans(sapply(res, function(x) x$coef)) - true_psi
get_pu <- function(res, truth) rowMeans(sapply(res, function(x) x$coef < true_psi))

sapply(results, get_bias, truth = true_psi)
sapply(results, get_pu, truth = true_psi)

