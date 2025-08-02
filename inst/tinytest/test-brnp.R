library("MASS")

salmonella <- data.frame(freq = c(15, 16, 16, 27, 33, 20,
                                  21, 18, 26, 41, 38, 27,
                                  29, 21, 33, 60, 41, 42),
                         dose = rep(c(0, 10, 33, 100, 333, 1000), 3),
                         observation = rep(1:3, each = 6))


salmonella_fm <- freq ~ dose + log(dose+10)

fitML_glmnb <- glm.nb(salmonella_fm, data = salmonella)
fitML <- brnb(salmonella_fm, data = salmonella, link = "log", transformation = "identity", type = "ML")
fitBR_mean <- update(fitML, type = "AS_mean")
fitBR_median <- update(fitML, type = "AS_median")
fitBR_mixed <- update(fitML, type = "AS_mixed")
fitBC_mean <- update(fitML, type = "correction")

## brnp returns the same estimates as MASS::glm.nb
expect_equal(coef(fitML), coef(fitML_glmnb), tol = 1e-06)

## brnp returns mean BR and median BR estimates as in Kenne Pagui et al (2021+)
expect_equal(coef(fitBR_mean, model = "full"),
             c(2.216, -0.001, 0.309, 0.065),
             check.attributes = FALSE, tol = 1e-03)
expect_equal(coef(fitBR_median, model = "full"),
             c(2.211, -0.001, 0.309, 0.069),
             check.attributes = FALSE, tol = 1e-03)
expect_equal(coef(fitBC_mean, model = "full"),
             c(2.210, -0.001, 0.311, 0.063),
             check.attributes = FALSE, tol = 1e-03)

## fit using weights

duptimes <- c(2,1,3,5,rep(1,14))
idx <- rep(1:nrow(salmonella), duptimes)
dupsalmonella <- salmonella[idx,]

fitMLe <- update(fitML, data = dupsalmonella)
fitMLw <- update(fitML, weights = duptimes, data = salmonella)

fitBR_mediane <- update(fitML, data = dupsalmonella, type = "AS_median")
fitBR_medianw <- update(fitML, weights = duptimes, data = salmonella, type = "AS_median")

fitBR_meane <- update(fitML, data = dupsalmonella, type = "AS_mean")
fitBR_meanw <- update(fitML, weights = duptimes, data = salmonella, type = "AS_mean")

fitBR_mixede <- update(fitML, data = dupsalmonella, type = "AS_mixed")
fitBR_mixedw <- update(fitML, weights = duptimes, data = salmonella, type = "AS_mixed")

fitBCe <- update(fitML, data = dupsalmonella, type = "AS_mixed")
fitBCw <- update(fitML, weights = duptimes, data = salmonella, type = "AS_mixed")

fitJe <- update(fitML, data = dupsalmonella, type = "MPL_Jeffreys",
                transformation = "inverse")
fitJw <- update(fitML, weights = duptimes, data = salmonella, type = "MPL_Jeffreys",
                transformation = "inverse")

## prior weights work as expected
## all numerical results are the same
expect_equal(coef(fitMLe, "full"), coef(fitMLw, "full"), tolerance = 1e-10)
expect_equal(coef(fitBR_mediane, "full"), coef(fitBR_medianw, "full"), tolerance = 1e-10)
expect_equal(coef(fitBR_mixede, "full"), coef(fitBR_mixedw, "full"), tolerance = 1e-10)
expect_equal(coef(fitBR_meane, "full"), coef(fitBR_meanw, "full"), tolerance = 1e-10)
expect_equal(coef(fitJe, "full"), coef(fitJw, "full"), tolerance = 1e-10)

## Dispersion transformations
## dispersion transformations work as expected for ML/mixed BR/median BR
for (f0 in list(fitML, fitBR_median, fitBR_mixed)) {
    fsqrt <- update(f0, transformation = "sqrt")
    flog <- update(f0, transformation = "log")
    expect_equal(coef(fsqrt, model = "dispersion"),
                 sqrt(coef(f0, model = "dispersion")),
                 tol = 1e-05, check.attributes = FALSE)
    expect_equal(coef(flog, model = "dispersion"),
                 log(coef(f0, model = "dispersion")),
                 tol = 1e-05, check.attributes = FALSE)
}


## dispersion transformations work as expected for mean BR
f0 <- fitBR_mean
fsqrt <- update(f0, transformation = "sqrt")
flog <- update(f0, transformation = "log")
expect_false(isTRUE(all.equal(coef(fsqrt, model = "dispersion"),
                              sqrt(1 / coef(f0, model = "dispersion")),
                              tol = 1e-04, check.attributes = FALSE)))
expect_false(isTRUE(all.equal(coef(flog, model = "dispersion"),
                              -log(coef(f0, model = "dispersion")),
                              tol = 1e-04, check.attributes = FALSE)))

## error is produced for not implemented transformations
expect_error(update(fitBR_mean, transformation = "asd"))



## From https://github.com/ikosmidis/brglm2/issues/31
## Infinite estimates
## df1 <- data.frame(production=c(15, 12, 10, 9, 6, 8, 9, 5, 3, 3, 2, 1, 0, 0, 0, 0),
##                   Treatment_Num=c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4),
##                   Treatment=c("Fresh", "Fresh", "Fresh", "Fresh", "low", "low", "low", "low", "med", "med", "med", "med", "high", "high", "high", "high"),
##                   Genotype=c(1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2),
##                   Source_Salinity=c("Fresh", "Fresh", "Brackish", "Brackish", "Fresh", "Fresh", "Brackish", "Brackish", "Fresh", "Fresh", "Brackish", "Brackish", "Fresh", "Fresh", "Brackish", "Brackish"),
##                   Days_to_death <- c(500,  500,  500,  500,  400,  350,  300,  500,  200,  202,  260,  280,  150, 150, 160, 140),
##                   censored <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0))


## salmonella_fm <- production ~ Treatment
## mbr_inv <- brglm2::brnb(salmonella_fm, data = df1, link = "log", transformation = "inverse", type = "AS_mean")
## mbr_log <- brglm2::brnb(salmonella_fm, data = df1, link = "log", transformation = "log", type = "AS_mean")
## transformation = "identity" fails with default starting values
## mbr_sqrt <- brglm2::brnb(salmonella_fm, data = df1, link = "log", transformation = "identity", type = "AS_mean")

