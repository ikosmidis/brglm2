## source(system.file("inst", "brglm0/brglm0.R", package = "brglm2"))
data("lizards", package = "brglm2")


types <- list("ML" = "(maximum likelihood)",
               "correction" = "(bias correction)",
               "AS_mean" = "(mean bias-reducing adjusted score equations)",
               "AS_median" = "(median bias-reducing adjusted score equations)",
               "AS_mixed" = "(mixed bias-reducing adjusted score equations)",
               "MPL_Jeffreys" = "(maximum Jeffreys' prior penalized likelihood)")


for (type in names(types)) {
    liz <- glm(cbind(grahami, opalinus) ~ height + diameter + light + time, family = binomial(),
               data = lizards,
               method = "brglmFit", type = type)
    expect_stdout(print(liz), "Coefficients")
    expect_stdout(print(liz), "Degrees of Freedom")
    expect_stdout(print(liz), "Null Deviance")
    expect_equal(liz$type, type)
    summ <- summary(liz)
    expect_true(all(class(summ) %in% c("summary.brglmFit", "summary.glm")))
    expect_stdout(print(summ), "Type of estimator:")
    expect_stdout(print(summ), type)
    expect_stdout(print(summ), types[[type]])
}


## brnb
salmonella <- data.frame(freq = c(15, 16, 16, 27, 33, 20,
                                   21, 18, 26, 41, 38, 27,
                                   29, 21, 33, 60, 41, 42),
                         dose = rep(c(0, 10, 33, 100, 333, 1000), 3),
                         observation = rep(1:3, each = 6))
salmonella_fm <- freq ~ dose + log(dose + 10)
fit_brnb <- brnb(salmonella_fm, data = salmonella,
                 link = "log", transformation = "inverse", type = "ML")
summ <- summary(fit_brnb)
expect_stdout(print(summ), "Type of estimator:")
expect_stdout(print(summ), "ML")
expect_stdout(print(summ), "(maximum likelihood)")

## brmultinom
data("housing", package = "MASS")
fit_brmultinom <- brmultinom(Sat ~ Infl + Type + Cont, weights = Freq,
                             data = housing, type = "ML", ref = 1)
summ <- summary(fit_brmultinom)
expect_stdout(print(summ), "Type of estimator:")
expect_stdout(print(summ), "ML")
expect_stdout(print(summ), "(maximum likelihood)")


## bracl
data("stemcell", package = "brglm2")
fit_bracl <- bracl(research ~ as.numeric(religion) + gender, weights = frequency,
                    data = stemcell, type = "ML")
summ <- summary(fit_bracl)
expect_stdout(print(summ), "Type of estimator:")
expect_stdout(print(summ), "ML")
expect_stdout(print(summ), "(maximum likelihood)")
