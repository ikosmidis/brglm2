## Tolerance for comparisons
tolerance <- 1e-05

data("lizards")

lizards_grahami <- lizards[, c("grahami", "height", "diameter", "light", "time")]
lizards_grahami <- lizards_grahami[rep(seq.int(nrow(lizards_grahami)), lizards_grahami$grahami), ]
lizards_grahami$species <- "grahami"
lizards_grahami$grahami <- NULL
lizards_opalinus <- lizards[, c("opalinus", "height", "diameter", "light", "time")]
lizards_opalinus <- lizards_opalinus[rep(seq.int(nrow(lizards_opalinus)), lizards_opalinus$opalinus), ]
lizards_opalinus$species <- "opalinus"
lizards_opalinus$opalinus <- NULL
lizards1 <- rbind(lizards_grahami, lizards_opalinus)
lizards1$species <- factor(lizards1$species, levels = c("opalinus", "grahami"))

model1 <- glm(formula = species ~ height + diameter + light + time, family = binomial(logit),
              data = lizards1,
              method = "brglmFit")
model2 <- brmultinom(species ~ height + diameter + light + time, data = lizards1)

tolerance <- 1e-06
## coefficients from the binomial fit match those from the multinomial model with the poisson trick
expect_equal(coef(model1), coef(model2)[1,], tol = tolerance)

