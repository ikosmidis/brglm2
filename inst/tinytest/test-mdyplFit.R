library("brglm")
data("lizards", package = "brglm2")

a <- 0.3

liz_fm <- cbind(grahami, opalinus) ~ height + diameter + light + time
expect_error(
    glm(liz_fm, family = binomial("probit"), data = lizards,
        method = "mdyplFit", alpha = a)
)

expect_error(
    glm(liz_fm, family = gaussian(), data = lizards,
        method = "mdyplFit", alpha = a)
)

liz_fit_DY <- glm(liz_fm, family = binomial(), data = lizards,
                  method = "mdyplFit", alpha = a)

adj_lizards <- lizards |>
    transform(totals = grahami + opalinus) |>
    transform(g_y = grahami / totals) |>
    transform(grahami = alpha * grahami + totals * (1 - alpha) / 2) |>
    transform(opalinus = totals - grahami)
temp_fit <- glm(cbind(grahami, opalinus) ~ height + diameter +
                        light + time, family = binomial(),
                data = adj_lizards)

## Correct alpha
expect_equal(a, liz_fit_DY$alpha)

## Correct coefficients
expect_equal(coef(liz_fit_DY), coef(temp_fit))

## Correct deviance
yy <- adj_lizards$g_y
mm <- fitted(temp_fit)
tt <- adj_lizards$totals
expect_equal(deviance(liz_fit_DY),
             sum(binomial()$dev.resids(yy, mm, tt)))

## Correct AIC
expect_equal(AIC(liz_fit_DY),
             binomial()$aic(yy, tt, mm, tt) + 2 * length(coef(temp_fit)))

## Correct log-likelihood
expect_equal(as.vector(logLik(liz_fit_DY)),
             - binomial()$aic(yy, tt, mm, tt) / 2)

## Correct estimated standard errors
v <- tt * mm * (1 - mm)
xx <- model.matrix(liz_fit_DY)
vv <- solve(crossprod(xx * sqrt(v)))
expect_equal(sqrt(diag(vv)), coef(summary(liz_fit_DY))[, "Std. Error"])

##
