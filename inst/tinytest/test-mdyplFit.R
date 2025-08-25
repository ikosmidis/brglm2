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
    transform(grahami = a * grahami + totals * (1 - a) / 2) |>
    transform(opalinus = totals - grahami)
expect_warning(
    temp_fit <- glm(cbind(grahami, opalinus) ~ height + diameter +
                        light + time, family = binomial(),
                    data = adj_lizards)
)

## Correct alpha
expect_equal(a, liz_fit_DY$alpha)

## Correct coefficients
expect_equal(coef(liz_fit_DY), coef(temp_fit))

## Correct deviance
yy <- adj_lizards$g_y
mm <- fitted(temp_fit)
tt <- adj_lizards$totals
expect_equal(deviance(liz_fit_DY),
             sum(binomial()$dev.resids(a * yy + (1- a)/2, mm, tt)))

## Correct estimated standard errors
v <- tt * mm * (1 - mm)
xx <- model.matrix(liz_fit_DY)
vv <- solve(crossprod(xx * sqrt(v)))
expect_equal(sqrt(diag(vv)), coef(summary(liz_fit_DY))[, "Std. Error"])

## Naming conventions
liz_fit_DY1 <- update(liz_fit_DY, method = "mdypl_fit")
objs <- names(liz_fit_DY)
objs <- objs[!(objs %in% c("call", "method"))]
for (what in objs) {
    expect_equal(liz_fit_DY1[[what]],
                 liz_fit_DY[[what]])
}


## taus
xx <- model.matrix(liz_fit_DY)
expect_equal(
    sapply(2:ncol(xx), function(j) sum( lm.fit(xx[, -c(1, j)], xx[, j])$residuals^2 ) / (nrow(xx) - ncol(xx) + 2)) |> sqrt(),
    brglm2:::taus(liz_fit_DY))


## logist_aic
set.seed(111)
tots <- c(11, 15, 5, 5, 15, 111, 8)
succ <- c(0, 11, 3, 5, 2, 5, 7)
fail <- tots - succ
probs <- c(0, runif(length(fail) - 2), 1)
expect_equal(brglm2:::logist_aic(succ / tots, tots, probs, tots),
             binomial()$aic(succ / tots, tots, probs, tots))
expect_equal(-brglm2:::logist_aic(succ / tots, tots, probs, tots) / 2,
             sum(dbinom(succ, tots, probs, log = TRUE)))


expect_equal(brglm2:::dbinom2(succ, tots, probs, log = FALSE),
             dbinom(succ, tots, probs, log = FALSE))

expect_equal(brglm2:::dbinom2(succ, tots, probs, log = TRUE),
             dbinom(succ, tots, probs, log = TRUE))

expect_equal(brglm2:::dbinom2(1, 1, 0), dbinom(1, 1, 0))
expect_equal(brglm2:::dbinom2(1, 1, 1), dbinom(1, 1, 1))
expect_equal(brglm2:::dbinom2(0, 1, 0), dbinom(0, 1, 0))
expect_equal(brglm2:::dbinom2(0, 1, 1), dbinom(0, 1, 1))
expect_equal(brglm2:::dbinom2(1, 0, 0), dbinom(1, 0, 0))
expect_equal(brglm2:::dbinom2(1, 0, 1), dbinom(1, 0, 1))
expect_equal(brglm2:::dbinom2(0, 0, 1), dbinom(0, 0, 1))
expect_equal(brglm2:::dbinom2(0, 0, 0), dbinom(0, 0, 0))
expect_equal(brglm2:::dbinom2(1, 1, 0), dbinom(1, 1, 0), log = TRUE)
expect_equal(brglm2:::dbinom2(1, 1, 1), dbinom(1, 1, 1), log = TRUE)
expect_equal(brglm2:::dbinom2(0, 1, 0), dbinom(0, 1, 0), log = TRUE)
expect_equal(brglm2:::dbinom2(0, 1, 1), dbinom(0, 1, 1), log = TRUE)
expect_equal(brglm2:::dbinom2(1, 0, 0), dbinom(1, 0, 0), log = TRUE)
expect_equal(brglm2:::dbinom2(1, 0, 1), dbinom(1, 0, 1), log = TRUE)
expect_equal(brglm2:::dbinom2(0, 0, 1), dbinom(0, 0, 1), log = TRUE)
expect_equal(brglm2:::dbinom2(0, 0, 0), dbinom(0, 0, 0), log = TRUE)

