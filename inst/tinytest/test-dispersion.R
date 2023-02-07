data("anorexia", package = "MASS")

anorexML <- glm(Postwt ~ Prewt + Treat+ offset(Prewt),
                family = gaussian, data = anorexia)
anorexBR <- update(anorexML, method = "brglmFit", type = "AS_mean")

tol <- sqrt(.Machine$double.eps)
## dispersion_ML is the usual biased estimate for the residual variance
expect_equal(anorexBR$dispersion_ML, sum((anorexia$Postwt - fitted(anorexML))^2)/nrow(anorexia),
             tolerance = 1e-06)

## dispersion is the usual bias-corrected estimate for the residual variance
expect_equal(anorexBR$dispersion, sum((anorexia$Postwt - fitted(anorexML))^2)/(nrow(anorexia) - length(coef(anorexML))), tolerance = 1e-06)


## context("dispersion parameter estimation")

## set.seed(123)
## N <- 20
## x <- matrix(rnorm(2*N), N, 2)
## y <- 0.3 + drop(x %*% c(-3, 1)) + rnorm(N, 0.4)

## fitML <- glm(y ~ x, family = gaussian())
## fitBR <- update(fitML, method = "brglmFit")

## tol <- sqrt(.Machine$double.eps)
## test_that("dispersionML is the usual biased estimate for the residual variance",
##           expect_equal(fitBR$dispersionML, sum((fitML$y - fitted(fitML))^2)/N))

## test_that("dispersion is the usual bias-corrected estimate for the residual variance",
##           expect_equal(fitBR$dispersion, sum((fitML$y - fitted(fitML))^2)/(N - 3)))






