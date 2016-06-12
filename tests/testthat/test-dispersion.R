context("dispersion parameter estimation")

utils::data(anorexia, package = "MASS")
anorexML <- glm(Postwt ~ Prewt + Treat + offset(Prewt),
                family = gaussian, data = anorexia)
anorexBR <- update(anorexML, method = "brglmFit")

tol <- sqrt(.Machine$double.eps)
test_that("dispersionML is the usual biased estimate for the residual variance",
          expect_equal(anorexBR$dispersionML, sum((anorexia$Postwt - fitted(anorexML))^2)/nrow(anorexia)))

test_that("dispersion is the usual bias-corrected estimate for the residual variance",
          expect_equal(anorexBR$dispersion, sum((anorexia$Postwt - fitted(anorexML))^2)/(nrow(anorexia) - length(coef(anorexML)))))



