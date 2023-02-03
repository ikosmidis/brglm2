data("anorexia", package = "MASS")

anorexML <- glm(Postwt ~ Prewt + Treat+ offset(Prewt),
                family = gaussian, data = anorexia)
anorexBR <- update(anorexML, method = "brglmFit")
anorexMBR <- update(anorexML, method = "brglmFit", control = list(type="AS_median"))

tol <- sqrt(.Machine$double.eps)

## dispersion is RSS over residual degrees of freedom minus 2/3"
expect_equal(anorexMBR$dispersion, sum((anorexia$Postwt - fitted(anorexML))^2)/(nrow(anorexia) - length(coef(anorexML)) - 2/3))

