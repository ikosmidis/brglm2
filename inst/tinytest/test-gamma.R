## A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
    lot = factor(c(rep(1, 9), rep(2, 9))))
mod <- glm(conc ~ lot*log(u), data = clotting, family = Gamma)

tol <- sqrt(.Machine$double.eps)
## ML estimate of gamma dispersion from brglmFit and from MASS::gamma.dispersion are the same
expect_equal(update(mod, method = "brglmFit", epsilon = 1e-10)$dispersion_ML, MASS::gamma.dispersion(mod), tolerance = tol)

## ML estimate of gamma shape from brglmFit and from MASS::gamma.dispersion are the same
expect_equal(1/update(mod, method = "brglmFit", epsilon = 1e-10, transformation = "inverse")$dispersion_ML, MASS::gamma.shape(mod)$alpha, tolerance = tol)
