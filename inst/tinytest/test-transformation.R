## A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
    lot = factor(c(rep(1, 9), rep(2, 9))))
mod_identity <- glm(conc ~ lot*log(u), data = clotting, family = Gamma,
                    method = "brglm_fit", type = "ML",
                    transformation = "identity")
mod_log <- glm(conc ~ lot*log(u), data = clotting, family = Gamma,
                    method = "brglm_fit", type = "ML",
                    transformation = "log")
mod_sqrt <- glm(conc ~ lot*log(u), data = clotting, family = Gamma,
                    method = "brglm_fit", type = "ML",
                    transformation = "sqrt")
mod_inverse <- glm(conc ~ lot*log(u), data = clotting, family = Gamma,
                    method = "brglm_fit", type = "ML",
                    transformation = "inverse")

c_identity <- coef(mod_identity, model = "full")
c_log <- coef(mod_log, model = "full")
c_log[5] <- exp(c_log[5])
c_sqrt <- coef(mod_sqrt, model = "full")
c_sqrt[5] <- c_sqrt[5]^2
c_inverse <- coef(mod_inverse, model = "full")
c_inverse[5] <- 1/c_inverse[5]

tol <- sqrt(.Machine$double.eps)
## ML estimate of gamma dispersion from brglmFit is invariant to trasnformation
expect_equal(c_identity, c_log, tolerance = tol, check.attributes = FALSE)
expect_equal(c_identity, c_sqrt, tolerance = tol, check.attributes = FALSE)
expect_equal(c_identity, c_inverse, tolerance = tol, check.attributes = FALSE)
