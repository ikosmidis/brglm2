library("numDeriv")
library("brglm")

## A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
    lot = factor(c(rep(1, 9), rep(2, 9))))
mod <- enrich(glm(conc ~ lot*log(u), data = clotting, family = Gamma()))

ifun <- get_information_function(mod)

penloglik <- function(pars, X, y) {
    p <- ncol(X)
    beta <- pars[1:p]
    dispersion <- pars[p + 1]
    eta <- X %*% beta
    mu <- 1/eta
    sum(dgamma(y, shape = 1/dispersion, scale = mu * dispersion, log = TRUE)) + 0.5 * log(det(ifun(beta, dispersion)))
}

modjef <- update(mod, method = "brglmFit", type = "MPL_Jeffreys", epsilon = 1e-15)

## the numerical gradient of the penalized log-likelihood is almost zero when evaluated at the estimates from type = 'MPL_Jeffreys
expect_equal(grad(penloglik,
                  x = coef(modjef, model = "full"),
                  X = model.matrix(mod),
                  y = model.response(mod$model)), rep(0, 5), tolerance = 1.5 * 1e-05)

## the numerical gradient of the penalized log-likelihood matches that from type = 'MPL_Jeffreys'
expect_warning(g1 <- update(mod, method = "brglmFit", type = "MPL_Jeffreys", maxit = 0, start = coef(mod, model = "full"))$grad)
g2 <- grad(penloglik,
           x = coef(mod, model = "full"),
           X = model.matrix(mod),
           y = model.response(mod$model))
expect_equal(unname(g1), g2, tolerance = 1e-05)

## source(system.file("inst", "brglm0/brglm0.R", package = "brglm2"))
data("lizards", package = "brglm2")

links <- lapply(c("logit", "probit", "cloglog", "cauchit"), make.link)


tol <- 1e-08
for (l in seq_along(links)) {
    expect_warning(
        lizardsBRlegacy <- brglm(cbind(grahami, opalinus) ~ height + diameter +
                                     light + time, family = binomial(links[[l]]), data=lizards,
                                 method = "brglm.fit", br.epsilon = 1e-10, br.maxit = 1000, pl = TRUE)
    )

    expect_warning(
        lizardsBR <- glm(cbind(grahami, opalinus) ~ height + diameter +
                             light + time, family = binomial(links[[l]]), data=lizards,
                         method = "brglmFit", epsilon = 1e-10, maxit = 1000, type = "MPL_Jeffreys")
    )
    ## glm with brglm.fit method and brglm_0 return the same coefficients for the lizards when link is links[[l]]$name)
    expect_equal(coef(lizardsBR), coef(lizardsBRlegacy), tolerance = tol)
}
