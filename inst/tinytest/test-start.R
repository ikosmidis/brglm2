## A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
    lot = factor(c(rep(1, 9), rep(2, 9))))
mod <- glm(conc ~ lot*log(u), data = clotting, family = Gamma, epsilon = 1e-10, maxit = 1000)

mod1 <- update(mod, method = "brglmFit", start = c(coef(mod)*0.9, 1), epsilon = 1e-10, maxit = 1000)
mod2 <- update(mod, method = "brglmFit", start = c(coef(mod)*0.9), epsilon = 1e-10, maxit = 1000)

tol <- 1e-03
## start argument is passed correctly in brglmFit"
expect_equal(coef(mod), coef(mod1), tolerance = tol)
expect_equal(coef(mod), coef(mod2), tolerance = tol)
