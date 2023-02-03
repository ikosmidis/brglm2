## source(system.file("brpr", "brpr.R", package = "brglm2"))
source("brpr.R")

## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
dobson <- data.frame(counts, outcome, treatment)

tol <- 1e-06
## BR estimates and std. errors from brglmFit and from brpr are the same for poisson
expect_warning({
    br1 <- summary(glm(counts ~ outcome + treatment, family = poisson(), method = "brglmFit"))
    br2 <- summary(brpr(counts ~ outcome + treatment, data = dobson))
})
c1 <- coef(br1)
c2 <- coef(br2)
c2 <- c2[rownames(c1), ]
expect_equal(c1,c2, tolerance = tol)


## ## mbest::firthglm.fit crashes for this example
## ## A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
## clotting <- data.frame(
##     u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
##     conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
##     lot = factor(c(rep(1, 9), rep(2, 9))))
## mod <- glm(conc ~ lot*log(u), data = clotting, family = Gamma())

## tol <- sqrt(.Machine$double.eps)
## test_that("BR estimates from brglmFit and from mbest::firthglm.fit are the same for gamma", {
##     br1 <- update(mod, method = "firthglm.fit")
##     br2 <- update(mod, method = "brglmFit")
##     c1 <- coef(br1)
##     c2 <- coef(br2)
##     c2 <- c2[names(c1)]
##     expect_equal(c1,c2, tolerance = tol)
## })


