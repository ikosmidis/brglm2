context("comparison for Poisson with mbrpr from the unpublished undergraduate 
         thesis of M. Zambelli (2016) of  the University of Padova ")

source(system.file("mbrpr", "mbrpr.R", package = "brglm2"))

## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
dobson <- data.frame(counts, outcome, treatment)

tol <- sqrt(.Machine$double.eps)

test_that("MBR estimates and std. errors from brglmFit and from mbrpr are the same for poisson", {
  br1 <- summary(fit<-glm(counts ~ outcome + treatment, family = poisson(), method = "brglmFit",
                     type="AS_median"))
  X <- model.matrix(fit)
  y <- fit$y
  br2 <- mbrpr(par=rep(0,ncol(X)),y=y,X=X,eps = 1e-10,maxit=500)
  c1 <- coef(br1)[,c(1,2)]
  dimnames(c1)=NULL
  c2 <- cbind(br2$coefficients,sqrt(diag(br2$InfoInv)))
  expect_equal(c1,c2, tolerance = tol)
})

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


