context("tests for brmultinom and its agreement with other methods")

library("pmlr")
library("nnet")

#####################################################################
## Analysis of the enzymes data set in ?pmlr
#####################################################################
data("hepatitis", package = "pmlr")
## Construct a variable with the multinomial categories according to
## the HCV and nonABC columns
hepat <- hepatitis
hepat$type <- with(hepat, factor(1 - HCV * nonABC + HCV + 2 * nonABC))
hepat$type <- factor(hepat$type, labels = c("noDisease", "C", "nonABC"))
contrasts(hepat$type) <- contr.treatment(3, base = 1)


##############
heppmlr <- pmlr(type ~ group * time,
                data = hepat, weights = counts, method = "wald",
                penalized = TRUE)
hepbr <- brmultinom(type ~ group * time,
                    data = hepat, weights = counts)

tol <- 1e-05
test_that("brmultinom returns the same estimates as pmlr", {
    expect_equal(coef(hepbr), t(drop(coef(heppmlr))), tolerance = tol)
})

##############
hepbr_mat <- brmultinom(counts*nnet::class.ind(type) ~ group * time,
                        data = hepat)

test_that("brmultinom returns the same estimates if counts are supplied as a matrix", {
    expect_equal(coef(hepbr), coef(hepbr_mat), tolerance = tol)
})

#####################################################################
## Analysis of the enzymes data set in ?pmlr
#####################################################################
data(enzymes)
## Exclude patients in Group 4 (post-necrotic cirrhosis)
enzymes <- subset(enzymes, Group != 4)
## Center and scale covariates
AST <- scale(log(enzymes$AST))
ALT <- scale(log(enzymes$ALT))
GLDH <- scale(log(enzymes$GLDH))
OCT <- scale(log(enzymes$OCT))
enzymes <- data.frame(Patient = enzymes$Patient,
                      Group = enzymes$Group, AST, ALT, GLDH, OCT)
## Multinomial: acute viral hepatitis and aggressive chronic hepatits
## vs. persistent chronic hepatitis
## Assign Group 2 (persistent chronic hepatitis) as baseline category
enzymes$Group <- factor(enzymes$Group, levels=c("2","1","3"))
enzymes$counts <- rep(1, nrow(enzymes))
## enzpmlr <- pmlr(Group ~ AST + GLDH, weights = counts, data = enzymes, method = "wald")
## enzbrmultinom <- brmultinom(Group ~ AST + GLDH, weights = counts, data = enzymes)
enzbrmultinom_ml <- brmultinom(Group ~ AST + GLDH, weights = counts, data = enzymes, type = "maximum_likelihood")
enzmultinom <- nnet::multinom(Group ~ AST + GLDH, weights = counts, data = enzymes)

###############
test_that("brmultinom returns the same estimates as nnet::multinom if type = 'maximum_likelihood'", {
    expect_equal(coef(enzbrmultinom_ml), coef(enzmultinom), tolerance = 1e-04)
})

test_that("brmultinom returns the same fitted values as nnet::multinom if type = 'maximum_likelihood'", {
    expect_equal(fitted(enzbrmultinom_ml), fitted(enzmultinom), tolerance = 1e-04)
})

test_that("brmultinom returns the same maximized loglikelihood as nnet::multinom if type = 'maximum_likelihood'", {
    expect_equal(logLik(enzbrmultinom_ml), logLik(enzmultinom), tolerance = tol)
})

test_that("brmultinom returns the same standard errors as nnet::multinom if type = 'maximum_likelihood'", {
    expect_equal(summary(enzbrmultinom_ml)$standard.errors, summary(enzmultinom)$standard.errors, tolerance = 1e-04)
})


## Artificial example to check the effect of ignoring to explicitly
## pass the covariate classes with zero counts, when a category is not
## observed.

## library("nnet")

## ## Coefficients form glm
## report <- function(fit1, fit2, tol = 1e-03) {
##     if (inherits(fit1, "glm") & inherits(fit2, "glm")) {
##         coefs1 <- coef(fit1)
##         coefs2 <- coef(fit2)
##     }
##     if (inherits(fit1, "glm") & inherits(fit2, "multinom")) {
##         coefs2 <- as.matrix(coef(fit2))
##         q <- nrow(coefs2)
##         parNames <- c(paste0("resp", 2:(q+1)), paste0("resp", 2:(q+1), ":x1"), paste0("resp", 2:(q+1), ":x2"))
##         coefs1 <- matrix(coef(fit1)[parNames], nrow = nrow(coefs2))
##         dimnames(coefs1) <- dimnames(coefs2)
##     }
##     if (inherits(fit1, "multinom") & inherits(fit2, "glm")) {
##         coefs1 <- coef(fit1)
##         q <- nrow(coefs1)
##         parNames <- c(paste0("resp", 2:(q+1)), paste0("resp", 2:(q+1), ":x1"), paste0("resp", 2:(q+1), ":x2"))
##         coefs2 <- matrix(coef(fit2)[parNames], nrow = nrow(coefs1))
##         dimnames(coefs2) <- dimnames(coefs1)
##         out <- list(multinom1 = coefs1, multinom2 = coefs2, comparison = all.equal(coefs1, coefs2, tol = tol))
##     }
##     if (inherits(fit1, "multinom") & inherits(fit2, "multinom")) {
##         coefs1 <- coef(fit1)
##         coefs2 <- coef(fit2)
##         out <- list(multinom1 = coefs1, multinom2 = coefs2, comparison = all.equal(coefs1, coefs2, tol = tol))
##     }
##     ## if (inherits(fit1, "glm")) {
##     ##     parNames <- c(paste0("resp", 2:(q+1)), paste0("resp", 2:(q+1), ":x1"), paste0("resp", 2:(q+1), ":x2"))
##     ##     coefs1 <- matrix(coef(fit1)[parNames], nrow = nrow(coefs2))
##     ##     dimnames(coefs1) <- dimnames(coefs2)
##     ##     out <- list(glm = coefs1, multinom = coefs2, comparison = all.equal(coefs1, coefs2, tol = tol))
##     ## }
##     ## if (inherits(fit1, "multinom")) {
##     ##     coefs1 <- coef(fit1)
##     ##     out <- list(multinom1 = coefs1, multinom2 = coefs2, comparison = all.equal(coefs1, coefs2, tol = tol))
##     ## }
##     list(fit1 = coefs1,
##          fit2 = coefs2,
##          comparison = all.equal(coefs1, coefs2, tol = tol))
## }

## ## All covariate class - covariate combinations present
## set.seed(123)
## ncat <- 3
## N <- ncat*5
## eps <- 0.01

## set.seed(123)
## probs <- runif(ncat, eps, ncat - eps)
## probs <- probs/sum(probs)
## totals <- 10
## counts <- c(rmultinom(N, totals, probs))
## resp <- factor(rep(seq.int(ncat), N))
## ## Simulate covariate values and repeat each value ncat times
## x1 <- rep(rnorm(N/ncat, 2, 2), each = ncat)
## x2 <- rep(rexp(N/ncat, 10), each = ncat)
## ## Indicator of covariate classes
## inds <- factor(rep(seq.int(N), each = ncat))
## ## Collect everything in a data frame
## artificialData0 <- data.frame(resp = resp, x1 = x1, x2 = x2, counts = counts, inds = inds)
## ## Clean up
## rm(list = c("counts", "resp", "x1", "x2", "inds"))

## ## Model fitting
## fit_multinom0 <- multinom(resp ~ x1 + x2, weights = counts, data = artificialData0)
## fit_glm0 <- glm(counts ~ -1 + inds + resp * (x1 + x2), data = artificialData0, family = poisson())

## ## Effect of removing a few observations
## set.seed(222)
## del <- sample(seq.int(nrow(artificialData0)), 5)
## artificialData1 <- artificialData0[-del, ]
## ## Model fitting
## fit_multinom1 <- multinom(resp ~ x1 + x2, weights = counts, data = artificialData1)
## fit_glm1 <- glm(counts ~ -1 + inds + resp * (x1 + x2), data = artificialData1, family = poisson())

## ## Now set zero counts for the missing category - covariate class combinations and refit
## artificialData1a <- within(artificialData0, {
##     counts[del] <- 0
## })
## fit_multinom1a <- multinom(resp ~ x1 + x2, weights = counts, data = artificialData1a)
## fit_glm1a <- glm(counts ~ -1 + inds + resp * (x1 + x2), data = artificialData1a, family = poisson())




## ## An alternative that does not require a specific structure from the data!
## cdat <- artificialData1
## ## Construct model matrix with the covariates
## Xinterest <- model.matrix(~ x1 + x2, data = cdat)
## ## Construct the matrix with the nuisances
## Xnuisance <- Diagonal(nrow(Xinterest))
## ## Construct the
## Xextended <- cbind(kronecker(rep(1, nlevels(cdat$resp)), Xnuisance),
##                    kronecker(Diagonal(nlevels(cdat$resp))[, -1], Xinterest))

## int <- seq.int(nrow(Xinterest))
## nd <- paste0("%0", nchar(max(int)), "d")
## colnames(Xextended) <- c(paste0(".nuisance", sprintf(nd, int)),
##                          paste(rep(levels(cdat$resp)[-1], each = ncol(Xinterest)),
##                                rep(colnames(Xinterest), nlevels(cdat$resp) - 1), sep = ":"))
## countsExtended <- c(nnet::class.ind(artificialData1$resp) * artificialData1$counts)
## fit_glm2a <- glm.fit(Xextended, countsExtended, family = poisson())


## ## fit_glm0 and fit_multinom0 are fits based on artificialData0 which
## ## has all covariate class - response combinations
## report(fit_glm0, fit_multinom0)

## ## fit_glm1 and fit_multinom1 are fits after removing the observations
## ## in del from artificialData0 (FALSE)
## report(fit_glm1, fit_multinom1)

## ## fit_glm1a and fit_multinom1a are fits after adding the covariate
## ## class - response combinations with zero counts (TRUE)
## report(fit_glm1a, fit_multinom1a)

## ## fit_multinom1a and fit_multinom1 (should be TRUE)
## report(fit_multinom1, fit_multinom1a)


## cglm also breaks because of the zero count
## source("~/Repositories/cglm/cglm5.R")
## source("~/Repositories/cglm/multinom.fit.R")
## source("~/Repositories/cglm/cglmContrastMatrix.R")
## source("~/Repositories/cglm/deviance.R")
## cglm(counts*nnet::class.ind(resp) ~ x1 + x2, data = artificialData1, method = "multinom.fit")
