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
data("enzymes", package = "pmlr")
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
enzmultinom <- nnet::multinom(Group ~ AST + GLDH, weights = counts, data = enzymes, trace = FALSE)

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


hepbr1 <- brmultinom(type ~ group * time,
                    data = hepat, weights = counts, ref = 1)
hepbr2 <- brmultinom(type ~ group * time,
                    data = hepat, weights = counts, ref = 2)
hepbr3 <- brmultinom(type ~ group * time,
                     data = hepat, weights = counts, ref = 3)
test_that("brmultinom fits are invariant to the value of ref (ref1 vs ref2)'", {
    expect_equal(hepbr1$fitted.values, hepbr2$fitted.values, tolerance = 1e-04)
})

test_that("brmultinom fits are invariant to the value of ref  (ref1 vs ref3)'", {
    expect_equal(hepbr1$fitted.values, hepbr3$fitted.values, tolerance = 1e-04)
})

