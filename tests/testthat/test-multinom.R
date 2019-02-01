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
expect_warning(
    hepbr <- brmultinom(type ~ group * time,
                        data = hepat, weights = counts)
)

tol <- 1e-05
test_that("brmultinom returns the same estimates as pmlr", {
    expect_equal(coef(hepbr), t(drop(coef(heppmlr))), tolerance = tol)
})

##############
expect_warning(
    hepbr_mat <- brmultinom(counts*nnet::class.ind(type) ~ group * time,
                            data = hepat)
)

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
expect_warning(
    enzbrmultinom_ml <- brmultinom(Group ~ AST + GLDH, weights = counts, data = enzymes, type = "ML")
)
enzmultinom <- nnet::multinom(Group ~ AST + GLDH, weights = counts, data = enzymes, trace = FALSE)

###############
test_that("brmultinom returns the same estimates as nnet::multinom if type = 'ML'", {
    expect_equal(coef(enzbrmultinom_ml), coef(enzmultinom), tolerance = 1e-04)
})

test_that("brmultinom returns the same fitted values as nnet::multinom if type = 'ML'", {
    expect_equal(fitted(enzbrmultinom_ml), fitted(enzmultinom), tolerance = 1e-04)
})

test_that("brmultinom returns the same model matrix as nnet::multinom", {
    expect_equal(model.matrix(enzbrmultinom_ml) ,model.matrix(enzmultinom))
})

test_that("brmultinom returns the same maximized loglikelihood as nnet::multinom if type = 'ML'", {
    expect_equal(logLik(enzbrmultinom_ml), logLik(enzmultinom), tolerance = tol)
})

test_that("brmultinom returns the same standard errors as nnet::multinom if type = 'ML'", {
    expect_equal(summary(enzbrmultinom_ml)$standard.errors, summary(enzmultinom)$standard.errors, tolerance = 1e-04)
})


expect_warning({
    hepbr1 <- brmultinom(type ~ group * time,
                         data = hepat, weights = counts, ref = 1)
    hepbr2 <- brmultinom(type ~ group * time,
                         data = hepat, weights = counts, ref = 2)
    hepbr3 <- brmultinom(type ~ group * time,
                         data = hepat, weights = counts, ref = 3)
})
test_that("brmultinom fits are invariant to the value of ref (ref1 vs ref2)'", {
    expect_equal(hepbr1$fitted.values, hepbr2$fitted.values, tolerance = 1e-04)
})

test_that("brmultinom fits are invariant to the value of ref  (ref1 vs ref3)'", {
    expect_equal(hepbr1$fitted.values, hepbr3$fitted.values, tolerance = 1e-04)
})

newdata <- data.frame(group = c("no-withhold", "withhold", "no-withhold", "withhold"),
                      time = c("pre", "pre", "post", "post"))

hepnnet <- multinom(type ~ group * time, data = hepat, weights = counts, trace = FALSE)
expect_warning({
    hepml1 <- brmultinom(type ~ group * time,
                         data = hepat, weights = counts, ref = 1, type = "ML", maxit = 5)
    hepml2 <- brmultinom(type ~ group * time,
                         data = hepat, weights = counts, ref = 2, type = "ML", maxit = 5)
    hepml3 <- brmultinom(type ~ group * time,
                         data = hepat, weights = counts, ref = 3, type = "ML", maxit = 5)
})

test_that("predict.brmultinom returns the right result", {
    expect_equal(predict(hepnnet), predict(hepml1), tolerance = 1e-04)
    expect_equal(predict(hepnnet), predict(hepml2), tolerance = 1e-04)
    expect_equal(predict(hepnnet), predict(hepml3), tolerance = 1e-04)
    expect_equal(predict(hepnnet, newdata = newdata), predict(hepml1, newdata = newdata), tolerance = 1e-04)
    expect_equal(predict(hepml1, newdata = newdata), predict(hepml2, newdata = newdata), tolerance = 1e-04)
    expect_equal(predict(hepml1, newdata = newdata), predict(hepml3, newdata = newdata), tolerance = 1e-04)
})

## Aligator data
## data("alligators", package = "brglm2")
## k <- 3
## all_ml <- brmultinom(foodchoice ~ size + lake , weights = round(freq/k),
##                      data = alligators, type = "ML", ref = 1)

## all_mean <- brmultinom(foodchoice ~ size + lake , weights = round(freq/k),
##                        data = alligators, type = "AS_mean", ref = 1)

## all_median <- brmultinom(foodchoice ~ size + lake , weights = round(freq/k),
##                          data = alligators, type = "AS_median", ref = 1)

## library("ggplot2")

## ## Collect probabilities
## pml <- data.frame(stack(data.frame(fitted(all_ml))), method = "ML")
## pmean <- data.frame(stack(data.frame(fitted(all_mean))), method = "AS_mean")
## pmedian <- data.frame(stack(data.frame(fitted(all_median))), method = "AS_median")
## probabilities <- rbind(pml, pmean, pmedian)
## names(probabilities) <- c("probability", "category", "method")
## probabilities1 <- NULL
## for (j in levels(probabilities$category)) {
##     cdat <- subset(probabilities, category == j)
##     cdat$id <- rep(seq.int(nrow(pml)/5), 3)
##     cdat <- reshape(cdat, timevar = "method", v.names = "probability", direction = "wide")
##     probabilities1 <- rbind(probabilities1, cdat)
## }

## ## Shrinkage plots
## ggplot(probabilities1) +
##     geom_point(aes(x = probability.ML, y = probability.AS_mean)) +
##     geom_abline(aes(intercept = 0, slope = 1), alpha = 0.5) +
##     lims(x = c(0,1), y = c(0,1)) +
##     facet_grid(~ category) +
##     theme_bw()

## ggplot(probabilities1) +
##     geom_point(aes(x = probability.ML, y = probability.AS_median)) +
##     geom_abline(aes(intercept = 0, slope = 1), alpha = 0.5) +
##     lims(x = c(0,1), y = c(0,1)) +
##     facet_grid(~ category) +
##     theme_bw()


