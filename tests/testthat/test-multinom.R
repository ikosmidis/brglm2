context("tests for multinomBR and its agreement with other methods")

source(system.file("inst", "brpr/brpr.R", package = "brglm2"))

library("pmlr")
library("nnet")


## coefinds takes as input the glm object that brpr results and
## returns the location of the coefficients for the fitted multinomial
## model.
coefinds <- function(obj) {
  require(gnm)
  cats <- as.character(obj$call$formula[3])
  cats <- strsplit(gsub(" ", "", cats), "\\*\\(")[[1]][1]
  cats <- rev(strsplit(cats, "\\+")[[1]])[1]
  cc <- pickCoef(obj, cats)
  cc
}

## expandMult will expand the multinomial data set so that all
## available categories are represented in each covariate class. If
## there are non-available categories in the original data set they
## are included in the resulted data frame and are assigned frequency
## zero.
expandMult <- function(data, catvar, countsvar) {
  require(gnm)
  temp <- expandCategorical(data = data, catvar = catvar,
                            sep = ".", countvar = "CcCounts",
                            idvar = "id", as.ordered = FALSE,
                            group = TRUE)
  ## Get the correct category counts
  temp[[countsvar]] <- temp[[countsvar]] * temp$CcCounts
  temp$CcCounts <-   temp$id <- NULL
  temp
}


#####################################################################
## Analysis of the hepatitis data set in Bull et. al. (2007)
#####################################################################
data("hepatitis", package = "pmlr")
## Construct a variable with the multinomial categories according to
## the HCV and nonABC columns
hepat <- hepatitis
hepat$type <- with(hepat, factor(1 - HCV * nonABC + HCV + 2 * nonABC))
hepat$type <- factor(hepat$type, labels = c("noDisease", "C", "nonABC"))
contrasts(hepat$type) <- contr.treatment(3, base = 1)
## compare the result with the one that pmlr gives
hepatnew <- expandMult(data = hepat, catvar = "type", countsvar = "counts")
contrasts(hepatnew$type) <- contr.treatment(3, base = 1)
heppmlr <- pmlr(type ~ group * time,
                data = hepat, weights = counts, method = "wald",
                penalized = TRUE)
hepbrpr <- brpr(counts ~ type * (group * time),
                fixed.totals = group:time, data = hepatnew)
hepbrglm <- brglmFit(model.matrix(~ group*time + type*(group*time), data = hepatnew), hepatnew$counts, family = poisson(), fixedTotals = with(hepatnew, group:time))

coef(hepglm)[coefinds(hepbrpr)]
coef(hepbrpr)[coefinds(hepbrpr)]
heppmlr$coefficients

parNames <- names(coefinds(hepbrpr))
tol <- 1e-08
test_that("estimates from multinomial regression using brglmFit are the same as those coming form brpr", {
    expect_equal(coef(hepbrglm)[parNames], coef(hepbrpr)[parNames], tol = tol)
})
test_that("coefficient standard errors from multinomial regression using brglmFit are the same as those coming form brpr", {
    expect_equal(sqrt(diag(tcrossprod(solve(qr.R(hepbrpr$qr)))))[parNames], coef(summary(hepbrpr))[parNames, 2], tol = tol)
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
## Re-express data set in an appropriate form for brpr
enzymes$counts <- rep(1, nrow(enzymes))
enzymes$ind <- factor(seq.int(nrow(enzymes)))
enz <- expandMult(enzymes, "Group", "counts")
enz$ind <- factor(enz$ind)
enzpmlr <- pmlr(Group ~ AST + GLDH, weights = counts,
                data = enzymes, method = "wald")
enzbrpr <- brpr(counts ~ -1 + ind + Group * (AST + GLDH),
                 fixed.totals = ind, data = enz,
                epsilon = 1e-14)
enzbrglm <- brglmFit(x = model.matrix(counts ~ -1 + ind + Group * (AST + GLDH), data = enz), y = enz$counts,
                     fixedTotals = enz$ind, family = poisson(), control = list(epsilon = 1e-14))

## brpr appears slower than pmlr in this case, most probably because
## of the large number of nuisances.
coef(enzbrpr)[names(coefinds(enzbrpr))]
coef(enzbrglm)[names(coefinds(enzbrpr))]
coef(enzpmlr)


enzmultinom <- multinom(Group ~ AST + GLDH, weights = counts,
                    data = enzymes)
enzglm <- glm(counts ~ -1 + ind + Group * (AST + GLDH), data = enz, family =poisson())



## Artificial example to check the effect of ignoring to explicitly
## pass the covariate classes with zero counts, when a category is not
## observed.

library("nnet")

## All covariate class - covariate combinations present
set.seed(123)
N <- 50
ncat <- 5
eps <- 0.01
if (N %% ncat == 0) {
    set.seed(123)
    probs <- runif(ncat, eps, ncat - eps)
    probs <- probs/sum(probs)
    totals <- 10
    counts <- c(rmultinom(N, totals, probs))
    resp <- factor(rep(seq.int(ncat), N))
    ## Simulate covariate values and repeat each value ncat times
    x1 <- rep(rnorm(N/ncat, 2, 2), each = ncat)
    x2 <- rep(rexp(N/ncat, 10), each = ncat)
    ## Indicator of covariate classes
    inds <- factor(rep(seq.int(N), each = ncat))
    ## Collect everything in a data frame
    artificialData0 <- data.frame(resp = resp, x1 = x1, x2 = x2, counts = counts, inds = inds)
    ## Clean up
    rm(list = c("counts", "resp", "x1", "x2", "inds"))

    ## Model fitting
    fit_multinom0 <- multinom(resp ~ x1 + x2, weights = counts, data = artificialData0)
    fit_glm0 <- glm(counts ~ -1 + inds + resp * (x1 + x2), data = artificialData0, family = poisson())

    ## Effect of removing a few observations
    set.seed(222)
    del <- sample(seq.int(nrow(artificialData0)), 5)
    artificialData1 <- artificialData0[-del, ]
    ## Model fitting
    fit_multinom1 <- multinom(resp ~ x1 + x2, weights = counts, data = artificialData1)
    fit_glm1 <- glm(counts ~ -1 + inds + resp * (x1 + x2), data = artificialData1, family = poisson())

    ## Now set zero counts for the missing category - covariate class combinations and refit
    artificialData1a <- within(artificialData0, {
        counts[del] <- 0
    })
    fit_multinom1a <- multinom(resp ~ x1 + x2, weights = counts, data = artificialData1a)
    fit_glm1a <- glm(counts ~ -1 + inds + resp * (x1 + x2), data = artificialData1a, family = poisson())
}



## Coefficients form glm
report <- function(fit1, fit2, tol = 1e-03) {
    if (inherits(fit1, "glm") & inherits(fit2, "glm")) {
        coefs1 <- coef(fit1)
        coefs2 <- coef(fit2)
    }
    if (inherits(fit1, "glm") & inherits(fit2, "multinom")) {
        coefs2 <- as.matrix(coef(fit2))
        q <- nrow(coefs2)
        parNames <- c(paste0("resp", 2:(q+1)), paste0("resp", 2:(q+1), ":x1"), paste0("resp", 2:(q+1), ":x2"))
        coefs1 <- matrix(coef(fit1)[parNames], nrow = nrow(coefs2))
        dimnames(coefs1) <- dimnames(coefs2)
    }
    if (inherits(fit1, "multinom") & inherits(fit2, "glm")) {
        coefs1 <- coef(fit1)
        q <- nrow(coefs1)
        parNames <- c(paste0("resp", 2:(q+1)), paste0("resp", 2:(q+1), ":x1"), paste0("resp", 2:(q+1), ":x2"))
        coefs2 <- matrix(coef(fit2)[parNames], nrow = nrow(coefs1))
        dimnames(coefs2) <- dimnames(coefs1)
        out <- list(multinom1 = coefs1, multinom2 = coefs2, comparison = all.equal(coefs1, coefs2, tol = tol))
    }
    if (inherits(fit1, "multinom") & inherits(fit2, "multinom")) {
        coefs1 <- coef(fit1)
        coefs2 <- coef(fit2)
        out <- list(multinom1 = coefs1, multinom2 = coefs2, comparison = all.equal(coefs1, coefs2, tol = tol))
    }
    ## if (inherits(fit1, "glm")) {
    ##     parNames <- c(paste0("resp", 2:(q+1)), paste0("resp", 2:(q+1), ":x1"), paste0("resp", 2:(q+1), ":x2"))
    ##     coefs1 <- matrix(coef(fit1)[parNames], nrow = nrow(coefs2))
    ##     dimnames(coefs1) <- dimnames(coefs2)
    ##     out <- list(glm = coefs1, multinom = coefs2, comparison = all.equal(coefs1, coefs2, tol = tol))
    ## }
    ## if (inherits(fit1, "multinom")) {
    ##     coefs1 <- coef(fit1)
    ##     out <- list(multinom1 = coefs1, multinom2 = coefs2, comparison = all.equal(coefs1, coefs2, tol = tol))
    ## }
    list(fit1 = coefs1,
         fit2 = coefs2,
         comparison = all.equal(coefs1, coefs2, tol = tol))
}

## fit_glm0 and fit_multinom0 are fits based on artificialData0 which
## has all covariate class - response combinations
report(fit_glm0, fit_multinom0)

## fit_glm1 and fit_multinom1 are fits after removing the observations
## in del from artificialData0
report(fit_glm1, fit_multinom1)

## fit_glm1a and fit_multinom1a are fits after adding the covariate
## class - response combinations with zero counts
report(fit_glm1a, fit_multinom1a)

## fit_multinom1a and fit_multinom1 should be the same
report(fit_multinom1, fit_multinom1a)
