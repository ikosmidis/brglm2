context("test that detect_separation works as expected")

## endometrial data from Heinze \& Schemper (2002) (see ?endometrial)
data("endometrial", package = "brglm2")
endometrial_separation <- glm(HG ~ I(-NV) + PI + EH, data = endometrial,
                              family = binomial("cloglog"),
                              method = "detect_separation")

## The lizards example from ?brglm::brglm
data("lizards", package = "brglm2")
lizards_separation <- glm(cbind(grahami, opalinus) ~ height + diameter +
                              light + time, family = binomial(logit), data = lizards,
                          method = "detect_separation")


test_that("infinte estimates have been found as expected", {
    expect_equal(endometrial_separation$betas, c(0, -Inf, 0, 0), check.attributes = FALSE)
    expect_equal(lizards_separation$betas, numeric(6), check.attributes = FALSE)
})


## ## hepatitis
## data("hepatitis", package = "pmlr")
## hepat <- hepatitis
## hepat$type <- with(hepat, factor(1 - HCV * nonABC + HCV + 2 * nonABC))
## hepat$type <- factor(hepat$type, labels = c("noDisease", "C", "nonABC"))
## y <- rowSums(hepat$counts*nnet::class.ind(hepat$type)[,c(1, 3)])
## glm(y/counts ~ group * time, data = hepat, weights = counts, family = binomial(),
##     method = "detect_separation")
## dd <- data.frame(y = c(1,1,1,0,0), x = c(1,2,4,4,5), off = c(1,2, 2,4,3))
## summary(glm(y ~ x + offset(off), family = binomial(), data = dd))
