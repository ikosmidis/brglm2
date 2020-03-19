context("post fit check for infinite estimates")

## endometrial data from Heinze \& Schemper (2002) (see ?endometrial)
data("endometrial", package = "brglm2")
expect_warning({
    endometrial_ml <- glm(HG ~ I(-NV) + PI + EH, data = endometrial,
                          family = binomial("cloglog"))
})

expect_warning(cie <- check_infinite_estimates(endometrial_ml))

test_that("infinte estimates have been found as expected", {
    expect_true(any(cie[, "I(-NV)"] > 1e+06))
})
