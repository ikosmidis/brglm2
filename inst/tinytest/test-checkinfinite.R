## endometrial data from Heinze \& Schemper (2002) (see ?endometrial)
data("endometrial", package = "brglm2")
expect_warning({
    endometrial_ml <- glm(HG ~ I(-NV) + PI + EH, data = endometrial,
                          family = binomial("cloglog"))
})

expect_error(cie <- check_infinite_estimates(endometrial_ml))
