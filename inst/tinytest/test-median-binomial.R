library("mbrglm")

data("lizards", package = "brglm2")
data("endometrial",package = "brglm2")

links <- lapply(c("logit", "probit", "cloglog", "cauchit"), make.link)

tol <- 1e-10
mbrglmControl <- mbrglm.control(mbr.epsilon = 1e-10, mbr.maxit = 1000)

for (l in seq_along(links)) {
    ## Lizards
    lizardsFormula <- cbind(grahami, opalinus) ~ height + diameter + light + time
    lizardsMBRlegacy <- mbrglm(lizardsFormula, family = binomial(links[[l]]), data = lizards,
                                       method = "mbrglm.fit",
                                       control.mbrglm = mbrglmControl)
    lizardsMBR <- glm(lizardsFormula, family = binomial(links[[l]]), data = lizards,
                      method = "brglmFit", type= "AS_median", epsilon = 1e-10, maxit = 1000)

    ## Endometrial
    endoFormula <- HG ~ NV + PI + EH
    endoMBRlegacy <- mbrglm(endoFormula, family = binomial(links[[l]]), data = endometrial,
                            method = "mbrglm.fit",
                            control.mbrglm = mbrglmControl)
    endoMBR <- glm(endoFormula, family = binomial(links[[l]]), data = endometrial,
                   method = "brglmFit", type = "AS_median", epsilon = 1e-10, maxit = 1000)

    c1 <- coef(summary(endoMBRlegacy))
    c2 <- coef(summary(endoMBR))

    ## glm with brglmFit method and mbrglm return the same coefficients when link is" links[[l]]$name
    expect_equal(coef(lizardsMBR), coef(lizardsMBRlegacy), tolerance = tol)
    expect_equal(c1, c2, tolerance = tol)
}

