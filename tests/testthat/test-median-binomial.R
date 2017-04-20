context("agreement with mbrglm when estimating binomial response models")

data("lizards", package = "brglm2")
data("endo",package = "mbrglm")

links <- lapply(c("logit", "probit", "cloglog", "cauchit"), make.link)


tol <- 1e-10
for (l in seq_along(links)) {
    lizardsMBRlegacy <- mbrglm::mbrglm(cbind(grahami, opalinus) ~ height + diameter +
                                      light + time, family = binomial(links[[l]]), data=lizards,
                                    method = "mbrglm.fit", 
                                    control.mbrglm = mbrglm::mbrglm.control(mbr.epsilon = 1e-10, mbr.maxit = 1000))
    lizardsMBR <- glm(cbind(grahami, opalinus) ~ height + diameter +
                       light + time, family = binomial(links[[l]]), data=lizards,
                     method = "brglmFit", epsilon = 1e-10, maxit = 1000,type="AS_median")
    
    endoMBRlegacy <- mbrglm::mbrglm(HG~NV+PI+EH, family = binomial(links[[l]]), data=endo,method = "mbrglm.fit", 
                                    control.mbrglm = mbrglm::mbrglm.control(mbr.epsilon = 1e-10, mbr.maxit = 1000))
    c1 <- coef(summary(endoMBRlegacy))
    endoMBR <- glm(HG~NV+PI+EH, family = binomial(links[[l]]), data=endo,
                     method = "brglmFit", epsilon = 1e-10, maxit = 1000,type="AS_median")
    c2 <- coef(summary(endoMBR))
    test_that(paste("glm with brglm.fit method  and mbrglm return the same coefficients for the lizards when link is", links[[l]]$name), {
      expect_equal(coef(lizardsMBR), coef(lizardsMBRlegacy), tolerance = tol)
      expect_equal(c1,c2, tolerance = tol)
    })
}

## Performance comparisons BR versus ML
## link1 <- "cauchit"
## system.time(replicate(100, {lizardsBR <- glm(cbind(grahami, opalinus) ~ height + diameter + light + time, family = binomial(link1), data=lizards, method = "brglmFit", epsilon = 1e-10, maxit = 1000, type = "ML")}))

## system.time(replicate(100, {lizardsBR <- glm(cbind(grahami, opalinus) ~ height + diameter + light + time, family = binomial(link1), data=lizards, method = "glm.fit", epsilon = 1e-10, maxit = 1000)}))

## link1 <- "cauchit"
## system.time(replicate(100, {lizardsBR <- glm(cbind(grahami, opalinus) ~ height + diameter + light + time, family = binomial(link1), data=lizards, method = "brglmFit", epsilon = 1e-10, maxit = 1000)}))

## system.time(replicate(100, {lizardsBR <- brglm(cbind(grahami, opalinus) ~ height + diameter + light + time, family = binomial(link1), data=lizards, br.epsilon = 1e-10, br.maxit = 1000)}))
