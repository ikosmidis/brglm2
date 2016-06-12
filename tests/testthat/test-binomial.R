context("agreement with brglm when estimation binomial resposne models")

data("lizards")

links <- lapply(c("logit", "probit", "cloglog", "cauchit"), make.link)


tol <- 1e-10
for (l in seq_along(links)) {
    lizardsBRlegacy <- brglm0(cbind(grahami, opalinus) ~ height + diameter +
                                        light + time, family = binomial(links[[l]]), data=lizards,
                                    method = "brglm.fit", br.epsilon = 1e-10, br.maxit = 1000)

    lizardsBR <- glm(cbind(grahami, opalinus) ~ height + diameter +
                         light + time, family = binomial(links[[l]]), data=lizards,
                     method = "brglmFit", epsilon = 1e-10, maxit = 1000)
    test_that(paste("glm with brglm.fit method and brglm_0 return the same coefficients for the lizards when link is", links[[l]]$name), {
        expect_equal(coef(lizardsBR), coef(lizardsBRlegacy), tolerance = tol)
    })
}
