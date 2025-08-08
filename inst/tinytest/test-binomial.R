library("brglm")
data("lizards", package = "brglm2")
links <- lapply(c("logit", "probit", "cloglog", "cauchit"), make.link)

tol <- 1e-08
for (l in seq_along(links)) {
    lizardsBRlegacy <- brglm(cbind(grahami, opalinus) ~ height + diameter +
                                 light + time, family = binomial(links[[l]]),
                             data = lizards, method = "brglm.fit",
                             br.epsilon = 1e-10, br.maxit = 1000)

    lizardsBR <- glm(cbind(grahami, opalinus) ~ height + diameter +
                         light + time, family = binomial(links[[l]]), data=lizards,
                     method = "brglmFit", epsilon = 1e-10, maxit = 1000)
    ## glm with brglm.fit method and brglm_0 return the same coefficients for the lizards when link is links[[l]]$name)
    expect_equal(coef(lizardsBR), coef(lizardsBRlegacy), tolerance = tol)
}


## Performance comparisons BR versus ML
## link1 <- "cauchit"
## system.time(replicate(100, {lizardsBR <- glm(cbind(grahami, opalinus) ~ height + diameter + light + time, family = binomial(link1), data=lizards, method = "brglmFit", epsilon = 1e-10, maxit = 1000, type = "ML")}))

## system.time(replicate(100, {lizardsBR <- glm(cbind(grahami, opalinus) ~ height + diameter + light + time, family = binomial(link1), data=lizards, method = "glm.fit", epsilon = 1e-10, maxit = 1000)}))

## link1 <- "cauchit"
## system.time(replicate(100, {lizardsBR <- glm(cbind(grahami, opalinus) ~ height + diameter + light + time, family = binomial(link1), data=lizards, method = "brglmFit", epsilon = 1e-10, maxit = 1000)}))

## system.time(replicate(100, {lizardsBR <- brglm(cbind(grahami, opalinus) ~ height + diameter + light + time, family = binomial(link1), data=lizards, br.epsilon = 1e-10, br.maxit = 1000)}))
