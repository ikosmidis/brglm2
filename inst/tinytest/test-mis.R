## Dobson (1990) Page 93: Randomized Controlled Trial :

logit_mis <- mis(link = "logit", sensitivity = 1, specificity = 1)
probit_mis <- mis(link = "probit", sensitivity = 1, specificity = 1)

expect_identical(class(logit_mis), "link-glm")
expect_identical(class(probit_mis), "link-glm")

lizards_f <- cbind(grahami, opalinus) ~ height + diameter + light + time
lizards_logit <- glm(lizards_f, family = binomial(logit), data = lizards)
lizards_probit <- glm(lizards_f, family = binomial(probit), data = lizards)

lizards_logit_mis <- update(lizards_logit, family = binomial(logit_mis),
                            start = coef(lizards_logit))
lizards_probit_mis <- update(lizards_probit, family = binomial(probit_mis),
                             start = coef(lizards_probit))

## mis link with sensitivity and specificity 1 are the same as original links
expect_equal(coef(lizards_logit), coef(lizards_logit_mis), tolerance = 1e-06)
expect_equal(coef(lizards_probit), coef(lizards_probit_mis), tolerance = 1e-06)

