library("brglm2")

data("lizards", package = "brglm2")
liz_fm <- cbind(grahami, opalinus) ~ height + diameter + light + time
## ML fit = MDYPL fit with `alpha = 1`
liz_ml <- glm(liz_fm, family = binomial(), data = lizards,
              method = "mdyplFit", alpha = 1)
liz_ml0 <- glm(liz_fm, family = binomial(), data = lizards)
## liz_ml is the same fit as liz_ml0
summ_liz_ml <- summary(liz_ml)
summ_liz_ml0 <- summary(liz_ml0)
all.equal(coef(summ_liz_ml), coef(summ_liz_ml0))

to_check <- names(summ_liz_ml)
to_check <- to_check[!(to_check %in% c("hd_correction", "alpha", "type", "call"))]
for (nam in to_check) {
    expect_equal(summ_liz_ml[nam], summ_liz_ml0[nam])
}




liz_mdypl <- glm(liz_fm, family = binomial(), data = lizards,
                 method = "mdyplFit", alpha = 0.8761)
summ_mdypl <- summary(liz_mdypl)
summ_mdypl_c <- summary(liz_mdypl, hd_correction = TRUE)

eta_sloe <- sloe(liz_mdypl)
nn <- sum(weights(liz_mdypl))
p <- length(coef(liz_mdypl)) - 1
se_pars <- solve_se(p / nn, eta_sloe, liz_mdypl$alpha, intercept = coef(liz_mdypl)["(Intercept)"], start = c(0.5, 1, 1, unname(coef(liz_mdypl)["(Intercept)"])), corrupted = TRUE)
expect_equal(se_pars, summ_mdypl_c$se_parameters)

cc <- coef(liz_mdypl)
cc["(Intercept)"] <- se_pars[4]
cc[-1] <- cc[-1] / se_pars[1]
expect_equal(cc, coef(summ_mdypl_c)[, "Estimate"])

ttt <- brglm2:::taus(liz_mdypl)
ses_c <- c(NA, se_pars[3] / (sqrt(nn) * ttt * se_pars[1]))
expect_equal(unname(coef(summ_mdypl_c)[, "Std. Error"]), ses_c)
expect_equal(coef(summ_mdypl_c)[, "z value"], cc / ses_c)

tots <- lizards$grahami + lizards$opalinus
lizX <- model.matrix(liz_mdypl)
expect_equal(summ_mdypl_c$aic,
             brglm2:::logist_aic(liz_mdypl$y_adj, weights(liz_ml), drop(plogis(lizX %*% cc)), weights(liz_ml)) + 2 * liz_mdypl$rank)
