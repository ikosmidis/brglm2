## simulation study: infant birth weights

library("MASS")
library("brglm2")

## Prepare the birth weight data set
bwt <- with(birthwt, {
    age <- age
    racewhite <- ifelse(race==1,1,0)
    smoke <- smoke
    ptl <- ifelse(ptl>0,1,0)
    ptd <- factor(ptl > 0)
    ht <- ht
    loglwt <- log(lwt)
    data.frame(normwt = 1-low, age, racewhite, smoke, ptl,ht,loglwt,ftv)
})
bwt <- subset(bwt, subset = (ftv==0), select= -c(ftv))

## true models
bwt_ml <- glm(normwt ~ ., family = binomial, data = bwt)
beta <- coef(bwt_ml)

## random seed
set.seed(435)

## simulated response variable
Nsim <- 10000
simdata <- simulate(bwt_ml, nsim = Nsim)

ml <- br <- bc <- mbr <- ml_se <- br_se <- bc_se <- mbr_se <- separation <-
    matrix(NA, nrow = Nsim, ncol = length(beta))

for (i in 1:Nsim) {
    current_data <- within(bwt, { normwt <- simdata[[i]] })
    if (i%%100 == 0) print(i)
    ml_fit <- update(bwt_ml, data = current_data)
    br_fit <- update(bwt_ml, method = "brglmFit", type = "AS_mean", data = current_data)
    bc_fit <- update(bwt_ml, method = "brglmFit", type = "correction", data = current_data)
    mbr_fit <- update(bwt_ml, method = "brglmFit", type = "AS_median", data = current_data)
    sep_fit <- update(ml_fit, method = "detect_separation")

    sum_ml <- summary(ml_fit)
    sum_brmean <- summary(br_fit)
    sum_bcorr <- summary(bc_fit)
    sum_brmedian <- summary(mbr_fit)

    ml[i, ] <- sum_ml$coef[, 1]
    ml_se[i, ] <- sum_ml$coef[, 2]
    br[i, ] <- sum_brmean$coef[, 1]
    br_se[i, ] <- sum_brmean$coef[, 2]
    bc[i, ] <- sum_bcorr$coef[, 1]
    bc_se[i, ] <- sum_bcorr$coef[, 2]
    mbr[i, ] <- sum_brmedian$coef[, 1]
    mbr_se[i, ] <- sum_brmedian$coef[, 2]
    separation[i, ] <- sep_fit$betas
}

ml.inc <- apply(separation, 1, function(b) all(b == 0))

## bias in beta parameterization
bias.beta <- data.frame(ml = colMeans(ml[ml.inc, ]),
                        br = colMeans(br),
                        mbr = colMeans(mbr)) - beta
## bias in psi parameterization
bias.psi <- data.frame(ml = colMeans(exp(ml[ml.inc, ])),
                       br = colMeans(exp(br)),
                       mbr = colMeans(exp(mbr))) - exp(beta)
## Coverage of 95% Wald confidence intervals
ml.ci.l <- ml[ml.inc, ] - qnorm(0.975) * ml_se[ml.inc, ]
ml.ci.u <- ml[ml.inc, ] + qnorm(0.975) * ml_se[ml.inc, ]
br.ci.l <- br - qnorm(0.975) * br_se
br.ci.u <- br + qnorm(0.975) * br_se
mbr.ci.l <- mbr - qnorm(0.975) * mbr_se
mbr.ci.u <- mbr + qnorm(0.975) * mbr_se
coverage <- data.frame(ml = rowMeans(t(ml.ci.l) < beta &  t(ml.ci.u) > beta),
                       br = rowMeans(t(br.ci.l) < beta &  t(br.ci.u) > beta),
                       mbr = rowMeans(t(mbr.ci.l) < beta &  t(mbr.ci.u) > beta))
## probability of underestimation
PU <- data.frame(ml = rowMeans(t(ml[ml.inc, ]) < beta),
                 br = rowMeans(t(br) < beta),
                 mbr = rowMeans(t(mbr) < beta))

save(bias.beta, bias.psi, PU, coverage, file = "birth_weight_simulation_results.rda")
