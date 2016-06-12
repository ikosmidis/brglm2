library("plyr")
library("doMC")
library("ggplot2")
registerDoMC(3)

## Simulated data gamma
simuFun <- function(nsimu, X, beta, disp) {
    set.seed(123)
    ## What is the transformation from shape a, scale s to mean mu dispersion disp? Check out McCullagh & Nelder
    ## mean = a*s and var = a*s^2
    ## mean = mu and var = mu^2*disp
    ## So: a = mu/s and mu*disp = s
    ## So: a = disp and s = mu/disp
    mu <- 1/drop(X%*%beta)
    N <- nrow(X)
    samples <- lapply(seq.int(nsimu), function(i) rgamma(N, shape = mu/disp, scale = disp))
    ## Test that we are doing the right thing here:
    ## all.equal(rowMeans(samples), mu, tol = 0.001)
    res <- ldply(samples, function(y) {
        mML <- glm(y ~ -1 + X, family = Gamma(), method = "glm.fit", start = beta)
        mBC <- glm(y ~ -1 + X, family = Gamma(), method = "brglmFit", start = beta, correction = TRUE)
        mBR <- glm(y ~ -1 + X, family = Gamma(), method = "brglmFit", start = beta)
        c0 <- c(coef(mML), summary(mML)$dispersion)
        c1 <- c(coef(mML), mBR$dispersionML)
        c2 <- c(coef(mBR), mBR$dispersion)
        c3 <- c(coef(mBC), mBC$dispersion)
        out <- data.frame(rbind(c0, c1, c2, c3))
        colnames(out) <- c("beta1", "beta2", "dispersion")
        out$method <- c("MLp", "ML", "BR", "BC")
        out
    }, .parallel = TRUE)
    biases <- ddply(res, .(method), function(dat) {
        colMeans(dat[, c("beta1", "beta2", "dispersion")]) - c(beta, disp)
    })
    sds <- ddply(res, .(method), function(dat) {
        apply(dat[, c("beta1", "beta2", "dispersion")], 2, sd)
    })
    list(biases = biases, sds = sds)
}

set.seed(123)
X <- cbind(1, rexp(10, 3))
ssize <- numeric(10)
res <- as.list(numeric(length(ssize)))
for (j in seq_along(ssize)) {
    res[[j]] <- simuFun(nsimu = 1000, X = X, beta = c(0.5, 1), disp = 0.05)
    ssize[j] <- nrow(X)
    cat("Simulation for sample size:", ssize[j], "is complete\n")
    print(res[[j]])
    X <- rbind(X, X)
}

biases <- data.frame(do.call("rbind", lapply(res, function(x) x$biases)), size = rep(ssize, each = 4))
sds <- data.frame(do.call("rbind", lapply(res, function(x) x$sds)), size = rep(ssize, each = 4))

ggplot(data = biases) + geom_point(aes(x = size, y = abs(beta2))) + facet_grid(~ method) + scale_y_log10() + scale_x_log10()

lm(log(abs(beta2)) ~ log(size), data = biases, subset = method == "ML")
lm(log(abs(beta2)) ~ log(size), data = biases, subset = method == "BR")


## Example from useR!2011 presentation
library(MASS)
clot <- data.frame(u = c(5,10,15,20,30,40,60,80,100),
                   lot1 = c(118,58,42,35,27,25,21,19,18),
                   lot2 = c(69,35,26,21,18,16,13,12,12))
mML <- glm(lot1 ~ log(u), data=clot, family=Gamma())

trueCoefs <- coef(mML)

samples <- simulate(mML, 10000)

pvalues <- adply(samples, 2, function(y) {
    dd <- clot
    dd$lot1 <- unlist(y)
    modML <- glm(lot1 ~ log(u), data = dd, method = "glm.fit", start = trueCoefs, family = Gamma())
    modBRinverse <- glm(lot1 ~ log(u), data = dd, method = "brglmFit", start = trueCoefs, family = Gamma(), dispTrans = "inverse")
    modBRidentity <- glm(lot1 ~ log(u), data = dd, method = "brglmFit", start = trueCoefs, family = Gamma(), dispTrans = "identity")
    modBRsqrt <- glm(lot1 ~ log(u), data = dd, method = "brglmFit", start = trueCoefs, family = Gamma(), dispTrans = "sqrt")
    modBRlog <- glm(lot1 ~ log(u), data = dd, method = "brglmFit", start = trueCoefs, family = Gamma(), dispTrans = "log")
    sdML <- summary(modML, dispersion = MASS::gamma.dispersion(modML))$coef[, 2]
    sdMLdeviance <- summary(modML)$coef[, 2]
    sdBRinverse <- summary(modBRinverse)$coef[, 2]
    sdBRidentity <- summary(modBRidentity)$coef[, 2]
    sdBRlog <- summary(modBRlog)$coef[, 2]
    sdBRsqrt <- summary(modBRsqrt)$coef[, 2]
    ## Wlad statistics for the hypothesis that the coefficients take the true value
    zML <- (coef(modML) - trueCoefs)/sdML
    zMLdeviance <- (coef(modML) - trueCoefs)/sdMLdeviance
    zBRinverse <- (coef(modBRinverse) - trueCoefs)/sdBRinverse
    zBRsqrt <- (coef(modBRsqrt) - trueCoefs)/sdBRsqrt
    zBRidentity <- (coef(modBRidentity) - trueCoefs)/sdBRidentity
    zBRlog <- (coef(modBRlog) - trueCoefs)/sdBRlog
    out <- rbind(2 * pnorm(-abs(zML)),
                 2 * pnorm(-abs(zMLdeviance)),
                 2 * pnorm(-abs(zBRinverse)),
                 2 * pnorm(-abs(zBRsqrt)),
                 2 * pnorm(-abs(zBRidentity)),
                 2 * pnorm(-abs(zBRlog)))
    out <- data.frame(out, method = c("ML", "MLdeviance", "BRinverse", "BRsqrt", "BRidentity", "BRlog"))
    out
}, .parallel = FALSE)

alphas <- c(0.01, 0.05, 0.1)
res <- lapply(alphas, function(a) ddply(pvalues, .(method), function(dat) 1 - mean(dat$log.u. < a)))
names(res) <- alphas
res
