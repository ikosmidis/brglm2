# Copyright (C) 2021 Ioannis Kosmidis

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/



#' @rdname ordinal_superiority
#' @export
ordinal_superiority.bracl <- function(object, formula, data,
                                      measure = c("gamma", "Delta"),
                                      level = 0.95,
                                      bc = FALSE) {
    measure <- match.arg(measure)
    ## If bc is TRUE and object is not a reduced mean-bias fit then
    ## compute reduced mean-bias estimators.
    if (!inherits(object, "bracl")) {
        stop("ordinal superiority measures are not available for objects of class ", class(object)[1])
    }
    if (!isTRUE(object$parallel)) {
        stop("ordinal superiority measures are available only for \"bracl\" objects with `parallel = TRUE`")
    }
    if (isTRUE(bc) & !(object$type %in% c("AS_mean", "AS_mixed"))) {
        object <- update(object, type = "AS_mean")
    }
    mf <- model.frame(formula, model.frame(object))
    Terms <- attr(mf, "terms")
    z <- model.matrix(Terms, mf, object$contrasts[all.vars(formula)])
    znames <- colnames(z)[-match("(Intercept)", colnames(z), nomatch = 0L)]
    if (isTRUE(length(znames) > 1) | !isTRUE(all(sort(unique(z[, znames])) %in% c(0, 1)))) {
        stop("`formula` can have only one grouping explanatory variable with two levels")
    }
    X <- model.matrix(object, data = data)
    z_ind <- match(znames, colnames(X), nomatch = 0L)
    if (isTRUE(z_ind == 0)) {
        stop("the grouping explanatory variable must be one of the explanatory variables in `object`")
    }
    Xnoz <- unique(X[, -z_ind, drop = FALSE])
    X <- matrix(NA, nrow = nrow(Xnoz), ncol = ncol(X), dimnames = list(NULL, colnames(X)))
    X[, -z_ind] <- Xnoz
    nx <- nrow(X)
    gammas <- .ordsup(coef(object),
                      X = X, group_id = z_ind, ncat = object$ncat,
                      ref = object$ref, lev = object$lev, measure = "gamma")
    coef_vcov <- vcov(object)
    ## mean bias reduction of gammas
    bias_gammas <- numeric(length(gammas))
    if (bc) {
        for (j in seq.int(nx)) {
            hess <- numDeriv::hessian(function(theta, ...) .ordsup(theta, ...)[j],
                                      x = coef(object),
                                      X = X, group_id = z_ind, ncat = object$ncat,
                                      ref = object$ref, lev = object$lev, measure = "gamma")
            bias_gammas[j] <- sum(diag(coef_vcov %*% hess)) / 2
        }
    }
    gammas <- gammas - bias_gammas
    ## mean_gammas <- mean(gammas)
    ## compute standard error for gamma
    grads <- numDeriv::jacobian(.ordsup, x = coef(object),
                                X = X, group_id = z_ind, ncat = object$ncat,
                                ref = object$ref, lev = object$lev, measure = "gamma")
    ## grad_mean <- apply(grads, 2, mean)
    se <- apply(grads, 1, function(x) sqrt(crossprod(x, (coef_vcov %*% x))))
    ## se_mean <- sqrt(crossprod(grad_mean, (coef_vcov %*% grad_mean)))
    ## Confidence intervals as in Agresti and Kateri
    a <- 1/2 + level/2
    pct <- paste(format(100 * c(1 - a, a), trim = TRUE, scientific = FALSE, digits = 3), "%")
    lsd <- drop(qnorm(a) * se / (gammas * (1 - gammas)))
    ci <- qlogis(gammas) + cbind(rep(-1, nx),  1) * lsd
    ## lsd_mean <- drop(qnorm(a) * se_mean / (mean_gammas * (1 - mean_gammas)))
    ## ci_mean <- qlogis(mean_gammas) + c(-1, 1) * lsd_mean
    if (isTRUE(measure == "Delta")) {
        out <- cbind(Xnoz, 2 * gammas - 1, 2 * se, 2 * plogis(ci) - 1)
        ## out_mean <- c(2 * mean_gammas - 1, 2 * se_mean, 2 * plogis(ci_mean) - 1)
        colnames(out)[ncol(Xnoz) + 1:4] <- c("Delta", "se", pct)
        ## names(out_mean) <- c("Delta*", "se", pct)
    }
    else {
        out <- cbind(Xnoz, gammas, se, plogis(ci))
        ## out_mean <- c(mean_gammas, se_mean, plogis(ci_mean))
        colnames(out)[ncol(Xnoz) + 1:4] <- c("gamma", "se", pct)
        ## names(out_mean) <- c("gamma*", "se", pct)
    }
    ## list(individual = out, mean = out_mean)
    out
}

## X should be the covariate values where the ordsup is computed and a column for z. X is duplicated internally with z == 1 and z == 0
## only for parallel = TRUE
## group_varianble: indicator of x
.ordsup <- function(coef, X, group_id, ncat, ref, lev, measure = "gamma") {
    X <- rbind(X, X)
    nX <- nrow(X)
    X[, group_id] <- z <- rep(c(0, 1), each = nX / 2)
    nams <- names(coef)
    int <- (ncat - 1):1
    sl <- nams[-int]
    coefs <- cbind(rev(cumsum(coef[int])),
                   int * matrix(coef[sl], nrow = ncat - 1, ncol = length(sl), byrow = TRUE))
    rownames(coefs) <- lev[-ref]
    fits <- matrix(0, nrow = nrow(X), ncol = ncat, dimnames = list(rownames(X), lev))
    fits1 <- apply(coefs, 1, function(b) X %*% b)
    fits[, rownames(coefs)] <- fits1
    Y <- t(apply(fits, 1, function(x) exp(x) / sum(exp(x))))
    probs0 <- Y[z == 0, , drop = FALSE]
    probs1 <- Y[z == 1, , drop = FALSE]
    gamma_fun <- function(p0, p1) {
        out <- outer(p0, p1, "*")
        sum(out[upper.tri(out)]) + sum(diag(out)) / 2
    }
    out <- numeric(nX / 2)
    for (i in seq.int(nX / 2)) {
        out[i] <- gamma_fun(probs0[i, ], probs1[i, ])
    }
    if (measure == "gamma") out else 2 * out - 1
}
