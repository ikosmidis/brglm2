# Copyright (C) 2025- Ioannis Kosmidis

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

#' Fitting function for [glm()] for maximum Diaconis-Ylvisaker prior
#' penalized likelihood estimation of logistic regression models
#'
#' [mdyplFit()] is a fitting method for [glm()] that fits logistic
#' regression models using maximum Diaconis-Ylvisaker prior penalized
#' likelihood estimation.
#'
#' @inheritParams stats::glm.fit
#' @aliases mdypl_fit
#' @param x a design matrix of dimension `n * p`.
#' @param y a vector of observations of length `n`.
#' @param control a list of parameters controlling the fitting
#'     process. See [mdyplControl()] for details.
#'
#' @details
#'
#' [mdyplFit()] uses [stats::glm.fit()] to fit a logistic regression
#' model on responses `alpha * y + (1 - alpha) / 2`, where `y` are the
#' original binomial responses scaled by the binomial totals. This is
#' equivalent to penalizing the likelihood by the Diaconis-Ylvisaker
#' prior with shrinkage parameter \eqn{\alpha} and regression parameters
#' set to zero. See Rigon & Aliverti (2023) and Sterzinger & Kosmidis
#' (2024).
#'
#' By default, `alpha = n / (p + n)` is used, where `n` is the sum of
#' the binomial totals. Alternative values of `alpha` can be passed to
#' the `control` argument; see [mdyplControl()] for setting up the
#' list passed to `control`. If `alpha = 1` then [mdyplFit()] will
#' simply do maximum likelihood estimation.
#'
#' Note that `null.deviance`, `deviance` and `aic` in the resulting
#' object are computed at the adjusted responses. Hence, methods such
#' as [logLik()][stats::logLik()] and [AIC()][stats::AIC()] use the
#' penalized log-likelihood. With the default `alpha`, the inferential
#' procedures based on penalized likelihood are asymptotically
#' equivalent to the ones that use the unpenalized likelihood when
#' `p/n` is vanishing asymptotically.
#'
#' For high-dimensionality corrected estimates, standard errors and z
#' statistics, use the [`summary`][summary.mdyplFit()] method for
#' [`"mdyplFit"`][mdyplFit()] objects with `hd_correction = TRUE`.
#'
#' [mdypl_fit()] is an alias to [mdyplFit()].
#'
#' @return
#'
#' An object inheriting from [`"mdyplFit"`][mdyplFit()] object, which
#' is a list having the same elements to the list that
#' [stats::glm.fit()] returns, with a few extra arguments.
#'
#' @author Ioannis Kosmidis `[aut, cre]` \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso [mdyplControl()], [summary.mdyplFit()], [plrtest.mdyplFit()], [glm()]
#'
#' @references
#'
#' Sterzinger P, Kosmidis I (2024). Diaconis-Ylvisaker prior
#' penalized likelihood for \eqn{p/n \to \kappa \in (0,1)} logistic
#' regression. *arXiv*:2311.07419v2, \url{https://arxiv.org/abs/2311.07419}.
#'
#' Rigon T, Aliverti E (2023). Conjugate priors and bias reduction for
#' logistic regression models. *Statistics & Probability Letters*,
#' **202**, 109901. \doi{10.1016/j.spl.2023.109901}.
#'
#' @examples
#'
#' data("lizards", package = "brglm2")
#' liz_fm <- cbind(grahami, opalinus) ~ height + diameter + light + time
#' ## ML fit = MDYPL fit with `alpha = 1`
#' liz_ml <- glm(liz_fm, family = binomial(), data = lizards,
#'               method = "mdyplFit", alpha = 1)
#' liz_ml0 <- glm(liz_fm, family = binomial(), data = lizards)
#'
#' ## liz_ml is the same fit as liz_ml0
#' summ_liz_ml <- summary(liz_ml)
#' summ_liz_ml0 <- summary(liz_ml0)
#' all.equal(coef(summ_liz_ml), coef(summ_liz_ml0))
#'
#' ## MDYPL fit with default `alpha` (see `?mdyplControl`)
#' liz_fm <- cbind(grahami, opalinus) ~ height + diameter + light + time
#' liz_mdypl <- glm(liz_ml, family = binomial(), data = lizards,
#'                  method = "mdyplFit")
#'
#' ## Comparing outputs from ML and MDYPL, with and without
#' ## high-dimensionality corrections.
#' summary(liz_mdypl)
#' summary(liz_mdypl, hd_correction = TRUE)
#' summ_liz_ml
#' summary(liz_ml, hd_correction = TRUE)
#' ## Not much difference in fits here as this is a low dimensional
#' ## problem with dimensionality constant
#' (liz_ml$rank - 1) / sum(weights(liz_ml))
#'
#'
#'
#' ## The case study in Section 8 of Sterzinger and
#' ## Kosmidis (2024)
#' data("MultipleFeatures", package = "brglm2")
#'
#' ## Center the fou.* and kar.* features
#' vars <- grep("fou|kar", names(MultipleFeatures), value = TRUE)
#' train_id <- which(MultipleFeatures$training)
#' MultipleFeatures[train_id, vars] <- scale(MultipleFeatures[train_id, vars], scale = FALSE)
#' ## Compute the MDYPL fits
#' kappa <- length(vars) / sum(MultipleFeatures$training)
#' full_fm <- formula(paste("I(digit == 7) ~", paste(vars, collapse = " + ")))
#' nest_vars <- grep("fou", vars, value = TRUE)
#' nest_fm <- formula(paste("I(digit == 7) ~", paste(nest_vars, collapse = " + ")))
#' full_m <- glm(full_fm, data = MultipleFeatures, family = binomial(),
#'               method = mdyplFit, alpha = 1 / (1 + kappa), subset = training)
#' nest_m <- update(full_m, nest_fm)
#'
#' ## With a naive penalized likelihood ratio test we get no evidence
#' ## against the hypothesis that the model with only `fou` features
#' ## is an as good descrition of `7` as the model with both `fou` and
#' ## `kar` features.
#' plrtest(nest_m, full_m)
#'
#' ## With a high-dimensionality correction theres is strong evidence
#' ## against the model with only `fou` features
#' plrtest(nest_m, full_m, hd_correction = TRUE)
#'
#'
#' \dontrun{
#' ## A simulated data set as in Rigon & Aliverti (2023, Section 4.3)
#'
#' set.seed(123)
#' n <- 1000
#' p <- 500
#' gamma <- sqrt(5)
#' X <- matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
#' betas0 <- rep(c(-1, -1/2, 0, 2, 3), each = p / 5)
#' betas <- gamma * betas0 / sqrt(sum(betas0^2))
#' probs <- plogis(drop(X %*% betas))
#' y <- rbinom(n, 1, probs)
#' fit_mdypl <- glm(y ~ -1 + X, family = binomial(), method = "mdyplFit")
#'
#' ## The default value of `alpha` is `n / (n + p)` here
#' identical(n / (n + p), fit_mdypl$alpha)
#'
#' ## Aggregate bias of MDYPL and rescaled MDYPL estimators
#' ag_bias <- function(estimates, beta) mean(estimates - beta)
#' ag_bias(coef(summary(fit_mdypl))[, "Estimate"], betas)
#' ag_bias(coef(summary(fit_mdypl, hd_correction = TRUE))[, "Estimate"], betas)
#'
#' }
#' @export
mdyplFit <- function(x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
                     mustart = NULL, offset = rep(0, nobs), family = binomial(),
                     control = list(), intercept = TRUE,
                     singular.ok = TRUE) {

    nobs <- NROW(y)
    if (!isTRUE(family$family == "binomial" && family$link == "logit")) {
        stop('`mdyplFit` currently supports only `binomial` family with `"logit"` link')
    }

    control <- do.call("mdyplControl", control)

    missing_offset <- is.null(offset)

    if (is.null(weights)) {
        weights <- rep.int(1, nobs)
    }

    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }

    if (missing_offset) {
        offset <- rep.int(0, nobs)
    }

    alpha <- unless_null(control$alpha, sum(weights) / (sum(weights) + ncol(x) - intercept))

    ## adjust responses as per MDYPL with beta_P = 0
    y_adj <- alpha * y + (1 - alpha) / 2
    control_glm <- glm.control(epsilon = control$epsilon,
                               maxit = control$maxit, trace = control$trace)
    out <- glm.fit(x = x, y = y_adj, weights = weights,
                   etastart = etastart, mustart = mustart,
                   offset = offset, family = quasibinomial(),
                   control = control_glm,
                   intercept = intercept, singular.ok = singular.ok)

    mus <- out$fitted.values
    if (intercept & missing_offset) {
        nullmus <- mdyplFit(x = x[, "(Intercept)", drop = FALSE], y = y, weights = weights,
                            offset = rep(0, nobs), family = family, intercept = TRUE,
                            control = control,
                            start = family$linkfun(mean(y)))$fitted.values
    }

    if (!intercept) {
        nullmus <- family$linkinv(offset)
    }
    ## If there is an intercept and an offset then, for calculating
    ## the null deviance glm will make a call to the fitter to fit the
    ## glm with intercept and the offset
    if (intercept & !missing_offset) {
        nullmus <- mus
        ## doen't really matter what nullmus is set to. glm will make
        ## a new call to mdyplFit and use the deviance from that call
        ## as null
    }

    out$family <- family

    ## Reset quantities in terms of original responses wherever needed
    dev.resids <- family$dev.resids
    out$null.deviance <- sum(dev.resids(y_adj, nullmus, weights))
    out$deviance <- sum(dev.resids(y_adj, mus, weights))
    out$aic <- logist_aic(y_adj, n, mus, weights, deviance) + 2 * out$rank
    out$residuals <- (y - mus) / (mus * (1 - mus))
    out$y_adj <- y_adj
    out$y <- y
    out$alpha <- alpha
    out$type <- "MPL_DY"
    out$control <- control
    out$class <- c("mdyplFit")
    out$n_init <- n ## needed when `hd_correction = TRUE` in summary where aic is recomputed
    out
}

## Similar to binomial()$aic but works with y in (0, 1)
dbinom2 <- function(x, size, prob, log = FALSE) {
    su <- x
    fa <- size - x
    db <- prob^su * (1 - prob)^fa / beta(su + 1, fa + 1) / (size + 1)
    db[size < su] <- 0
    if (isTRUE(log)) log(db) else db
}

logist_aic <- function(y, n, mu, wt, dev) {
    m <- if (any(n > 1)) n else wt
    -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom2(m * y, m, mu, log = TRUE))
}

#' Auxiliary function for [glm()] fitting using the [brglmFit()]
#' method.
#'
#' Typically only used internally by [brglmFit()], but may be used to
#' construct a `control` argument.
#'
#' @aliases mdypl_control
#' @param alpha the shrinkage parameter (in `[0, 1]`) in the
#'     Diaconis-Ylvisaker prior penalty. Default is \code{NULL}, which
#'     results in `alpha = n / (n + p)`, where `n` is the sum of the
#'     binomial totals and `p` is the number of model
#'     parameters. Setting `alpha = 1` corresponds to using maximum
#'     likelihood, i.e. no penalization. See Details.
#' @param epsilon positive convergence tolerance epsilon. Default is
#'     `1e-08`.
#' @param maxit integer giving the maximal number of iterations
#'     allowed. Default is `25`.
#' @param trace logical indicating if output should be produced for
#'     each iteration. Default is `FALSE`.
#'
#' @details
#'
#' Internally, [mdyplFit()] uses [stats::glm.fit()] to fit a logistic
#' regression model on responses `alpha * y + (1 - alpha) / 2`, where
#' `y` are the original binomial responses scaled by the binomial
#' totals. `epsilon`, `maxit` and `trace` control the
#' [stats::glm.fit()] call; see [stats::glm.control()].
#'
#' @return
#'
#' A list with components named as the arguments.
#'
#' @author Ioannis Kosmidis `[aut, cre]` \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso [mdyplFit()], [glm.control()]
#'
#' @export
mdyplControl <- function(alpha = NULL, epsilon = 1e-08, maxit = 25, trace = FALSE) {
    out <- glm.control(epsilon, maxit, trace)
    if (!is.null(alpha)) {
        if (!(is.numeric(alpha)) || isTRUE(alpha < 0) || isTRUE(alpha > 1))
            stop("`alpha` should be in [0, 1]")
    }
    out$alpha <- alpha
    out
}


#' Method for computing confidence intervals for one or more
#' regression parameters in a [`"mdyplFit"`][mdyplFit()] object
#'
#' @inheritParams stats::confint
#' @inheritParams summary.mdyplFit
#'
#' @author Ioannis Kosmidis `[aut, cre]` \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso [mdyplFit()], [summary.mdyplFit()]
#'
#' @examples
#'
#' \dontrun{
#'
#' set.seed(123)
#' n <- 2000
#' p <- 800
#' set.seed(123)
#' betas <- c(rnorm(p / 4, mean = 7, sd = 1), rep(0, 3 * p / 4))
#' X <- matrix(rnorm(n * p, 0, 1/sqrt(n)), nrow = n, ncol = p)
#' probs <- plogis(drop(X %*% betas))
#' y <- rbinom(n, 1, probs)
#' fit_mdypl <- glm(y ~ -1 + X, family = binomial(), method = "mdyplFit")
#'
#' wald_ci <- confint(fit_mdypl)
#' adj_wald_ci <- confint(fit_mdypl, hd_correction = TRUE)
#' ag_coverage <- function(cis, beta) mean((cis[, 1] < beta) & (cis[, 2] > beta))
#' ag_coverage(wald_ci, betas)
#' ag_coverage(adj_wald_ci, betas)
#'
#' }
#'
#' @method confint mdyplFit
#' @export
confint.mdyplFit <- function(object, parm, level = 0.95, hd_correction = FALSE, ...) {
    ## A modification of confint.default to use summary objects
    summ <- summary(object, hd_correction = hd_correction, ...)
    coefs <- coef(summ)
    cf <- coefs[, "Estimate"]
    pnames <- rownames(coefs)
    if (missing(parm))
        parm <- pnames
    else if (is.numeric(parm))
        parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    ## Directly from stats:::.format_perc
    pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
    ses <- coefs[, "Std. Error"]
    ci[] <- cf[parm] + ses %o% fac
    ci
}

