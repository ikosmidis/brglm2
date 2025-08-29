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
#' By default, `alpha = m / (p + m)` is used, where `m` is the sum of
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

    ## Reset quantities in terms of original responses
    dev.resids <- family$dev.resids
    out$null.deviance <- sum(dev.resids(y_adj, nullmus, weights))
    out$deviance <- sum(dev.resids(y_adj, mus, weights))
    out$aic <- logist_aic(y_adj, n, mus, weights, deviance) + 2 * out$rank
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
#'     results in `alpha = m / (m + p)`, where `m` is the sum of the
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

#' Estimate the corrupted signal strength in a model with
#' (sub-)Gaussian covariates
#'
#' @param object an [`"mdyplFit"`][mdyplFit()] object.
#'
#'
#' @details
#'
#' The Signal Strength Leave-One-Out Estimator (SLOE) is defined in
#' Yadlowsky et al. (2021) when the model is estimated using maximum
#' likelihood (i.e. when `object$alpha = 1`; see [mdyplControl()] for
#' what `alpha` is). The SLOE adaptation when estimation is through
#' maximum Diaconis-Ylvisaker prior penalized likelihood
#' ([mdypl_fit()]) has been put forward in Sterzinger & Kosmidis
#' (2025).
#'
#' In particular, [sloe()] computes an estimate of the corrupted
#' signal strength which is the limit \deqn{\nu^2} of \eqn{var(X
#' \hat\beta(\alpha))}, where \eqn{\hat\beta(\alpha)} is the maximum
#' Diaconis-Ylvisaker prior penalized likelihood (MDYPL) estimator as
#' computed by [mdyplFit()] with shrinkage parameter \eqn{alpha}.
#'
#' @return
#'
#' A scalar.
#'
#' @author Ioannis Kosmidis `[aut, cre]` \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso [summary.mdyplFit()]
#'
#' @references
#'
#' Sterzinger P, Kosmidis I (2024). Diaconis-Ylvisaker prior
#' penalized likelihood for \eqn{p/n \to \kappa \in (0,1)} logistic
#' regression. *arXiv*:2311.07419v2, \url{https://arxiv.org/abs/2311.07419}.
#'
#' Yadlowsky S, Yun T, McLean C Y, D' Amour A (2021). SLOE: A Faster
#' Method for Statistical Inference in High-Dimensional Logistic
#' Regression. In M Ranzato, A Beygelzimer, Y Dauphin, P Liang, J W
#' Vaughan (eds.), *Advances in Neural Information Processing
#' Systems*, **34**, 29517–29528. Curran Associates,
#' Inc. \url{https://proceedings.neurips.cc/paper_files/paper/2021/file/f6c2a0c4b566bc99d596e58638e342b0-Paper.pdf}.
#'
#' @export
sloe <- function(object) {
    mu <- fitted(object)
    v <- mu * (1 - mu)
    h <- hatvalues(object)
    S <- object$linear.predictors - (object$y_adj - mu) / v * h / (1 - h)
    sd(S)
}

taus <- function(object) {
    X <- model.matrix(object)
    has_intercept <- attr(terms(object), "intercept")
    if (has_intercept) X <- X[, colnames(X) != "(Intercept)"]
    L <- qr.R(qr(X))
    rss <- 1 / colSums(backsolve(L, diag(ncol(X)), transpose = TRUE)^2)
    sqrt(rss / (nrow(X) - ncol(X) + 1))
}


#' Summary method for [`"mdyplFit"`][mdyplFit()] objects
#'
#' @inheritParams stats::summary.glm
#' @inheritParams solve_se
#' @param hd_correction if `FALSE` (default), then the summary
#'     corresponding to standard asymptotics is computed. If `TRUE`
#'     then the high-dimensionality corrections in Sterzinger &
#'     Kosmidis (2024) are employed to updates estimates, estimated
#'     standard errors, z-statistics. See Details.
#' @param se_start a vector of starting values for the state evolution
#'     equations. See the `start` argument in [solve_se()].
#' @param ... further arguments to be passed to [summary.glm()].
#'
#' @details
#'
#' If `hd_correction = TRUE`, the [sloe()] estimator of the square
#' root of the corrupted signal strength is estimated from `object`,
#' as are the conditional variances of each covariate given the others
#' (excluding the intercept). The latter are estimated using residual
#' sums of squares from the linear regression of each covariate on all
#' the others, as proposed in Zhao et al (2021, Section 5.1). Then the
#' appropriate state evolution equations are solved using [solve_se()]
#' with `corrupted = TRUE`, and the obtained constants are used to
#' rescale the estimates, and adjust estimated standard errors and
#' z-statistics as in Sterzinger & Kosmidis (2024).
#'
#' The key assumptions under which the rescaled estimates and corrected
#' standard errors and z-statistics are asymptotically valid are that
#' the covariates have sub-Gaussian distributions, and that the signal
#' strength, which is the limit \deqn{\gamma^2} of \eqn{var(X \beta)}
#' is finite as \eqn{p / n \to \kappa \in (0, 1)}, with \eqn{\kappa \in
#' (0, 1)}. See Sterzinger & Kosmidis (2024).
#'
#' If `hd_correction = TRUE`, and the model has an intercept, then the
#' result provides only a corrected estimate of the intercept with no
#' accompanying standard error, z-statistic, and p-value. Also,
#' `vcov(summary(object, hd_correction = TRUE))` is always
#' `NULL`. Populating those objects with appropriate estimates is the
#' subject of current work.
#'
#' @return
#'
#' A list with objects as in the result of [stats::summary.glm()],
#' with extra component `se_parameters`, which is the vector of the
#' solution to the state evolution equations with extra attributes
#' (see [solve_se()]).
#'
#' @author Ioannis Kosmidis `[aut, cre]` \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso [mdyplFit()], [solve_se()]
#'
#' @references
#'
#' Zhao Q, Sur P, Cand\`es E J (2022). The asymptotic distribution of
#' the MLE in high-dimensional logistic models: Arbitrary
#' covariance. *Bernoulli*, **28**,
#' 1835–1861. \doi{10.3150/21-BEJ1401}.
#'
#' Sterzinger P, Kosmidis I (2024). Diaconis-Ylvisaker prior
#' penalized likelihood for \eqn{p/n \to \kappa \in (0,1)} logistic
#' regression. *arXiv*:2311.07419v2, \url{https://arxiv.org/abs/2311.07419}.
#'
#' @examples
#'
#' \dontrun{
#'
#' set.seed(123)
#' n <- 2000
#' p <- 400
#' set.seed(123)
#' betas <- c(rnorm(p / 2, mean = 7, sd = 1), rep(0, p / 2))
#' X <- matrix(rnorm(n * p, 0, 1/sqrt(n)), nrow = n, ncol = p)
#' probs <- plogis(drop(X %*% betas))
#' y <- rbinom(n, 1, probs)
#' fit_mdypl <- glm(y ~ -1 + X, family = binomial(), method = "mdyplFit")
#'
#' st_summary <- summary(fit_mdypl)
#' hd_summary <- summary(fit_mdypl, hd_correction = TRUE)
#'
#' cols <- hcl.colors(3, alpha = 0.2)
#' par(mfrow = c(1, 2))
#' plot(betas, type = "l", ylim = c(-3, 14),
#'      main = "MDYPL estimates",
#'      xlab = "Parameter index", ylab = NA)
#' points(coef(st_summary)[, "Estimate"], col = NA, bg = cols[1], pch = 21)
#'
#' plot(betas, type = "l", ylim = c(-3, 14),
#'      main = "rescaled MDYPL estimates",
#'      xlab = "Parameter index", ylab = NA)
#' points(coef(hd_summary)[, "Estimate"], col = NA, bg = cols[2], pch = 21)
#'
#' ## z-statistics
#' z_mdypl <- coef(st_summary)[betas == 0, "z value"]
#' qqnorm(z_mdypl, col = NA, bg = cols[1], pch = 21, main = "z value")
#' abline(0, 1, lty = 2)
#' z_c_mdypl <- coef(hd_summary)[betas == 0, "z value"]
#' qqnorm(z_c_mdypl, col = NA, bg = cols[2], pch = 21, main = "corrected z value")
#' abline(0, 1, lty = 2)
#'
#'}
#'
#' @method summary mdyplFit
#' @export
summary.mdyplFit <- function(object, hd_correction = FALSE, se_start,
                             gh = NULL, prox_tol = 1e-10, transform = TRUE,
                             init_iter = 50,
                             init_method = "Nelder-Mead", ...) {
    ## Get summary object
    summ <- summary.glm(object, ...)
    if (isTRUE(hd_correction)) {
        coefs <- coef(summ)
        nobs <- sum(pw <- weights(object))
        has_intercept <- attr(terms(object), "intercept")
        p <- nrow(coefs) - has_intercept
        nu_sloe <- sloe(object)
        theta <- if (has_intercept) coefs["(Intercept)", "Estimate"] else NULL
        if (missing(se_start)) {
            se_start <- c(0.5, 1, 1, theta)
        }
        ka <- p / nobs
        se_pars <- try(solve_se(kappa = ka, ss = nu_sloe, alpha = object$alpha,
                                intercept = if (has_intercept) theta else NULL,
                                start = se_start,
                                corrupted = TRUE, gh = gh, prox_tol = prox_tol,
                                transform = transform, init_method = init_method,
                                init_iter = init_iter), silent = TRUE)
        if (inherits(se_pars, "try-error")) {
            msg <- paste("Unable to solve the state evolution equations. Try to supply an alternative vector of ", 3 + has_intercept, " to `se_start`, with starting values for `mu` (in (0, 1)), `b` (> 0), `sigma` (> 0)", if (has_intercept) ", `intercept.`" else ".")
            stop(msg)
        }
        tt <- taus(object)
        no_int <- !(rownames(coefs) %in% "(Intercept)")
        coefs[no_int, "Estimate"] <- coefs[no_int, "Estimate"] / se_pars[1]
        coefs[no_int, "Std. Error"] <- se_pars[3] / (sqrt(nobs) * tt * se_pars[1])
        coefs[, "z value"] <- coefs[, "Estimate"] / coefs[, "Std. Error"]
        coefs[, "Pr(>|z|)"] <- 2 * pnorm(-abs(coefs[, "z value"]))
        if (has_intercept) {
            coefs["(Intercept)", ] <- c(se_pars[4], NA, NA, NA)
        }
        summ$coefficients <- coefs

        family <- object$family
        dev.resids <- family$dev.resids
        mus <- family$linkinv(drop(model.matrix(object) %*% coefs[, "Estimate"]))
        y <- object$y
        ## Null deviance is not updated
        d_res <- sqrt(pmax(family$dev.resids(y, mus, pw), 0))
        summ$deviance.resid <- ifelse(y > mus, d_res, -d_res)
        summ$deviance <- sum(summ$deviance.resid^2)
        summ$aic <- logist_aic(object$y_adj, object$n_init, mus, pw, summ$deviance) + 2 * object$rank
        summ$cov.scaled <- summ$cov.unscaled <- NULL
        summ$se_parameters <- se_pars
        if (!isTRUE(all(abs(attr(se_pars, "funcs")) < 1e-04))) {
            msg <- paste0("Potentially unstable solution of the state evolution equations. Try to supply an alternative vector of ", 3 + has_intercept, " to `se_start`, with starting values for `mu` (in (0, 1)), `b` (> 0), `sigma` (> 0)", if (has_intercept) ", `intercept.`" else ".")
            warning(msg)
        }
        summ$signal_strength <- (nu_sloe^2 - ka * se_pars[3]^2) / se_pars[1]^2
        summ$kappa <- ka
    }
    summ$hd_correction <- hd_correction
    summ$alpha <- object$alpha
    summ$type <- object$type
    class(summ) <- c("summary.mdyplFit", class(summ))
    summ
}

## Almost all code in print.summary.mdyplFit is from
## stats:::print.summary.glm apart from minor modifications
#' @rdname summary.mdyplFit
#' @method print summary.mdyplFit
#' @export
print.summary.mdyplFit <- function (x, digits = max(3L, getOption("digits") - 3L),
                                    symbolic.cor = x$symbolic.cor,
                                    signif.stars = getOption("show.signif.stars"), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    cat("Deviance Residuals: \n")
    if (x$df.residual > 5) {
        x$deviance.resid <- setNames(quantile(x$deviance.resid,
                                              na.rm = TRUE), c("Min", "1Q", "Median", "3Q", "Max"))
    }
    xx <- zapsmall(x$deviance.resid, digits + 1L)
    print.default(xx, digits = digits, na.print = "", print.gap = 2L)
    if (length(x$aliased) == 0L) {
        cat("\nNo Coefficients\n")
    } else {
        df <- if ("df" %in% names(x))
                  x[["df"]]
              else NULL
        if (!is.null(df) && (nsingular <- df[3L] - df[1L]))
            cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
                sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients
        if (!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn,
                                                                     colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                     na.print = "NA", ...)
    }
    cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ",
        format(x$dispersion), ")\n\n", apply(cbind(paste(format(c("Null",
                                                                  "Residual"), justify = "right"), "deviance:"), format(unlist(x[c("null.deviance",
                                                                                                                                   "deviance")]), digits = max(5L, digits + 1L)), " on",
                                                   format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"),
                                             1L, paste, collapse = " "), sep = "")
    if (nzchar(mess <- naprint(x$na.action)))
        cat("  (", mess, ")\n", sep = "")
    cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)))
    cat("\n\nType of estimator:", x$type, get_type_description(x$type), "with alpha =", round(x$alpha, 2))
    cat("\n", "Number of Fisher Scoring iterations: ", x$iter, "\n", sep = "")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            } else {
                correl <- format(round(correl, 2L), nsmall = 2L,
                                 digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    if (x$hd_correction) {
        cat("\nHigh-dimensionality correction applied with")
        cat("\nDimentionality parameter (kappa) =", round(x$kappa, 2))
        cat("\nEstimated signal strength (gamma) =", round(x$signal_strength, 2))
        cat("\nState evolution parameters (mu, b, sigma) =", paste0("(", paste(round(x$se_parameters[1:3], 2), collapse = ", "), ")"), "with max(|funcs|) =", max(abs(attr(x$se_parameters, "funcs"))), "\n")
    }
    invisible(x)
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
confint.mdyplFit <- function(object, parm, level = 0.95, hd_correction = FALSE, se_start, ...) {
    ## A modification of confint.default to use summary objects
    summ <- summary(object, hd_correction = hd_correction, se_start, ...)
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
