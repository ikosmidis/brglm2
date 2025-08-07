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
#' orginal binomial responses scaled by the binomial totals. This is
#' equivalent to penalizing the likelihood by the Diaconis-Ylvisaker
#' prior with shirnkage parameter $\alpha$ and regression parameters
#' set to zero. See Rigon & Aliverti (2023) and Sterzinger & Kosmidis
#' (2024).
#'
#' [mdypl_fit()] is an alias to [mdyplFit()].
#'
#' @return
#'
#' An object inheriting from [`"mdyplFit"`][mdyplFit()] object, which
#' is a list having the same elements to the list that
#' [stats::glm.fit()] returns, with a few extra arguments. By default,
#' `alpha = m / (p + m)` is used, where `m` is the sum of the binomial
#' totals. Alternative values of `alpha` can be passed to the
#' `control` argument; see [mdyplControl()] for setting up the list
#' passed to `control`.
#'
#' @author Ioannis Kosmidis `[aut, cre]` \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso [mdyPLcontrol()], [glm.fit()], [glm()]
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
#' ## A simulated data set as in Rigon & Aliverti (2023, Section 4.3)
#'
#' set.seed(123)
#' n <- 1000
#' p <- 200
#' gamma <- 5
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
#' cols <- hcl.colors(3, alpha = 0.2)
#' par(mfrow = c(1, 2))
#' plot(betas, type = "l", ylim = c(-1, 1),
#'      main = "MDYPL estimates",
#'      xlab = "Parameter index", ylab = NA)
#' points(coef(fit_mdypl), col = NA, bg = cols[1], pch = 21)
#' sc_betas <- hd_summary.mdyplFit(fit_mdypl, se_start = c(0.5, 1, 1))
#' plot(betas, type = "l", ylim = c(-1, 1),
#'      main = "rescaled MDYPL estimates",
#'      xlab = "Parameter index", ylab = NA)
#' points(sc_betas[, "Rescaled-estimate"], col = NA, bg = cols[2], pch = 21)
#'
#' @export
mdyplFit <- function(x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
                     mustart = NULL, offset = rep(0, nobs), family = binomial(),
                     control = list(), intercept = TRUE,
                     fixed_totals = NULL, singular.ok = TRUE) {

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

    alpha <- unless_null(control$alpha, sum(weights) / (sum(weights) + ncol(x)))

    ## adjust responses as per MDYPL with beta_P = 0
    y_adj <- alpha * y + (1 - alpha) / 2

    out <- glm.fit(x = x, y = y_adj, weights = weights,
                   etastart = etastart, mustart = mustart,
                   offset = offset, family = quasibinomial(),
                   control = list(epsilon = control$epsilon,
                                  maxit = control$maxit, trace = control$trace),
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
    out$null.deviance <- sum(dev.resids(y, nullmus, weights))
    out$y_adj <- y_adj
    out$y <- y
    out$deviance <- sum(dev.resids(y, mus, weights))
    out$aic <- family$aic(y, n, mus, weights, deviance) + 2 * out$rank
    out$alpha <- alpha
    out$type <- "MPL_DY"

    out$class <- c("mdyplFit", "brglmFit")
    out
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
#' `y` are the orginal binomial responses scaled by the binomial
#' totals. `epsilon`, `maxit` and `trace` control the
#' [stats::glm.fit()] call; see [stats::glm.control()].
#'
#' @return
#'
#' A list with components named as the arguments.
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
#' In partiuclar, [sloe()] computes an estimate of the corrupteed
#' signal strength which is the limit \deqn{\nu^2} of \eqn{var(X
#' \hat\beta(\alpha))}, where \eqn{\hat\beta(\alpha)} is the maximimum
#' Diaconis-Ylvisaker prior penalized likelihood (MDYPL) estimator as
#' computed by [mdyplFit()] with shirnkage parameter \eqn{alpha}.
#'
#' @return
#'
#' A scalar.
#'
#' @references
#'
#' Sterzinger P, Kosmidis I (2024). Diaconis-Ylvisaker prior
#' penalized likelihood for \eqn{p/n \to \kappa \in (0,1)} logistic
#' regression. *arXiv*:2311.07419v2, \url{https://arxiv.org/abs/2311.07419}.
#'
#' Yadlowsky S, Yun T, McLean CY, D' Amour A (2021). SLOE: A Faster
#' Method for Statistical Inference in High-Dimensional Logistic
#' Regression. In M Ranzato, A Beygelzimer, Y Dauphin, P Liang, JW
#' Vaughan (eds.), *Advances in Neural Information Processing
#' Systems*, **34**, 29517–29528. Curran Associates, Inc. \url{https://proceedings.neurips.cc/paper_files/paper/2021/file/f6c2a0c4b566bc99d596e58638e342b0-Paper.pdf}.
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
    RSS <- 1 / colSums(backsolve(L, diag(ncol(X)), transpose = TRUE)^2)
    sqrt(RSS / (nrow(X) - ncol(X) + 1))
}

#' High-dimensionality correction of estimates and Wald statistics
#' from [mdyplFit()] objects.
#'
#' @param object an [`"mdyplFit"`][mdyplFit()] object.
#' @param se_start starting values for the parameters of the state
#'     evoultion equations; see [solve_se()].
#' @param null a vector of null values for the parameters estimated in
#'     `object`. Default is `0`.
#'
#' @return
#'
#' A table of corrected coefficients, Wald statistics and p-values
#'
#' @export
hd_summary.mdyplFit <- function(object, se_start, null = 0, ...) {
    coefs <- coef(object)
    nobs <- sum(object$prior.weights)
    has_intercept <- attr(terms(object), "intercept")
    p <- length(coefs) - has_intercept
    eta_sloe <- sloe(object)
    se_pars <- solve_se(kappa = p / nobs, ss = eta_sloe, alpha = object$alpha,
                        intercept = if (has_intercept) coef(object)["(Intercept)"] else NULL,
                        start = se_start,
                        corrupted = TRUE,
                        ...)
    tt <- taus(object)
    adj_z <- sqrt(nobs) * tt * (coef(object) - se_pars[1] * null) / se_pars[3]
    adj_coef <- coefs / se_pars[1]
    pv <- 2 * pnorm(-abs(adj_z))
    coef_table <- cbind(adj_coef, adj_z, pv)
    dimnames(coef_table) <- list(names(coefs), c("Rescaled-estimate", "z value", "Pr(>|z|)"))
    coef_table
}
