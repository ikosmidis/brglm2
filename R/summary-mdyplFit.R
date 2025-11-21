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
    inds <- is.infinite(S) | is.na(S)
    sd(S[!inds])
}

taus <- function(X) {
    X <- X[, colnames(X) != "(Intercept)"]
    R <- qr.R(qr(X))
    rss <- 1 / colSums(backsolve(R, diag(ncol(X)), transpose = TRUE)^2)
    sqrt(rss / (nrow(X) - ncol(X) + 1))
}

#' Summary method for [`"mdyplFit"`][mdyplFit()] objects
#'
#' @inheritParams stats::summary.glm
#' @inheritParams solve_se
#' @param hd_correction if `FALSE` (default), then the corresponding
#'     quantities are computed according to standard asymptotics. If
#'     `TRUE` then the high-dimensionality corrections in Sterzinger &
#'     Kosmidis (2024) are employed to updates estimates, estimated
#'     standard errors, z-statistics, etc. See Details.
#' @param solve_se_control a list of further arguments to be passed to
#'     [solve_se()]. Even if explicitly specified, the arguments
#'     `kappa`, `ss`, `alpha`, and `intercept` are always set
#'     according to `object`, and `corrupted` is set to `TRUE`.
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
#' strength, which is the limit \eqn{\gamma^2} of \eqn{var(X \beta)}
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
#' \donttest{
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
summary.mdyplFit <- function(object, hd_correction = FALSE,
                             solve_se_control = list(), ...) {
    ## Get summary object
    summ <- summary.glm(object, ...)
    if (isTRUE(hd_correction)) {
        coefs <- coef(summ)
        nobs <- sum(pw <- weights(object))
        has_intercept <- attr(terms(object), "intercept")
        p <- nrow(coefs) - has_intercept
        nu_sloe <- sloe(object)
        theta <- if (has_intercept) coefs["(Intercept)", "Estimate"] else NULL
        solve_se_control$kappa <- ka <- p / nobs
        solve_se_control$ss <- nu_sloe
        solve_se_control$alpha <- object$alpha
        solve_se_control$intercept <- theta
        if (is.null(solve_se_control$start)) {
            solve_se_control$start <- c(0.5, 1, 1, theta)
        }
        solve_se_control$corrupted <- TRUE
        se_pars <- try(do.call("solve_se", solve_se_control), silent = TRUE)
        if (inherits(se_pars, "try-error")) {
            msg <- paste(se_pars, "Unable to solve the state evolution equations. See `?solve_se` and use `solve_se_control` to control the optimizer, including supplying a vector of", 3 + has_intercept, " starting values for `mu` (in (0, 1)), `b` (> 0), `sigma` (> 0)", if (has_intercept) ", `intercept.`" else ".")
            stop(msg)
        }
        xx <- model.matrix(object)[, !summ$aliased]
        tt <- taus(xx)
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
        mus <- family$linkinv(drop(xx %*% coefs[, "Estimate"]))
        y <- object$y
        ## Null deviance is not updated
        d_res <- sqrt(pmax(family$dev.resids(y, mus, pw), 0))
        summ$deviance.resid <- ifelse(y > mus, d_res, -d_res)
        summ$deviance <- sum(summ$deviance.resid^2)
        summ$aic <- logist_aic(object$y_adj, object$n_init, mus, pw, summ$deviance) + 2 * object$rank
        summ$cov.scaled <- summ$cov.unscaled <- NULL
        summ$se_parameters <- se_pars
        if (!isTRUE(all(abs(attr(se_pars, "funcs")) < 1e-04))) {
            msg <- paste("Unable to solve the state evolution equations. See `?solve_se` and use `solve_se_control` to control the optimizer, including supplying a vector of", 3 + has_intercept, " starting values for `mu` (in (0, 1)), `b` (> 0), `sigma` (> 0)", if (has_intercept) ", `intercept.`" else ".")
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
        cat("\nEstimated signal strength (gamma^2) =", round(x$signal_strength, 2))
        cat("\nState evolution parameters (mu, b, sigma) =", paste0("(", paste(round(x$se_parameters[1:3], 2), collapse = ", "), ")"), "with max(|funcs|) =", max(abs(attr(x$se_parameters, "funcs"))), "\n")
    }
    invisible(x)
}
