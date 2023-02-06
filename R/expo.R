# Copyright (C) 2022- Ioannis Kosmidis

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

#' Estimate the exponential of parameters of generalized linear models
#' using various methods
#'
#' The [expo()] method uses the supplied [`"brglmFit"`][brglmFit] or
#' [`"glm"`][glm] object to estimate the exponential of parameters of
#' generalized linear models with maximum likelihood or various mean
#' and median bias reduction methods. [expo()] is useful for computing
#' (corrected) estimates of the multiplicative impact of a unit
#' increase on a covariate on the mean of a Poisson log-linear model
#' (`family = poisson("log")` in [glm()]) while adjusting for other
#' covariates, the odds ratio associated with a unit increase on a
#' covariate in a logistic regression model (`family =
#' binomial("logit")` in [glm()]) while adjusting for other
#' covariates, the relative risk associated with a unit increase on a
#' covariate in a relative risk regression model (`family =
#' binomial("log")` in [glm()]) while adjusting for other covariates,
#' among others.
#'
#' @aliases brglmFit_expo
#' @aliases expo
#' @param object an object of class [`"brglmFit"`][brglmFit] or
#'     [`"glm"`][glm].
#' @param type the type of correction to be used. The available
#'     options are `"correction*"` (explicit mean bias correction with
#'     a multiplicative adjustment), `"correction*"` (explicit mean
#'     bias correction with an additive adjustment), `"Lylesetal2012"`
#'     (explicit median bias correction using the multiplicative
#'     adjustment in Lyles et al., 2012), `"AS_median"` (median bias
#'     reduction), and `"ML"` (maximum likelihood). See Details.
#' @param level the confidence level required. Default is `0.95`.
#'
#' @details
#'
#' The supported methods through the `type` argument are:
#'
#' * `"ML"`: the estimates of the exponentiated parameters are
#' \eqn{\exp(\hat\theta_j)}, where \eqn{\theta_j} is the maximum
#' likelihood estimates for the \eqn{j}th regression parameter.
#'
#' * `"correction*"`: the estimates of the exponentiated parameters
#' are \eqn{\exp(\hat\theta_j) / (1 + \hat{v}_j / 2)}, where
#' \eqn{\hat\theta_j} is the estimate of the \eqn{j}th regression
#' parameter using `type = "AS_mixed"` in [brglmFit()].
#'
#' * `"correction+"`: the estimates of the exponentiated parameters
#' are \eqn{\exp(\hat\theta_j) (1 - \hat{v}_j / 2)}, where
#' \eqn{\hat\theta_j} is the estimate of the \eqn{j}th regression
#' parameter using `type = "AS_mixed"` in [brglmFit()].
#'
#' * `"Lylesetal2012"`: the estimates of the exponentiated parameters
#' are \eqn{\exp(\hat\theta_j) exp(- \hat{v}_j / 2)}, where
#' \eqn{\hat\theta_j} is the estimate of the \eqn{j}th regression
#' parameter using `type = "AS_mixed"` in [brglmFit()]. This estimator
#' has been proposed in Lyles et al. (2012).
#'
#' * `"AS_median"`: the estimates of the exponentiated parameters are
#' \eqn{\exp(\hat\theta_j)}, where \eqn{\hat\theta_j} is the estimate
#' of the \eqn{j}th regression parameter using `type = "AS_median"` in
#' [brglmFit()].
#'
#' `"correction*"` and `"correction+"` are based on multiplicative and
#' additive adjustments, respectively, of the exponential of a
#' reduced-bias estimator (like the ones coming from [brglmFit()] with
#' `type = "AS_mixed"`, `type = "AS_mean"`, and `type =
#' "correction"`). The form of those adjustments results from the
#' expression of the first-term in the mean bias expansion of the
#' exponential of a reduced-bias estimator. See, for example, Di
#' Caterina & Kosmidis (2019, expression 12) for the general form of
#' the first-term of the mean bias of a smooth transformation of a
#' reduced-bias estimator.
#'
#' The estimators from `"correction+"`, `"correction*"`,
#' `"Lylesetal2012"` have asymptotic mean bias of order smaller than
#' than of the maximum likelihood estimator. The estimators from
#' `"AS_median"` are asymptotically closed to being median unbiased
#' than the maximum likelihood estimator is.
#'
#' Estimated standard errors are computed using the delta method,
#' where both the Jacobin and the information matrix are evaluated at
#' the logarithm of the estimates of the exponentiated parameters.
#'
#' Confidence intervals results by taking the exponential of the
#' limits of standard Wald-type intervals computed at the logarithm of
#' the estimates of the exponentiated parameters.
#'
#'
#' @return a list inheriting from class [`"brglmFit_expo"`][brglmFit_expo] with
#'     components `coef` (the estimates of the exponentiated
#'     regression parameters), `se` (the corresponding estimated
#'     standard errors for the exponentiated parameters), `ci`
#'     (confidence intervals of level `level` for the exponentiated
#'     parameters), and `type` for the `type` of correction that has
#'     been requested.
#'
#' @author Ioannis Kosmidis `[aut, cre]` \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso [brglm_fit()] and and [brglm_control()]
#'
#' @references
#'
#' Di Caterina C, Kosmidis I (2019). Location-Adjusted Wald Statistics for Scalar
#' Parameters. *Computational Statistics & Data Analysis*, **138**,
#' 126-142. \doi{10.1016/j.csda.2019.04.004}.
#'
#' Kosmidis I, Kenne Pagui E C, Sartori N (2020). Mean and median bias
#' reduction in generalized linear models. *Statistics and Computing*,
#' **30**, 43-59. \doi{10.1007/s11222-019-09860-6}.
#'
#' Cordeiro G M, McCullagh P (1991). Bias correction in generalized
#' linear models. *Journal of the Royal Statistical Society. Series B
#' (Methodological)*, **53**, 629-643. \doi{10.1111/j.2517-6161.1991.tb01852.x}.
#'
#' Lyles R H, Guo Y, Greenland S (2012). Reducing bias and mean
#' squared error associated with regression-based odds ratio
#' estimators. *Journal of Statistical Planning and Inference*,
#' **142** 3235–3241. \doi{10.1016/j.jspi.2012.05.005}.
#'
#' @examples
#'
#' ## The lizards example from ?brglm::brglm
#' lizardsML <- glm(cbind(grahami, opalinus) ~ height + diameter +
#'                  light + time, family = binomial(logit), data = lizards,
#'                  method = "glm.fit")
#' # Get estimates, standard errors, and confidence intervals of odds
#' # ratios with various methods
#' expo(lizardsML, type = "ML")
#' expo(lizardsML, type = "correction*")
#' expo(lizardsML, type = "Lylesetal2012")
#' expo(lizardsML, type = "correction+")
#' expo(lizardsML, type = "AS_median")
#'
#' ## Example from ?glm
#' ## Dobson (1990) Page 93: Randomized Controlled Trial :
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
#' expo(glm.D93, type = "correction*")
#'
#' @export
expo.brglmFit <- function(object, type = c("correction*", "correction+", "Lylesetal2012", "AS_median", "ML"), level = 0.95) {
    type <- match.arg(type)
    to_correct <- type %in% c("correction*", "correction+", "Lylesetal2012")
    if (to_correct) {
        if (isTRUE(object$type != "AS_mixed")) {
            object <- update(object, type = "AS_mixed", data = object$data)
        }
    } else {
        ## AS_mean and ML are invariant to transformation
        if (isTRUE(object$type != type)) {
            object <- update(object, type = type, data = object$data)
        }
    }

    info_fun <- get_information_function(object)
    dispersion <- object$dispersion
    coefs <- coef(object, model = "mean")
    ses <- sqrt(diag(solve(info_fun(coefs, dispersion))))[seq_along(coefs)]
    trans_coefs <- exp(coefs) * switch(type,
                                       "ML" = 1,
                                       "AS_median" = 1,
                                       "Lylesetal2012" = exp(- ses^2/2),
                                       "correction*" = 1 / (1 + ses^2/2),
                                       "correction+" = (1 - ses^2/2))

    if (to_correct) {
        ## Recompute coefs and ses on the original scale from trans_coefs if method is not invariant
        coefs <- log(trans_coefs)
        ses <- sqrt(diag(solve(info_fun(coefs, dispersion))))[seq_along(coefs)]
    }

    trans_ses <- trans_coefs * ses
    a <- (1 - level)/2
    ci <- coefs + qnorm(a) * cbind(ses, -ses)
    colnames(ci) <- paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE, digits = 3), "%")
    out <- list(coef = trans_coefs, se = trans_ses, ci = exp(ci), type = type)
    out$call <- match.call()
    out$family <- object$family
    class(out) <- c("brglmFit_expo", class(out))
    out
}

#' @rdname expo.brglmFit
#' @export
expo.glm <- function(object, type = c("correction*", "correction+", "Lylesetal2012", "AS_median", "ML"), level = 0.95) {
    type <- match.arg(type)
    fit_type <- switch(type,
                       "correction*" = "AS_mixed",
                       "correction+" = "AS_mixed",
                       "Lylesetal2012" = "AS_mixed",
                       "AS_median" = "AS_median",
                       "ML" = "ML")
    object <- update(object, method = brglmFit, type = fit_type)
    out <- expo(object, type, level)
    out$call <- match.call()
    out
}


#' @method print brglmFit_expo
#' @export
print.brglmFit_expo <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cmat <- cbind(x$coef, x$se, x$ci)
    fam_link <- paste0(x$family$family, x$family$link)
    interpr <- switch(fam_link,
                      "binomiallogit" = "Odds ratios",
                      "poissonlog" = "Multiplicative effects to the mean",
                      "Gammalog" = "Multiplicative effects to the mean",
                      "binomiallog" = "Relative risks",
                      "inverse.gaussianlog"= "Multiplicative effects to the mean",
                      NA)
    colnames(cmat) <- c("Estimate", "Std. Error", colnames(x$ci))
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    if (!is.na(interpr)) {
        cat(interpr, "\n")
    }
    printCoefmat(cmat, digits = digits, ...)
    cat("\n\nType of estimator:", x$type, get_type_description(x$type), "\n")
}

#' Extract estimates from [`"brglmFit_expo"`][brglmFit_expo] objects
#'
#' @inheritParams stats::coef
#'
#' @method coef brglmFit_expo
#' @export
coef.brglmFit_expo <- function(object, ...) {
    object$coef
}
