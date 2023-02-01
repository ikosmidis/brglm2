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

#' Estimate exponentiated parameters of generalized linear models
#' using various methods
#'
#' [expo()] updates the supplied [brglmFit()]
#' object to estimate exponentiated parameters of generalized linear
#' models with maximum likelihood or various mean and median bias
#' reduction methods.
#'
#' @aliases brglmFit_expo
#' @param object an object of class [brglmFit()],
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
#' The supported methods are:
#'
#' COMPLETE ME
#'
#' @return a list inheriting from class [brglmFit_expo] with
#'     components `coef` (the estimates of the exponentiated
#'     regression parameters), `se` (the corresponding estimated
#'     stadnard errors for the exponentiated parameters), `ci`
#'     (confidence intervals of level `level` for the exponentiated
#'     parameters), and `type` for the `type` of correction that has
#'     been requested.
#'
#' @author Ioannis Kosmidis \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso [brglm_fit()] and and [brglm_control()]
#'
#' @references
#'
#'
#' Kosmidis I, Kenne Pagui E C, Sartori N (2020). Mean and median bias
#' reduction in generalized linear models. *Statistics and Computing*,
#' **30**, 43-59 \doi{10.1007/s11222-019-09860-6}
#'
#' Cordeiro G M, McCullagh P (1991). Bias correction in generalized
#' linear models. *Journal of the Royal Statistical Society. Series B
#' (Methodological)*, **53**, 629-643 \doi{10.1111/j.2517-6161.1991.tb01852.x}
#'
#' Lyles R H, Guo Y, Greenland S (2012). Reducing bias and mean
#' squared error associated with regression-based odds ratio
#' estimators’. *Journal of Statistical Planning and Inference*,
#' **142** 3235–3241. \doi{10.1016/j.jspi.2012.05.005}.
#'
#' @examples
#'
#' COMPLETE ME
#'
#'
#'
#' @export
expo.brglmFit <- function(object, type = c("correction*", "correction+", "Lylesetal2012", "AS_median", "ML"), level = 0.95) {
    type <- match.arg(type)
    object_type <- object$type
    to_correct <- type %in% c("correction*", "correction+", "Lylesetal2012")
    if (to_correct) {
        if (object_type != "AS_mixed") {
            object <- update(object, type = "AS_mixed", data = object$data)
        }
    } else {
        ## AS_mean and ML are invariant to transformation
        if (object_type != type) {
            object <- update(object, type = type, data = object$data)
        }
    }
    info_fun <- get_information_function(object)
    dispersion <- object$dispersion
    coefs <- coef(object, model = "mean")
    ses <- sqrt(diag(solve(info_fun(coefs, dispersion))))
    trans_coefs <- exp(coefs) * switch(type,
                                       "ML" = 1,
                                       "AS_median" = 1,
                                       "Lylesetal2012" = exp(- ses^2/2),
                                       "correction*" = 1 / (1 + ses^2/2),
                                       "correction+" = (1 - ses^2/2))
    if (to_correct) {
        ## Recompute coefs and ses on the original scale from trans_coefs if method is not invariant
        coefs <- log(trans_coefs)
        ses <- sqrt(diag(solve(info_fun(coefs, dispersion))))
    }
    trans_ses <- trans_coefs * ses
    a <- (1 - level)/2
    ci <- coefs + qnorm(a) * cbind(ses, -ses)
    colnames(ci) <- paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE, digits = 3), "%")
    out <- list(coef = trans_coefs, se = trans_ses, ci = exp(ci), type = type)
    out$call <- match.call()
    class(out) <- c("brglmFit_expo", class(out))
    out
}

#' @method print brglmFit_expo
#' @export
print.brglmFit_expo <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cmat <- cbind(x$coef, x$se, x$ci)
    colnames(cmat) <- c("Estimate", "Std. Error", colnames(x$ci))
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    printCoefmat(cmat, digits = digits, ...)
    cat("\n\nType of estimator:", x$type, get_type_description(x$type), "\n")
}

#' Extract estimates from [brglmFit_expo] objects
#'
#' @inheritParams stats::coef
#'
#' @method coef brglmFit_expo
#' @export
coef.brglmFit_expo <- function(object, ...) {
    object$coef
}
