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

#' Penalized likelihood ratio test for [`"mdyplFit"`][mdyplFit()] objects
#'
#' Computes the Diaconis-Ylvisaker prior penalized likelihood ratio
#' test statistic or its adjusted version using high-dimensionality
#' correction under proportional asymptotics. Associated p-values are
#' also computed using a chi squared distribution.
#'
#' @aliases plrtest
#' @inheritParams summary.mdyplFit
#' @param object1 a [`"mdyplFit"`][mdyplFit()] object
#' @param object2 a [`"mdyplFit"`][mdyplFit()] object
#' @param ... further arguments to be passed to methods. Currently not used.
#'
#' @details
#'
#' Both `object1` and `object2` should have been fitted using the
#' [mdyplFit()] method for [glm()], and the same shrinkage parameter
#' `alpha`; see [mdyplFit()] and [mdyplControl()] for setting `alpha`.
#'
#' If `hd_correction = TRUE` then the deviance and the associated
#' p-value are adjusted using a high-dimensionality correction under
#' proportional asymptotics as in Sterzinger & Kosmidis (2024); see
#' [summary.mdyplFit()].
#'
#' @author Ioannis Kosmidis `[aut, cre]` \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso [mdyplFit()], [summary.mdyplFit()], [mdypl_control()]
#'
#' @references
#'
#' Sterzinger P, Kosmidis I (2024). Diaconis-Ylvisaker prior
#' penalized likelihood for \eqn{p/n \to \kappa \in (0,1)} logistic
#' regression. *arXiv*:2311.07419v2, \url{https://arxiv.org/abs/2311.07419}.
#'
#' @method plrtest mdyplFit
#' @export
plrtest.mdyplFit <- function(object1, object2, hd_correction = FALSE, ...) {
    if (!inherits(object1, "mdyplFit") || !inherits(object2, "mdyplFit")) {
        stop("At least one of `object1` and `object2` is not of class `mdyplFit`")
    }
    if (!identical(object1$alpha, object2$alpha)) {
        stop("`object1` and `object2` have not been fitted using the same `alpha`. See ?mdypl_control.")
    }
    anv <- anova(object1, object2, test = "LRT")
    if (isTRUE(hd_correction)) {
        npars1 <- object1$rank
        npars2 <- object2$rank
        if (npars2 > npars1) {
            summ <- summary(object2, hd_correction = TRUE)
        } else {
            summ <- summary(object1, hd_correction = TRUE)
        }
        se_pars <- summ$se_parameters
        anv$Deviance <- anv$Deviance * se_pars[2] / (summ$kappa * se_pars[3]^2)
        anv[, "Pr(>Chi)"] <- pchisq(anv$Deviance, abs(anv$Df), lower.tail = FALSE)
        attr(anv, "hd_correction") <- TRUE
        attr(anv, "kappa") <- summ$kappa
        attr(anv, "signal_strength") <- summ$signal_strength
        attr(anv, "se_parameters") <- se_pars
    }
    class(anv) <- c("anova_mdyplFit", class(anv))
    anv
}

#' @method print anova_mdyplFit
#' @export
print.anova_mdyplFit <- function(x, ...) {
    print_anova(x, ...)
    if (isTRUE(attr(x, "hd_correction"))) {
        se_pars <- attr(x, "se_parameters")
        cat("\nHigh-dimensionality correction applied with")
        cat("\nDimentionality parameter (kappa) =", round(attr(x, "kappa"), 2))
        cat("\nEstimated signal strength (gamma) =", round(attr(x, "signal_strength"), 2))
        cat("\nState evolution parameters (mu, b, sigma) =", paste0("(", paste(round(se_pars[1:3], 2), collapse = ", "), ")"), "with max(|funcs|) =", max(abs(attr(se_pars, "funcs"))), "\n")
    }
}

print_anova <- utils::getFromNamespace("print.anova", "stats")
