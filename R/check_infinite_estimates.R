# Copyright (C) 2016-2020 Ioannis Kosmidis

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


#' A simple diagnostic of whether the maximum likelihood estimates are
#' infinite
#'
#' @aliases checkInfiniteEstimates
#'
#' @param object the result of a \code{\link{glm}} call.
#' @param nsteps starting from \code{maxit = 1}, the GLM is refitted
#'     for \code{maxit = 2}, \code{maxit = 3}, \ldots, \code{maxit =
#'     nsteps}. Default value is 30.
#' @param ... currently not used.
#'
#'
#' @details
#'
#' \code{check_infinite_estimates} attempts to identify the occurrence
#' of infinite estimates in GLMs with binomial responses by
#' successively refitting the model. At each iteration the maximum
#' number of allowed IWLS iterations is fixed starting from 1 to
#' \code{nsteps} (by setting \code{control = glm.control(maxit = j)},
#' where \code{j} takes values 1, \ldots, nsteps in
#' \code{\link{glm}}). For each value of \code{maxit}, the estimated
#' asymptotic standard errors are divided to the corresponding ones
#' from \code{control = glm.control(maxit = 1)}. Then, based on the
#' results in Lesaffre & Albert (1989), if the sequence of ratios in
#' any column of the resultant matrix diverges, then complete or
#' quasi-complete separation occurs and the maximum likelihood
#' estimate for the corresponding parameter has value minus or plus
#' infinity.
#'
#' @note
#'
#' \code{check_infinite_estimates} will be removed from \pkg{brglm2}
#' at version 0.8. An new version of \code{check_infinite_estimates}
#' is now maintained in the \pkg{detectseparation} R package at
#' \url{https://cran.r-project.org/package=detectseparation}. In order
#' to use the version in \code{detect_separation} load first
#' \pkg{brglm2} and then \pkg{detectseparation}, i.e.
#' \code{library(brglm2); library(detectseparation)}.
#'
#' @seealso \code{\link[nnet]{multinom}}, \code{\link{brmultinom}}
#'
#' @references
#'
#' Kosmidis I, Firth D (2020). Jeffreys-prior penalty, finiteness
#' and shrinkage in binomial-response generalized linear
#' models. *Biometrika* \doi{10.1093/biomet/asaa052}
#'
#' Lesaffre E, Albert A (1989). Partial Separation in Logistic
#' Discrimination. *Journal of the Royal Statistical Society. Series B
#' (Methodological)*, **51**, 109-116 \doi{10.1111/j.2517-6161.1989.tb01752.x}
#'
#' @examples
#'
#' ## endometrial data from Heinze \& Schemper (2002) (see ?endometrial)
#' data("endometrial", package = "brglm2")
#' endometrialML <- glm(HG ~ NV + PI + EH, data = endometrial,
#'                      family = binomial("probit"))
#' ## clearly the maximum likelihood estimate for the coefficient of
#' ## NV is infinite
#' check_infinite_estimates(endometrialML)
#'
#' \dontrun{
#' ## Aligator data (Agresti, 2002, Table~7.1)
#' data("alligator", package = "brglm2")
#' all_ml <- brmultinom(foodchoice ~ size + lake , weights = round(freq/3),
#'                      data = alligators, type = "ML", ref = 1)
#' ## Clearly some estimated standard errors diverge as the number of
#' ## Fisher scoring iterations increases
#' matplot(check_infinite_estimates(all_ml), type = "l", lty = 1,
#'         ylim = c(0.5, 1.5))
#' }
#' @export
check_infinite_estimates.glm <- function(object, nsteps = 20, ...)
{
    ## Deprecation warning comes from method in `brglm2-package.R`
    is_brmultinom <- inherits(object, "brmultinom")

    if ((class(object)[1] != "glm") & (!is_brmultinom)) {
        warning("check_infinite_estimates has been designed for objects of class 'glm'")
    }
    if ((object$family$family != "binomial") & (!is_brmultinom)) {
        warning("check_infinite_estimates has been designed for binomial-response models")
    }

    if (is_brmultinom) {
        betas <- coef(object)
        dims <- dim(betas)
        betasNames <- paste(rep(colnames(betas), dims[1]), rownames(betas), sep = ":")
        betas <- c(betas)
        names(betas) <- betasNames
    }
    else {
        betas <- coef(object)
        betasNames <- names(betas)
    }
    eps <- .Machine$double.eps
    noNA <- !is.na(betas)
    stdErrors <- matrix(0, nsteps, length(betas))
    start <- NULL
    for (i in 1:nsteps) {
        if (is_brmultinom) {
            suppressWarnings(temp.object <- update(object, control = list(maxit = i, epsilon = eps, type = "ML"), start = start))
            stdErrors[i, noNA] <- sqrt(diag(vcov(temp.object))[noNA])
        }
        else {
            suppressWarnings(temp.object <- update(object, control = list(maxit = i, epsilon = eps), start = start))

            stdErrors[i, noNA] <- summary(temp.object)$coef[betasNames[noNA], "Std. Error"]
        }
        start <- c(coef(temp.object))
    }
    res <- sweep(stdErrors, 2, stdErrors[1, ], "/")
    colnames(res) <- betasNames
    res
}


