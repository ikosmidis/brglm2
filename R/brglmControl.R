# Copyright (C) 2016- Ioannis Kosmidis

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


#' Auxiliary function for [glm()] fitting using the [brglmFit()]
#' method.
#'
#' Typically only used internally by [brglmFit()], but may be used to
#' construct a `control` argument.
#'
#' @aliases brglm_control
#' @param epsilon positive convergence tolerance epsilon. Default is
#'     `1e-06`.
#' @param check_aliasing logical indicating where a QR decomposition
#'     of the model matrix should be used to check for
#'     aliasing. Default is `TRUE`. See Details.
#' @param maxit integer giving the maximal number of iterations
#'     allowed. Default is `100`.
#' @param trace logical indicating if output should be produced for
#'     each iteration. Default is `FALSE`.
#' @param type the type of fitting method to be used. The options are
#'     `"AS_mean"` (mean-bias reducing adjusted scores), `"AS_median"`
#'     (median-bias reducing adjusted scores), `"AS_mixed"` (bias
#'     reduction using mixed score adjustments; default),
#'     `"correction"` (asymptotic bias correction), `"MPL_Jeffreys"`
#'     (maximum penalized likelihood with powers of the Jeffreys prior
#'     as penalty) and `"ML"` (maximum likelihood).
#' @param transformation the transformation of the dispersion to be
#'     estimated. Default is `"identity"`. See Details.
#' @param slowit a positive real used as a multiplier for the
#'     stepsize. The smaller it is the smaller the steps are. Default
#'     is `1`.
#' @param max_step_factor the maximum number of step halving steps to
#'     consider. Default is `12`.
#' @param response_adjustment a (small) positive constant or a vector
#'     of such. Default is `NULL`. See Details.
#' @param a power of the Jeffreys prior penalty. See Details.
#' @param ... further arguments passed to [brglmControl()]. Currently
#'     ignored in the output.
#'
#' @details
#'
#' [brglmControl()] provides default values and sanity checking for
#' the various constants that control the iteration and generally the
#' behaviour of [brglmFit()].
#'
#' When `trace = TRUE`, calls to [cat()] produce the output for each
#' iteration.  Hence, `options(digits = *)` can be used to increase
#' the precision.
#'
#' When `check_aliasing = TRUE` (default), a QR decomposition of the
#' model matrix is computed to check for aliasing. If the model matrix
#' is known to be of full rank, then `check_aliasing = FALSE` avoids
#' the extra computational overhead of an additional QR decomposition,
#' which can be substantial for large model matrices. However, setting
#' `check_aliasing = FALSE` tells [brglmFit()] that the model matrix
#' is full rank, and hard to trace back errors will result if it is
#' rank deficient.
#'
#' `transformation` sets the transformation of the dispersion
#' parameter for which the bias reduced estimates are computed. Can be
#' one of `"identity"`, `"sqrt"`, `"inverse"`, `"log"` and
#' `"inverseSqrt"`. Custom transformations are accommodated by
#' supplying a list of two expressions (transformation and inverse
#' transformation). See the examples for more details.
#'
#' The value of `response_adjustment` is only relevant if [brglmFit()]
#' is called with `start = NULL`, and `family` is [binomial()] or
#' [poisson()]. For those models, an initial maximum likelihood fit is
#' obtained on adjusted data to provide starting values for the
#' iteration in [brglmFit()]. The value of `response_adjustment`
#' governs how the data is adjusted. Specifically, if `family` is
#' [binomial()], then the responses and totals are adjusted by
#' `response_adjustment` and `2 * response_adjustment`, respectively;
#' if `family` is [poisson()], then the responses are adjusted by and
#' `response_adjustment`. `response_adjustment = NULL` (default)
#' is equivalent to setting it to
#' `"number of parameters" / "number of observations"`.
#'
#' When `type = "AS_mixed"` (default), mean bias reduction is used for
#' the regression parameters, and median bias reduction for the
#' dispersion parameter, if that is not fixed. This adjustment has
#' been developed based on equivariance arguments (see, Kosmidis et
#' al, 2020, Section 4) in order to produce regression parameter
#' estimates that are invariant to arbitrary contrasts, and estimates
#' for the dispersion parameter that are invariant to arbitrary
#' non-linear transformations. `type = "AS_mixed"` and `type =
#' "AS_mean"` return the same results if [brglmFit()] is called with
#' `family` [binomial()] or [poisson()] (i.e. families with fixed
#' dispersion).
#'
#' When `type = "MPL_Jeffreys"`, [brglmFit()] will maximize the
#' penalized log-likelihood \deqn{l(\beta, \phi) + a\log \det i(\beta,
#' \phi)}{l(beta, phi) + a log det i(beta, phi)} where \eqn{i(\beta,
#' \phi)}{i(beta, phi)} is the expected information matrix about the
#' regression parameters \eqn{\beta} and the dispersion parameter
#' \eqn{\phi}. See, `vignette("iteration", "brglm2")` for more
#' information. The argument `a` controls the amount of penalization
#' and its default value is `a = 1/2`, corresponding to maximum
#' penalized likelihood using a Jeffreys-prior penalty. See, Kosmidis
#' & Firth (2021) for proofs and discussion about the finiteness and
#' shrinkage properties of the maximum penalized likelihood estimators
#' for binomial-response generalized linear models.
#'
#' The estimates from `type = "AS_mean"` and `type =
#' "MPL_Jeffreys"` with `a = 1/2` (default) are identical for
#' Poisson log-linear models and logistic regression models, i.e. for
#' binomial and Poisson regression models with canonical links. See,
#' Firth (1993) for details.
#'
#' [brglm_control()] is an alias to [brglmControl()].
#'
#' @return a list with components named as the arguments, including
#'     symbolic expressions for the dispersion transformation
#'     (`Trans`) and its inverse (`inverseTrans`)
#'
#' @author Ioannis Kosmidis `[aut, cre]` \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso [brglm_fit()] and [glm.fit()]
#'
#' @references
#'
#' Kosmidis I, Firth D (2021). Jeffreys-prior penalty, finiteness
#' and shrinkage in binomial-response generalized linear
#' models. *Biometrika*, **108**, 71-82. \doi{10.1093/biomet/asaa052}.
#'
#' Kosmidis I, Kenne Pagui E C, Sartori N (2020). Mean and median bias
#' reduction in generalized linear models. *Statistics and Computing*,
#' **30**, 43-59. \doi{10.1007/s11222-019-09860-6}.
#'
#' Firth D (1993). Bias reduction of maximum likelihood estimates.
#' Biometrika, **80**, 27-38. \doi{10.2307/2336755}.
#'
#' @examples
#'
#' data("coalition", package = "brglm2")
#' ## The maximum likelihood fit with log link
#' coalitionML <- glm(duration ~ fract + numst2, family = Gamma, data = coalition)
#'
#' ## Bias reduced estimation of the dispersion parameter
#' coalitionBRi <- glm(duration ~ fract + numst2, family = Gamma, data = coalition,
#'                     method = "brglmFit")
#' coef(coalitionBRi, model = "dispersion")
#'
#' ## Bias reduced estimation of log(dispersion)
#' coalitionBRl <- glm(duration ~ fract + numst2, family = Gamma, data = coalition,
#'                     method = "brglmFit", transformation = "log")
#' coef(coalitionBRl, model = "dispersion")
#'
#' ## Just for illustration: Bias reduced estimation of dispersion^0.25
#' my_transformation <- list(expression(dispersion^0.25), expression(transformed_dispersion^4))
#' coalitionBRc <- update(coalitionBRi, transformation = my_transformation)
#' coef(coalitionBRc, model = "dispersion")
#'
#' @export
brglmControl <- function(epsilon = 1e-06, maxit = 100,
                         check_aliasing = TRUE,
                         trace = FALSE,
                         type = c("AS_mixed", "AS_mean", "AS_median", "correction", "MPL_Jeffreys", "ML"),
                         transformation = "identity",
                         slowit = 1,
                         response_adjustment = NULL,
                         max_step_factor = 12,
                         a = 1/2, ...) {
    type <- match.arg(type)

    if (is.character(transformation)) {
        Trans <- switch(transformation,
                        identity = expression(dispersion),
                        sqrt = expression(dispersion^0.5),
                        inverse = expression(1/dispersion),
                        log = expression(log(dispersion)),
                        inverseSqrt = expression(1/sqrt(dispersion)),
                        stop(transformation, " is not one of the implemented dispersion transformations"))
        inverseTrans <- switch(transformation,
                               identity = expression(transformed_dispersion),
                               sqrt = expression(transformed_dispersion^2),
                               inverse = expression(1/transformed_dispersion),
                               log = expression(exp(transformed_dispersion)),
                               inverseSqrt = expression(1/transformed_dispersion^2))
    } else {
        if (is.list(transformation) && (length(transformation) == 2)) {
            Trans <- transformation[[1]]
            inverseTrans <- transformation[[2]]
            transformation <- "custom_transformation"
        } else {
            stop("transformation can be either one of 'identity', 'sqrt', 'inverse', 'log' and 'inverseSqrt', or a list of two expressions")
        }
    }
    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("value of 'epsilon' must be > 0")

    if (!is.numeric(max_step_factor) || max_step_factor < 1) {
        warning("`max_step_factor = ", deparse(max_step_factor), "` is not a permissible value. Defaulting to 12")
        max_step_factor <- 12
    }
    list(epsilon = epsilon, maxit = maxit, trace = trace,
         check_aliasing = check_aliasing,
         response_adjustment = response_adjustment,
         type = type,
         Trans = Trans,
         inverseTrans = inverseTrans,
         transformation = transformation,
         slowit = slowit,
         max_step_factor = max_step_factor,
         a = a)
}

