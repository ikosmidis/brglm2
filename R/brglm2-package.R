# Copyright (C) 2020-2021 Ioannis Kosmidis

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

#'
#' brglm2: Bias Reduction in Generalized Linear Models
#'
#' Estimation and inference from generalized linear models using
#' implicit and explicit bias reduction methods (Kosmidis, 2014), and
#' other penalized maximum likelihood methods. Currently supported
#' methods include the mean bias-reducing adjusted scores approach in
#' Firth (1993) and Kosmidis & Firth (2009), the median bias-reduction
#' adjusted scores approach in Kenne Pagui et al. (2017), the
#' correction of the asymptotic bias in Cordeiro & McCullagh (1991),
#' the mixed bias-reduction adjusted scores approach in Kosmidis et al
#' (2020), maximum penalized likelihood with powers of the Jeffreys
#' prior as penalty, and maximum likelihood.
#'
#'
#' In the special case of generalized linear models for binomial,
#' Poisson and multinomial responses (both nominal and ordinal), mean
#' and median bias reduction and maximum penalized likelihood return
#' estimates with improved frequentist properties, that are also
#' always finite, even in cases where the maximum likelihood estimates
#' are infinite (e.g. complete and quasi-complete separation in
#' multinomial regression; see also \code{\link{detect_separation}}
#' and \code{\link{check_infinite_estimates}} for pre-fit and post-fit
#' methods for the detection of infinite estimates in binomial
#' response generalized linear models). Estimation in all cases takes
#' place via a modified Fisher scoring algorithm, and S3 methods for
#' the construction of confidence intervals for the reduced-bias
#' estimates are provided.
#'
#' The core model fitters are implemented by the functions
#' \code{\link{brglm_fit}} (univariate generalized linear models),
#' \code{\link{brmultinom}} (baseline category logit models for
#' nominal multinomial responses), and \code{\link{bracl}} (adjacent
#' category logit models for ordinal multinomial responses).
#'
#' @details
#'
#'
#' The similarly named **brglm** R package can only handle generalized
#' linear models with binomial responses. Special care has been taken
#' when developing **brglm2** in order not to have conflicts when the
#' user loads **brglm2** and **brglm** simultaneously. The development
#' and maintenance of the two packages will continue in parallel,
#' until **brglm2** incorporates all **brglm** functionality and gets
#' an appropriate wrapper to the \code{brglm::brglm} function.
#'
#' @author Ioannis Kosmidis \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso
#'
#' \code{\link{brglm_fit}}, \code{\link{brmultinom}}, \code{\link{bracl}}
#'
#' @references
#'
#' Kosmidis I, Firth D (2021). Jeffreys-prior penalty, finiteness
#' and shrinkage in binomial-response generalized linear
#' models. *Biometrika*, **108**, 71-82 \doi{10.1093/biomet/asaa052}
#'
#' Cordeiro G M, McCullagh P (1991). Bias correction in generalized
#' linear models. *Journal of the Royal Statistical Society. Series B
#' (Methodological)*, **53**, 629-643 \doi{10.1111/j.2517-6161.1991.tb01852.x}
#'
#' Firth D (1993). Bias reduction of maximum likelihood estimates,
#' Biometrika, **80**, 27-38 \doi{10.2307/2336755}
#'
#' Kenne Pagui E C, Salvan A, Sartori N (2017). Median bias
#' reduction of maximum likelihood estimates. *Biometrika*, **104**,
#' 923–938 \doi{10.1093/biomet/asx046}
#'
#' Kosmidis I, Kenne Pagui E C, Sartori N (2020). Mean and median bias
#' reduction in generalized linear models. *Statistics and Computing*,
#' **30**, 43-59 \doi{10.1007/s11222-019-09860-6}
#'
#' Kosmidis I, Firth D (2009). Bias reduction in exponential family
#' nonlinear models. *Biometrika*, **96**, 793-804 \doi{10.1093/biomet/asp055}
#'
#' Kosmidis I, Firth D (2010). A generic algorithm for reducing
#' bias in parametric estimation. *Electronic Journal of Statistics*,
#' **4**, 1097-1112 \doi{10.1214/10-EJS579}
#'
#' Kosmidis I (2014). Bias in parametric estimation: reduction and
#' useful side-effects. *WIRE Computational Statistics*, **6**,
#' 185-196 \doi{10.1002/wics.1296}
#'
#' @docType package
#' @name brglm2
#' @import stats
#' @import enrichwith
#' @import Matrix
#' @import MASS
#' @importFrom graphics plot
#' @importFrom nnet class.ind
#' @importFrom numDeriv grad
#' @useDynLib brglm2
#'
NULL

## NAMESPACE should have import(stats), import(Matrix)


## Suggestion by Kurt Hornik to avoid a warning related to the binding
## of n which is evaluated by family$initialize
if (getRversion() >= "2.15.1") globalVariables(c("n", "lambda"))

#' Ordinal superiority scores of Agresti and Kateri (2017)
#'
#' \code{\link{ordinal_superiority}} is a method for the estimation
#' and inference about model-based ordinal superiority scores
#' introduced in Agresti and Kateri (2017, Section 5) from fitted
#' objects. The mean bias of the estimates of the ordinal superiority
#' scores can be corrected.
#'
#' @param object a fitted object from an ordinal regression
#'     model. Currently only models from class \code{"bracl"} are
#'     supported.
#' @param formula a RHS formula indicating the group variable to use.
#' @param data an optional data frame in which to look for variables
#'     with which to compute ordinal superiority measures.  If
#'     omitted, an attempt is made to use the data that produced
#'     \code{object}.
#' @param measure either \code{"gamma"} (default) or \code{"Delta"},
#'     specifying the ordinal superiority measure to be returned.
#' @param level the confidence level required when computing
#'     confidence intervals for the ordinal superiority measures.
#' @param bc logical. If \code{FALSE} (default) then the ordinal
#'     superiority measures are computed using the estimates in
#'     \code{object}. If \code{TRUE} then the ordinal superiority
#'     measure estimates are corrected for mean bias.
#'
#' @examples
#' data("stemcell", package = "brglm2")
#'
#' # Adjacent category logit (proportional odds)
#' stem <- within(stemcell, {nreligion = as.numeric(religion)})
#' fit_bracl_p <- bracl(research ~ nreligion + gender, weights = frequency,
#'                      data = stem, type = "ML", parallel = TRUE)
#'
#' # Estimates and 95% confidence intervals for the probabilities that the response
#' # category for gender "female" is higher than the response category for gender "male",
#' # while adjusting for religion.
#' ordinal_superiority(fit_bracl_p, ~ gender)
#'
#' \dontrun{
#' # And their (very-similar in value here) bias corrected versions
#' # with 99% CIs
#' ordinal_superiority(fit_bracl_p, ~ gender, bc = TRUE, level = 0.99)
#' # Note that the object is refitted with type = "AS_mean"
#'
#' }
#'
#'
#' @references
#'
#' Agresti, A., Kateri, M. (2017). Ordinal probability effect measures
#' for group comparisons in multinomial cumulative link models.
#' *Biometrics*, **73** 214-219 #' \doi{10.1111/biom.12565}
#' @export
ordinal_superiority <- function(object, formula, data,
                                measure = c("gamma", "Delta"),
                                level = 0.95,
                                bc = FALSE) {
    UseMethod("ordinal_superiority")
}
