# Copyright (C) 2020- Ioannis Kosmidis

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
#' multinomial regression). Estimation in all cases takes place via a
#' modified Fisher scoring algorithm, and S3 methods for the
#' construction of confidence intervals for the reduced-bias estimates
#' are provided.
#'
#' The core model fitters are implemented by the functions
#' [brglm_fit()] (univariate generalized linear models),
#' [brmultinom()] (baseline category logit models for nominal
#' multinomial responses), [bracl()] (adjacent category logit models
#' for ordinal multinomial responses), and [brnb()] for negative
#' binomial regression.
#'
#' @details
#'
#'
#' The similarly named **brglm** R package can only handle generalized
#' linear models with binomial responses. Special care has been taken
#' when developing **brglm2** in order not to have conflicts when the
#' user loads **brglm2** and **brglm** simultaneously. The development
#' and maintenance of the two packages will continue in parallel,
#' until **brglm2** incorporates all **brglm** functionality and
#' provides an appropriate wrapper to the [brglm::brglm()] function.
#'
#' @author Ioannis Kosmidis `[aut, cre]` \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso
#'
#' [brglm_fit()], [brmultinom()], [bracl()]
#'
#' @references
#'
#' Kosmidis I, Firth D (2021). Jeffreys-prior penalty, finiteness
#' and shrinkage in binomial-response generalized linear
#' models. *Biometrika*, **108**, 71-82. \doi{10.1093/biomet/asaa052}.
#'
#' Cordeiro G M, McCullagh P (1991). Bias correction in generalized
#' linear models. *Journal of the Royal Statistical Society. Series B
#' (Methodological)*, **53**, 629-643. \doi{10.1111/j.2517-6161.1991.tb01852.x}.
#'
#' Firth D (1993). Bias reduction of maximum likelihood estimates,
#' Biometrika, **80**, 27-38. \doi{10.2307/2336755}.
#'
#' Kenne Pagui E C, Salvan A, Sartori N (2017). Median bias
#' reduction of maximum likelihood estimates. *Biometrika*, **104**,
#' 923–938. \doi{10.1093/biomet/asx046}.
#'
#' Kosmidis I, Kenne Pagui E C, Sartori N (2020). Mean and median bias
#' reduction in generalized linear models. *Statistics and Computing*,
#' **30**, 43-59. \doi{10.1007/s11222-019-09860-6}.
#'
#' Kosmidis I, Firth D (2009). Bias reduction in exponential family
#' nonlinear models. *Biometrika*, **96**, 793-804. \doi{10.1093/biomet/asp055}.
#'
#' Kosmidis I, Firth D (2010). A generic algorithm for reducing
#' bias in parametric estimation. *Electronic Journal of Statistics*,
#' **4**, 1097-1112. \doi{10.1214/10-EJS579}.
#'
#' Kosmidis I (2014). Bias in parametric estimation: reduction and
#' useful side-effects. *WIRE Computational Statistics*, **6**,
#' 185-196. \doi{10.1002/wics.1296}.
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

#' @export
ordinal_superiority <- function(object, formula, data,
                                measure = c("gamma", "Delta"),
                                level = 0.95,
                                bc = FALSE) {
    UseMethod("ordinal_superiority")
}

#' @export
expo <- function(object, type = c("ML", "correction", "AS_median", "Lylesetal2012"), level = 0.95) {
    UseMethod("expo")
}
