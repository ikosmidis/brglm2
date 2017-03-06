#' brglm2: Estimation and inference from generalized linear models
#'    using explicit and implicit methods for bias reduction
#'
#' Estimation and inference from generalized linear models based on
#' implicit methods for bias reduction (see Kosmidis, 2014, WIRE
#' computational statistics). brglm can achieve reduction of
#' estimation bias either through the adjusted score equations
#' approach in Firth (1993, Biometrika) and Kosmidis and Firth (2009,
#' Biometrika), or through the direct subtraction of an estimate of
#' the bias of the maximum likelihood estimator from the maximum
#' likelihood estimates. In the special case of generalized linear
#' models for binomial and multinomial responses, the adjusted score
#' equations approach returns estimates with improved frequentist
#' properties, that are also always finite, even in cases where the
#' maximum likelihood estimates are infinite (e.g. complete and
#' quasi-complete separation). Estimation in all cases takes place via
#' a modified Fisher scoring algorithm, and S3 methods for the
#' construction of confidence intervals for the reduced-bias estimates
#' are provided.
#'
#' @docType package
#' @name brglm2
#' @import stats
#' @import MASS
#' @import enrichwith
#' @importFrom graphics plot
#' @importFrom nnet class.ind
#' @import Matrix
#'
NULL
#> NULL

## Suggestion by Kurt Hornik to avoid a warning related to the binding
## of n which is evaluated by family$initialize
if (getRversion() >= "2.15.1") globalVariables(c("n", "customTrans", "lambda"))
