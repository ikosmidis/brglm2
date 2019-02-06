#' brglm2: Bias Reduction in Generalized Linear Models
#'
#' Estimation and inference from generalized linear models based on
#' implicit methods for bias reduction (see Kosmidis, 2014). brglm2
#' can achieve reduction of estimation bias by solving either the
#' mean-bias reducing adjusted score equations in Firth (1993) and
#' Kosmidis and Firth (2009), or the median-bias reducing adjusted
#' score equations in Kenne et al. (2016), or through the direct
#' subtraction of an estimate of the bias of the maximum likelihood
#' estimator from the maximum likelihood estimates as in Cordeiro and
#' McCullagh (1991),
#'
#'
#' In the special case of generalized linear models for binomial,
#' Poisson and multinomial responses, the adjusted score equations
#' approaches return estimates with improved frequentist properties,
#' that are also always finite, even in cases where the maximum
#' likelihood estimates are infinite (e.g. complete and quasi-complete
#' separation in multinomial regression; see also
#' \code{\link{detect_separation}} and
#' \code{\link{check_infinite_estimates}} for pre-fit and post-fit
#' methods for the detection of infinite estimates in binomial
#' response generalized linear models). Estimation in all cases takes
#' place via a modified Fisher scoring algorithm, and S3 methods for
#' the construction of confidence intervals for the reduced-bias
#' estimates are provided.
#'
#' @details
#'
#'
#' The similarly named **brglm** R package can only handle generalized
#' linear models with binomial responses. Special care has been taken
#' when developing **brglm2** in order not to have conflicts when the
#' user loads **brglm2** and **brglm** simultaneously. The development
#' and maintenance of the two packages will continue, until **brglm2**
#' incorporates all **brglm** functionality and gets an appropriate
#' wrapper to the \code{brglm::brglm} function.
#'
#' @author Ioannis Kosmidis \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @references
#'
#' Cordeiro G. M. & McCullagh, P. (1991). Bias correction in generalized
#' linear models. *Journal of the Royal Statistical Society. Series B
#' (Methodological)*, **53**, 629-643
#'
#' Firth D. (1993). Bias reduction of maximum likelihood estimates,
#' Biometrika, **80**, 27-38
#'
#' Kenne Pagui E C, Salvan A and Sartori N (2016). Median bias
#' reduction of maximum likelihood estimates. *arXiv*,
#' **arXiv:1604.04768**
#'
#' Kosmidis I and Firth D (2009). Bias reduction in exponential family
#' nonlinear models. *Biometrika*, **96**, 793-804
#'
#' Kosmidis I and Firth D (2010). A generic algorithm for reducing
#' bias in parametric estimation. *Electronic Journal of Statistics*,
#' **4**, 1097-1112
#'
#' Kosmidis I (2014). Bias in parametric estimation: reduction and
#' useful side-effects. *WIRE Computational Statistics*, **6**,
#' 185-196
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
#'
#'
NULL

## Suggestion by Kurt Hornik to avoid a warning related to the binding
## of n which is evaluated by family$initialize
if (getRversion() >= "2.15.1") globalVariables(c("n", "lambda"))

#' Generic method for checking for infinite estimates
#' @param object a fitted model object (e.g. the result of a
#'     \code{\link{glm}} call).
#' @param ... other options to be passed to the method.
#' @export
check_infinite_estimates <- function(object, ...) {
    UseMethod("check_infinite_estimates")
}
