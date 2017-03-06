#' Auxiliary function for \code{\link{glm}} fitting using the
#' \code{\link{brglmFit}} method.
#'
#' Typically only used internally by \code{\link{brglmFit}}, but may
#' be used to construct a \code{control} argument.
#'
#' @inheritParams stats::glm.control
#' @param epsilon positive convergence tolerance epsilon
#' @param maxit integer giving the maximal number of iterations
#'     allowed
#' @param trace logical indicating if output should be produced for
#'     each iteration
#' @param type type of estimation method. Supported methods include
#'     \code{adjusted_scores}, \code{correction} and
#'     \code{maximum_likelihood}
#' @param dispTrans the transformation of the dispersion to be
#'     estimated. Default is \code{identity}. See Details.
#' @param slowit a positive real used as a multiplier for the
#'     stepsize. The smaller it is the smaller the steps are
#' @param maxStepFactor the maximum number of step halving steps to
#'     consider
#'
#' @details \code{\link{brglmControl}} provides default values and
#'     sanity checking for the various constants that control the
#'     iteration and generally the behaviour of
#'     \code{\link{brglmFit}}.
#'
#'      When \code{trace} is true, calls to \code{cat} produce the
#'      output for each iteration.  Hence, \code{options(digits = *)}
#'      can be used to increase the precision.
#'
#'      dispTrans sets the transformation of the dispersion parameter
#'      for which the bias reduced estimates are computed. Can be one
#'      of "identity", "sqrt", "inverse", "log" and
#'      "inverseSqrt". Custom transformations are accommodated by
#'      supplying a list of two expressions (transformation and
#'      inverse transformation). See the examples for more details.
#'
#'
#' @return a list with components named as the arguments, including
#'     symbolic expressions for the dispersion transformation
#'     (\code{Trans}) and its inverse (\code{inverseTrans})
#'
#' @seealso \code{\link{brglmFit}} and \code{\link{glm.fit}}
#'
#' @examples
#'
#' data("coalition", package = "Zelig")
#' ## The maximum likelihood fit with log link
#' coalitionML <- glm(duration ~ fract + numst2, family = Gamma, data = coalition)
#'
#' ## Bias-reduced estimation of the dispersion parameter
#' coalitionBRi <- update(coalitionML, method = "brglmFit")
#' coef(coalitionBRi, model = "dispersion")
#'
#' ## Bias-reduced estimation of log(dispersion)
#' coalitionBRl <- update(coalitionML, method = "brglmFit",  dispTrans = "log")
#' coef(coalitionBRl, model = "dispersion")
#'
#' ## Just for illustration: Bias-reduced estimation of dispersion^0.25
#' coalitionBRc <- update(coalitionML, method = "brglmFit",
#'                        dispTrans = list(expression(disp^0.25), expression(transdisp^4)))
#' coef(coalitionBRc, model = "dispersion")
#'
brglmControl <- function(epsilon = 1e-10, maxit = 100,
                         trace = FALSE,
                         type = c("adjusted_scores", "correction", "maximum_likelihood"),
                         dispTrans = "identity",
                         slowit = 1,
                         maxStepFactor = 1) {
    type <- match.arg(type)

    if (is.character(dispTrans)) {
        Trans <- switch(dispTrans,
                        identity = expression(disp),
                        sqrt = expression(disp^0.5),
                        inverse = expression(1/disp),
                        log = expression(log(disp)),
                        inverseSqrt = expression(1/sqrt(disp)),
                        stop(dispTrans, " is not one of the implemented dispersion transformations"))
        inverseTrans <- switch(dispTrans,
                               identity = expression(transdisp),
                               sqrt = expression(transdisp^2),
                               inverse = expression(1/transdisp),
                               log = expression(exp(transdisp)),
                               inverseSqrt = expression(1/transdisp^2))
    }
    else {
        if (is.list(dispTrans) && (length(dispTrans) == 2)) {
            Trans <- dispTrans[[1]]
            inverseTrans <- dispTrans[[2]]
            dispTrans <- "custom_transformation"
        }
        else {
            stop("dispTrans can be either one of 'identity', 'sqrt', 'inverse', 'log' and 'inverseSqrt', or a list of two expressions")
        }
    }
    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("value of 'epsilon' must be > 0")
    list(epsilon = epsilon, maxit = maxit, trace = trace,
       type = type,
       Trans = Trans,
       inverseTrans = inverseTrans,
       dispTrans = dispTrans,
       slowit = slowit,
       maxStepFactor = maxStepFactor)
}


