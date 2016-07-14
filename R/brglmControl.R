#' Auxiliary function for \code{\link{glm}} fitting using the
#' \code{\link{brglmFit}} method.  Typically only used internally by
#' \code{\link{brglmFit}}, but may be used to construct a
#' \code{control} argument to either function.
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
#' @return a list with components named as the arguments, including
#'     symbolic expressions for the dispersion transformation
#'     (\code{Trans}) and its inverse (\code{inverseTrans})
#'
#' @seealso \code{\link{brglmFit}} and \code{\link{glm.fit}}
#'
brglmControl <- function(epsilon = 1e-10, maxit = 100,
                         trace = FALSE,
                         type = c("adjusted_scores", "correction", "maximum_likelihood"),
                         dispTrans = "identity",
                         slowit = 1,
                         maxStepFactor = 1) {
    type <- match.arg(type)
    Trans <- switch(dispTrans,
                    identity = expression(disp),
                    sqrt = expression(disp^0.5),
                    inverse = expression(1/disp),
                    log = expression(log(disp)),
                    inverseSqrt = expression(1/sqrt(disp)),
                    custom = customTrans[[1]])
    inverseTrans <- switch(dispTrans,
                           identity = expression(transdisp),
                           sqrt = expression(transdisp^2),
                           inverse = expression(1/transdisp),
                           log = expression(exp(transdisp)),
                           inverseSqrt = expression(1/transdisp^2),
                           custom = customTrans[[2]])
    if (!(dispTrans %in% c("identity",
                           "sqrt",
                           "inverse",
                           "log",
                           "custom",
                           "inverseSqrt"))) {
    stop(dispTrans, " is not a supported transofrmation of the dispersion")
    }
    ## if (dispTrans == "custom") {
    ##     warning("Please note that bo check has been made whether Trans and inverseTrans agree")
    ## }
    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("value of 'epsilon' must be > 0")
    ## if (!is.numeric(maxit) || maxit <= 0)
    ##     stop("maximum number of iterations must be > 0")
    list(epsilon = epsilon, maxit = maxit, trace = trace,
       type = type,
       Trans = Trans,
       inverseTrans = inverseTrans,
       dispTrans = dispTrans,
       slowit = slowit,
       maxStepFactor = maxStepFactor)
}


