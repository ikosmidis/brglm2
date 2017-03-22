#' A simple check of whether the maximum likelihood estimates are infinite
#'
#'
#' @param object the result of a \code{\link{glm}} call
#' @param nsteps starting from \code{maxit = 1}, the GLM is refitted
#'     for \code{maxit = 2}, \code{maxit = 3}, \ldots, \code{maxit =
#'     nsteps}. Default value is 30.
#' @param ... currently not used#'
#'
#'
#' @details
#'
#' \code{detect_infinite_estimates} attempts to identify the occurence
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
#' @seealso \code{\link[nnet]{multinom}}
#'
#' @references
#'
#' Lesaffre, E., & Albert, A. (1989). Partial Separation in Logistic
#' Discrimination. *Journal of the Royal Statistical Society. Series B
#' (Methodological)*, **51**, 109-116
#'
#' @examples
#'
#' ## endometrial data from Heinze \& Schemper (2002) (see ?endometrial)
#' data("endometrial", package = "brglm2")
#' endometrialML <- glm(HG ~ NV + PI + EH, data = endometrial,
#'                      family = binomial("probit"))
#' ## clearly the maximum likelihood estimate for the coefficient of
#' ## NV is infinite
#' detect_infinite_estimates(endometrialML)
#'
#' @export
detect_infinite_estimates.glm <- function (object, nsteps = 30, ...)
{
    if (class(object)[1] != "glm") {
        warning("detect_infinite_estimates has been designed for objects of class 'glm'")
    }
    if (object$family$family != "binomial") {
        warning("detect_infinite_estimates has been designed for binomial-response models")
    }
    eps <- .Machine$double.eps
    betasNames <- names(betas <- coef(object))
    noNA <- !is.na(betas)
    stdErrors <- matrix(0, nsteps, length(betas))
    for (i in 1:nsteps) {
        suppressWarnings(temp.object <- update(object, control = glm.control(maxit = i,
            epsilon = eps)))
        stdErrors[i, noNA] <- summary(temp.object)$coef[betasNames[noNA], "Std. Error"]
    }
    res <- sweep(stdErrors, 2, stdErrors[1, ], "/")
    colnames(res) <- names(coef(object))
    res
}


