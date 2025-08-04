#'
#' @examples
#'
#' ## A data set like that in Section 4.3 of
#' #https://doi.org/10.1016/j.spl.2023.109901
#' n <- 1000
#' p <- 400
#' gamma <- 3
#' X <- matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
#' betas0 <- rep(c(-3, -3/2, 0, 3/2, 3), each = p / 5)
#' betas <- gamma * betas0 / sqrt(sum(betas0^2))
#'
#' probs <- plogis(drop(X %*% betas))
#' y <- rbinom(n, 1, probs)
#'
#' fit_mdypl <- glm(y ~ -1 + X, family = binomial(), method = "mdyplFit")
#' plot(betas, type = "l", ylim = c(-0.5, 0.5))
#' points(coef(fit_mdypl), type ="l")
#'
mdyplFit <- function(x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
                     mustart = NULL, offset = rep(0, nobs), family = binomial(),
                     control = list(), intercept = TRUE,
                     fixed_totals = NULL, singular.ok = TRUE) {

    nobs <- NROW(y)
    if (!isTRUE(family$family == "binomial" && family$link == "logit")) {
        stop('`mdyplFit` only supports `binomial` family with `"logit"` link')
    }

    control <- do.call("mdyplControl", control)
    missing_offset <- is.null(offset)

    if (is.null(weights)) {
        weights <- rep.int(1, nobs)
    }

    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }

    if (missing_offset) {
        offset <- rep.int(0, nobs)
    }

    alpha <- unless_null(control$alpha, sum(weights) / (sum(weights) + ncol(x)))

    ## adjust responses as per MDYPL with beta_P = 0
    y_adj <- alpha * y + (1 - alpha) / 2

    out <- glm.fit(x = x, y = y_adj, weights = weights,
                   etastart = etastart, mustart = mustart,
                   offset = offset, family = family,
                   control = list(epsilon = control$epsilon,
                                  maxit = control$maxit, trace = control$trace),
                   intercept = intercept, singular.ok = singular.ok)

    mus <- out$fitted.values
    if (intercept & missing_offset) {
        nullmus <- mdyplFit(x = x[, "(Intercept)", drop = FALSE], y = y, weights = weights,
                            offset = rep(0, nobs), family = family, intercept = TRUE,
                            control = control,
                            start = family$linkfun(mean(y)))$fitted.values
    }

    if (!intercept) {
        nullmus <- family$linkinv(offset)
    }
    ## If there is an intercept and an offset then, for calculating
    ## the null deviance glm will make a call to the fitter to fit the
    ## glm with intercept and the offset
    if (intercept & !missing_offset) {
        nullmus <- mus
        ## doen't really matter what nullmus is set to. glm will make
        ## a new call to brglmFit and use the deviance from that call
        ## as null
    }

    ## Reset quantities in terms of original responses
    dev.resids <- family$dev.resids
    out$null.deviance <- sum(dev.resids(y, nullmus, weights))
    out$y_adj <- y_adj
    out$y <- y
    out$deviance <- sum(dev.resids(y, mus, weights))
    out$aic <- family$aic(y, n, mus, weights, deviance) + 2 * out$rank
    out$alpha <- alpha

    out$class <- c("mdyplFit", "brglmFit")
    out

}

#' Auxiliary function for [glm()] fitting using the [brglmFit()]
#' method.
#'
#' Typically only used internally by [brglmFit()], but may be used to
#' construct a `control` argument.
#'
#' @aliases mdypl_control
#' @param alpha the shrinkage parameter (in `[0, 1]`) in the
#'     Diaconis-Ylvisaker prior penalty. Default is \code{1}, which
#'     corresponds to maximum likelihood estimation. Setting `alpha =
#'     NULl` will result in `alpha = m / (m + p)`, where `m` is the
#'     sum of the binomial totals and `p` is the number of model
#'     parameters.
#' @param epsilon positive convergence tolerance epsilon. Default is
#'     `1e-08`.
#' @param maxit integer giving the maximal number of iterations
#'     allowed. Default is `25`.
#' @param trace logical indicating if output should be produced for
#'     each iteration. Default is `FALSE`.
#'
#' @export
mdyplControl <- function(alpha = NULL, epsilon = 1e-08, maxit = 25, trace = FALSE)
{
    out <- glm.control(epsilon, maxit, trace)
    if (!is.null(alpha)) {
        if (!(is.numeric(alpha)) || isTRUE(alpha < 0) || isTRUE(alpha > 1))
            stop("`alpha` should be in [0, 1]")
    }
    out$alpha <- alpha
    out
}

#' @export
sloe <- function(object, ...) {
    mu <- fitted(object)
    v <- mu * (1 - mu)
    h <- hatvalues(object)
    S <- object$linear.predictors - (object$y_adj - mu) / v * h / (1 - h)
    sd(S)
}


#' @export
correct.mdyplFit <- function(object, ...) {

}
