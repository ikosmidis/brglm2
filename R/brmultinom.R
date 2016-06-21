#' Bias reduction for multinomial response models using the
#' "Poisson trick".
#'
#' \code{brmultinom} is a wrapper of \code{\link{brglmFit}} that fits
#' multinomial regression models using implicit and explicit bias
#' reduction methods. See Kosmidis & Firth (2011) for details.
#'
#' @inheritParams nnet::multinom
#' @param control a list of parameters for controlling the fitting
#'     process. See \code{\link{brglmControl}} for details.
#' @param ... arguments to be used to form the default 'control'
#'     argument if it is not supplied directly.
#'
#' @details
#' \code{brmultinom} is a wrapper of \code{\link{brglmFit}} that fits
#' multinomial regression models through the "Poisson trick" (see, for
#' example, Palmgren, 1981 and Kosmidis & Firth, 2009). The
#' implementation relies on the construction of an "extended" model
#' matrix for the log-linear model and constraints on the sums of the
#' Poisson means. Specifically, a log-linear model is fitted on a
#' Kronecker product
#' (\url{https://en.wikipedia.org/wiki/Kronecker_product}) of the
#' original model matrix \code{X} implied by the formula, augmented by
#' \code{nrow(X)} dummy variables.
#'
#' The extended model matrix is sparse, and the \pkg{Matrix} package
#' is used for its effective storage.
#'
#' While \code{\link{brmultinom}} can be used for analyses using
#' multinomial regression models, the current implementation is more
#' of a "proof of concept" and is not expected to scale well with
#' either of \code{nrow(X)}, \code{ncol(X)} or the number of levels in
#' the cateogrical response.
#'
#' @seealso \code{\link[nnet]{multinom}}
#'
#' @export
brmultinom <- function(formula, data, weights, subset, na.action, contrasts = NULL, ref = 1, control = list(...), ...) {
    call <- match.call()
    if (missing(data)) {
        data <- environment(formula)
    }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(model.frame)
    mf <- eval.parent(mf)
    Terms <- attr(mf, "terms")
    X <- model.matrix(Terms, mf, contrasts)
    Xcontrasts <- attr(X, "contrasts")
    Y <- model.response(mf, "any")

    ## The chunk of code between +BEGIN and +END has been adopted from
    ## nnet::multinom
    ##+BEGIN
    if (!is.matrix(Y)) {
        Y <- as.factor(Y)
        lev <- levels(Y)
    }
    else {
        lev <- colnames(Y)
    }
    w <- model.weights(mf)
    if (length(w) == 0L)
        if (is.matrix(Y))
            w <- rep(1, dim(Y)[1L])
        else w <- rep(1, length(Y))
    if (is.factor(Y)) {
        counts <- table(Y)
        if (any(counts == 0L)) {
            empty <- lev[counts == 0L]
            warning(sprintf(ngettext(length(empty), "group %s is empty",
                                     "groups %s are empty"), paste(sQuote(empty),
                                                                   collapse = " ")), domain = NA)
            Y <- factor(Y, levels = lev[counts > 0L])
            lev <- lev[counts > 0L]
        }
        if (length(lev) < 2L)
            stop("need two or more classes to fit a multinomial logit model")
        if (length(lev) == 2L)
            Y <- as.integer(Y) - 1
        else Y <- nnet::class.ind(Y)
    }
    ##+END

    ncat <- if (is.matrix(Y)) ncol(Y) else length(lev)
    nvar <- ncol(X)

    keep <- w > 0

    ## if (any(!keep)) {
    ##     warning("Observations with non-positive weights have been omited from the computations")
    ## }

    nkeep <- sum(keep)
    ## Set up the model matrix for the poisson fit
    Xnuisance <- Diagonal(nkeep)
    Xextended <- cbind(kronecker(rep(1, ncat), Xnuisance),
                       kronecker(Diagonal(ncat)[, -ref], X[keep, ]))
    int <- seq.int(nkeep)
    ## Set up the extended response
    Yextended <- c(Y[keep] * w[keep])

    nd <- paste0("%0", nchar(max(int)), "d")
    colnames(Xextended) <- c(paste0(".nuisance", sprintf(nd, int)),
                             ## CHECK: lev[-1] contrasts?
                             ofInterest <- paste(rep(lev[-ref], each = nvar),
                                                 rep(colnames(X), ncat - 1), sep = ":"))


    fit <- brglmFit(x = Xextended, y = Yextended,
                    start = NULL,
                    family = poisson("log"), control = control, intercept = TRUE, fixedTotals = rep(seq.int(nkeep), ncat))


    ## ## Set up the model matrix for the poisson fit
    ## Xnuisance <- Matrix::Diagonal(nrow(X))
    ## Xextended <- cbind(kronecker(rep(1, ncat), Xnuisance),
    ##                    kronecker(Matrix::Diagonal(ncat)[, -1], X))
    ## int <- seq.int(nrow(X))
    ## ## Set up the extended response
    ## Yextended <- c(Y * w)

    ## nd <- paste0("%0", nchar(max(int)), "d")
    ## colnames(Xextended) <- c(paste0(".nuisance", sprintf(nd, int)),
    ##                          ## CHECK: lev[-1] contrasts?
    ##                          paste(rep(lev[-1], each = nvar),
    ##                                rep(colnames(X), ncat - 1), sep = ":"))

    ## fixed <- rep(int, ncat)
    ## keep <- rep(w > 0, ncat)


    ## ## TODO:
    ## ## + starting values
    ## ## + subset
    ## ## + na.action
    ## ## + control
    ## fit <- brglmFit(x = Xextended[keep, ], y = Yextended[keep],
    ##                 start = NULL,
    ##                 family = poisson("log"), control = control, intercept = TRUE, fixedTotals = fixed[keep])

    fit$fitted.values <- matrix(fit$fitted.values, ncol = ncat)/w[keep]
    rownames(fit$fitted.values) <- rownames(X)[keep]
    colnames(fit$fitted.values) <- lev
    class(fit) <- c("brmultinom", fit$class, "glm")
    fit$ofInterest <- ofInterest
    fit$ncat <- ncat
    fit$lev <- lev
    fit$ref <- ref
    fit$coefNames <- colnames(X)
    fit
}

#' @method coef brmultinom
#' @export
coef.brmultinom <- function(object, ...) {
    ncat <- object$ncat
    if (length(object$ofInterest)) {
        coefficients <- object$coefficients[object$ofInterest]
        coefs <- matrix(coefficients, nrow = ncat - 1, byrow = TRUE)
        dimnames(coefs) <- with(object, list(lev[-object$ref], coefNames))
    }
    else {
        coefs <- NULL
    }
    coefs
}

#' @method print brmultinom
#' @export
print.brmultinom <- function(x, digits = max(5L, getOption("digits") - 3L), ...) {
     if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control = NULL)
     }
     cat("\nCoefficients:\n")
     if (is.null(coef(x))) {
         print("No coefficients")
     }
     else {
         print(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
     }
     cat("\nResidual Deviance:", format(x$deviance, digits = digits), "\n")
}

#' @method logLik brmultinom
#' @export
logLik.brmultinom <- function(object, ...) {
    structure(-object$deviance/2,
              df = sum(!is.na(coef(object))),
              nobs = sum(object$weights),
              class = "logLik")
}

#' @method summary brmultinom
#' @export
summary.brmultinom <- function (object, correlation = FALSE, digits = options()$digits,
                                Wald.ratios = FALSE, ...) {
    ncat <- object$ncat
    coefficients <- coef.brmultinom(object)
    object$digits <- digits
    object$AIC <- AIC(object)
    if (is.null(coefficients)) {
        object$coefficients <- NULL
        object$standard.errors <- NULL
        if (Wald.ratios)
            object$Wald.ratios <- NULL
        if (correlation)
            object$correlation <- NULL
    }
    else {
        vc <- vcov.brglmFit(object)
        vc <- vc[object$ofInterest, object$ofInterest]
        se <- sqrt(diag(vc))
        ses <- matrix(se, nrow = ncat - 1, byrow = TRUE, dimnames = dimnames(coefficients))
        object$coefficients <- coefficients
        object$standard.errors <- ses
        object$AIC <- AIC(object)
        if (Wald.ratios)
            object$Wald.ratios <- coef/ses
        if (correlation)
            object$correlation <- vc/outer(se, se)
    }
    class(object) <- "summary.brmultinom"
    object
}

#' Print method for \code{\link{summary.brmultinom}} objects
#' @section Note:
#' Code adopted from \code{nnet:::print.summary.brmultinom}
#' @method print summary.brmultinom
#' @export
print.summary.brmultinom <- function (x, digits = x$digits, ...)
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control = NULL)
    }
    cat("\nCoefficients:\n")
    ## if (x$is.binomial) {
    ##     print(cbind(Values = x$coefficients, `Std. Err.` = x$standard.errors,
    ##         `Value/SE` = x$Wald.ratios), digits = digits)
    ## }
    print(x$coefficients, digits = digits)
    cat("\nStd. Errors:\n")
    print(x$standard.errors, digits = digits)
    if (!is.null(x$Wald.ratios)) {
        cat("\nValue/SE (Wald statistics):\n")
        print(x$coefficients/x$standard.errors, digits = digits)
    }
    cat("\nResidual Deviance:", format(x$deviance), "\n")
    cat("AIC:", format(x$AIC), "\n")
    if (!is.null(correl <- x$correlation)) {
        p <- dim(correl)[2L]
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            ll <- lower.tri(correl)
            correl[ll] <- format(round(correl[ll], digits))
            correl[!ll] <- ""
            print(correl[-1L, -p], quote = FALSE, ...)
        }
    }
    invisible(x)
}
