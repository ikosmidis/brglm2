#' Bias reduction for multinomial response models using the
#' "Poisson trick". See Kosmidis & Firth (2011) for details
#'
#' @details
#'
#' This function is a wrapper for brglmFit.R that can be used to get
#' bias reduced estimates for the parameters of multinomial regression
#' models. The implementation is based on constructing an appropriate
#' model matrix by taking a Kronecker product
#' (\url{https://en.wikipedia.org/wiki/Kronecker_product}) of the
#' model matrix X implied by the formula and augmenting it with
#' \code{nrow(X)} dummy variables, to indi
#'
#' uses the \pkg{Matrix} package to
#'
#'
#' @export
brmultinom <- function(formula, data, weights, subset, na.action, contrasts = NULL, control = list(...), ...) {
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
        ## Consider removing the below
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

    if (any(!keep)) {
        warning("Observations with non-positive weights have been omited from the computations")
    }

    nkeep <- sum(keep)
    ## Set up the model matrix for the poisson fit
    Xnuisance <- Matrix::Diagonal(nkeep)
    Xextended <- cbind(kronecker(rep(1, ncat), Xnuisance),
                       kronecker(Matrix::Diagonal(ncat)[, -1], X[keep, ]))
    int <- seq.int(nkeep)
    ## Set up the extended response
    Yextended <- c(Y[keep] * w[keep])

    nd <- paste0("%0", nchar(max(int)), "d")
    colnames(Xextended) <- c(paste0(".nuisance", sprintf(nd, int)),
                             ## CHECK: lev[-1] contrasts?
                             ofInterest <- paste(rep(lev[-1], each = nvar),
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
    fit$coefNames <- colnames(X)
    fit
}

coef.brmultinom <- function(object, ...) {
    coefs <- with(object, matrix(coefficients[ofInterest], nrow = ncat - 1, byrow = TRUE))
    dimnames(coefs) <- with(object, list(lev[-1], coefNames))
    coefs

}

print.brmultinom <- function(x, digits = max(5L, getOption("digits") - 3L), ...) {
     if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control = NULL)
    }
    cat("\nCoefficients:\n")
    print(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
    cat("\nResidual Deviance:", format(x$deviance, digits = digits), "\n")
}

logLik.brmultinom <- function(object, ...) {
    structure(-object$deviance/2,
              df = sum(!is.na(coef(object))),
              nobs = sum(object$weights),
              class = "logLik")
}
