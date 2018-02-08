# Copyright (C) 2016, 2017 Ioannis Kosmidis

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


#' Bias reduction for multinomial response models using the
#' Poisson trick.
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
#' @param ref the reference category to use for multinomial
#'     regression. Either an integer, in which case
#'     levels(response)[ref] is used as a baseline, or a character
#'     string. Default is 1.
#'
#' @details
#'
#' The models \code{\link{brmultinom}} handles are also known as
#' baseline-category logit models (see, Agresti, 2002, Section 7.1),
#' because they model the log-odds of every category against a
#' baseline category. The user can control which baseline (or
#' reference) category is used via the \code{ref}. By default
#' \code{\link{brmultinom}} uses the first category as reference.
#'
#' The maximum likelihood estimates for the parameters of
#' baseline-category logit models have infinite components with
#' positive probability, which can result in problems in their
#' estimation and the use of inferential procedures (e.g. Wad
#' tests). Albert and Andreson (1984) have categorised the possible
#' data patterns for such models into the exclusive and exhaustive
#' categories of complete separation, quasi-complete separation and
#' overlap, and showed that infinite maximum likelihood estimates
#' result when complete or quasi-complete separation occurs.
#'
#' The adjusted score approach to bias reduction that
#' \code{\link{brmultinom}} implements (\code{type = "AS_mean"}) is an
#' alternative to maximum likelihood that results in estimates with
#' smaller asymptotic bias that are also *always* finite, even in
#' cases of complete or quasi-complete separation.
#'
#' \code{brmultinom} is a wrapper of \code{\link{brglmFit}} that fits
#' multinomial logit regression models through the 'Poisson trick' (see, for
#' example, Palmgren, 1981; Kosmidis & Firth, 2011).
#'
#' The implementation relies on the construction of an 'extended'
#' model matrix for the log-linear model and constraints on the sums
#' of the Poisson means. Specifically, a log-linear model is fitted on
#' a Kronecker product
#' (\url{https://en.wikipedia.org/wiki/Kronecker_product}) of the
#' original model matrix \code{X} implied by the formula, augmented by
#' \code{nrow(X)} dummy variables.
#'
#' The extended model matrix is sparse, and the \pkg{Matrix} package
#' is used for its effective storage.
#'
#' While \code{\link{brmultinom}} can be used for analyses using
#' multinomial regression models, the current implementation is more
#' of a 'proof of concept' and is not expected to scale well with
#' either of \code{nrow(X)}, \code{ncol(X)} or the number of levels in
#' the cateogrical response.
#'
#' @author Ioannis Kosmidis \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso \code{\link[nnet]{multinom}}
#'
#' @references
#'
#' Agrest A (2002). Categorical data analysis (2nd Edition). Wiley. New York.
#'
#' Albert A and Anderson J A (1984). On the Existence of Maximum
#' Likelihood Estimates in Logistic Regression Models. *Biometrika*,
#' **71** 1--10.
#'
#' Kosmidis I and Firth D (2011). Multinomial logit bias reduction via
#' the Poisson log-linear model. *Biometrika*, **98**, 755-759.
#'
#' Palmgren, J. (1981). The Fisher Information Matrix for Log Linear
#' Models Arguing Conditionally on Observed Explanatory
#' Variables. *Biometrika*, **68**, 563-566.
#'
#' @examples
#'
#' data("housing", package = "MASS")
#'
#' # Maximum likelihood using nnet::multinom
#' houseML1nnet <- nnet::multinom(Sat ~ Infl + Type + Cont, weights = Freq,
#'                                data = housing)
#' # Maximum likelihood using brmultinom with baseline category 'Low'
#' houseML1 <- brmultinom(Sat ~ Infl + Type + Cont, weights = Freq,
#'                        data = housing, type = "ML", ref = 1)
#' # The estimates are numerically the same as houseML0
#' all.equal(coef(houseML1nnet), coef(houseML1), tolerance = 1e-04)
#'
#' # Maximum likelihood using brmultinom with 'High' as baseline
#' houseML3 <- brmultinom(Sat ~ Infl + Type + Cont, weights = Freq,
#'                       data = housing, type = "ML", ref = 3)
#' # The fitted values are the same as houseML1
#' all.equal(fitted(houseML1), fitted(houseML1), tolerance = 1e-10)
#'
#' # Bias reduction
#' houseBR3 <- update(houseML3, type = "AS_mean")
#' # Bias correction
#' houseBC3 <- update(houseML3, type = "correction")
#'
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
        ## if (length(lev) == 2L)
        ##     Y <- as.integer(Y) - 1
        else Y <- nnet::class.ind(Y)
    }
    ##+END

    ncat <- if (is.matrix(Y)) ncol(Y) else length(lev)
    nvar <- ncol(X)

    if (is.character(ref)) {
        refc <- ref
        ref <- match(ref, lev, nomatch = NA)
        if (is.na(ref)) stop("reference category ", refc, " does not exist")
    }

    keep <- w > 0

    if (any(!keep)) {
        warning("Observations with non-positive weights have been omited from the computations")
    }

    nkeep <- sum(keep)
    ## Set up the model matrix for the poisson fit
    Xnuisance <- Matrix::Diagonal(nkeep)
    Xextended <- cbind(Matrix::kronecker(rep(1, ncat), Xnuisance),
                       Matrix::kronecker(Matrix::Diagonal(ncat)[, -ref, drop = FALSE], X[keep, ]))
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
                    family = poisson("log"), control = control, intercept = TRUE, fixed_totals = rep(seq.int(nkeep), ncat))

    ## TODO:
    ## + starting values
    ## + subset
    ## + na.action
    ## + control

    fit$call <- call
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
     if (is.null(coef.brmultinom(x))) {
         print("No coefficients")
     }
     else {
         print(format(coef.brmultinom(x), digits = digits), print.gap = 2, quote = FALSE)
     }
     cat("\nResidual Deviance:", format(x$deviance, digits = digits), "\n")
}

#' @method logLik brmultinom
#' @export
logLik.brmultinom <- function(object, ...) {
    structure(-object$deviance/2,
              df = sum(!is.na(coef.brmultinom(object))),
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

#' @method vcov brmultinom
#' @export
vcov.brmultinom <- function(object, ...) {
    vc <- vcov.brglmFit(object, ...)
    vc <- vc[object$ofInterest, object$ofInterest]
    vc
}


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
