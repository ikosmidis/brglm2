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


#' Bias reduction for adjcacent category logit models using the
#' Poisson trick.
#'
#' \code{bracl} is a wrapper of \code{\link{brglmFit}} that fits
#' adjacent category logit models using implicit and explicit bias
#' reduction methods. See Kosmidis & Firth (2011) for details.
#'
#' @inheritParams nnet::multinom
#' @param control a list of parameters for controlling the fitting
#'     process. See \code{\link{brglmControl}} for details.
#' @param ... arguments to be used to form the default 'control'
#'     argument if it is not supplied directly.
#'
#' @details
#'
#' Complete me
#'
#' @author Ioannis Kosmidis \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso \code{\link[nnet]{multinom}}, \code{\link[nnet]{brmultinom}}
#'
#' @references
#'
#' Agresti A. (2002). Categorical data analysis (2nd
#' Edition). Wiley. New York.
#'
#' Albert A. and Anderson J. A. (1984). On the Existence of Maximum
#' Likelihood Estimates in Logistic Regression Models. *Biometrika*,
#' **71** 1--10.
#'
#' Kosmidis I. and Firth D. (2011). Multinomial logit bias reduction via
#' the Poisson log-linear model. *Biometrika*, **98**, 755-759.
#'
#' Palmgren, J. (1981). The Fisher Information Matrix for Log Linear
#' Models Arguing Conditionally on Observed Explanatory
#' Variables. *Biometrika*, **68**, 563-566.
#'
#' @examples
#'
#'
#' @export
bracl <- function(formula, data, weights, subset, na.action,
                  parallel = FALSE,
                  contrasts = NULL,
                  model = TRUE, x= TRUE,
                  control = list(...), ...) {
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
    ## Model matrix
    X <- model.matrix(Terms, mf, contrasts)
    xint <- match("(Intercept)", colnames(X), nomatch = 0L)
    ## If there is no intercept include it and issue warning
    if (xint == 0L) {
        X <- cbind(1L, X)
        xint <- 1
        colanmes(X) <- c("(Intercept)", colnames(X))
        warning("an intercept is needed and assumed")
    }
    Xcontrasts <- attr(X, "contrasts")
    ## Response
    Y <- model.response(mf, "any")
    n <- length(Y)
    if (!is.factor(Y)) {
        stop("response must be a factor")
    }
    lev <- levels(Y)
    ncat <- length(lev)
    if (ncat <= 2L) {
        stop("response must have 3 or more levels")
    }
    Y <- nnet::class.ind(Y)
    cons <- attr(X, "contrasts")
    ## Weights
    w <- model.weights(mf)
    if (!length(w)) {
        w <- rep(1, n)
    }

    ## Offset
    offset <- model.offset(mf)
    if (length(offset) <= 1L) {
        offset <- rep(0, n)
    }

    nvar <- ncol(X)
    keep <- w > 0
    nkeep <- sum(keep)

    cats0 <- apply(Y, 1, function(x) which(x == 1))
    cats <- rep(seq.int(ncat), each = nrow(X))

    fixed_totals <- rep(seq.int(nkeep), ncat)

    ## Set up the model matrix for the poisson fit
    Xnuisance <- Matrix::Diagonal(nkeep)
    int <- seq.int(nkeep)
    nd <- paste0("%0", nchar(max(int)), "d")
    if (parallel) {
        Xextended <- cbind(Matrix::kronecker(rep(1, ncat), Xnuisance),
                           Matrix::kronecker(Matrix::Diagonal(ncat)[, -ncat, drop = FALSE], X[keep, xint]),
        (ncat - cats) * Matrix::kronecker(c(rep(1, ncat - 1), 0), X[keep, -xint]))
        ofInterest <- c(paste(lev[-ncat], rep("(Intercept)", ncat - 1), sep = ":"), colnames(X)[-xint])
    }
    else {
        Xextended <- cbind(Matrix::kronecker(rep(1, ncat), Xnuisance),
                           Matrix::kronecker(Matrix::Diagonal(ncat)[, -ncat, drop = FALSE], X[keep, ]))
        ofInterest <- paste(rep(lev[-ncat], each = nvar),
                            rep(colnames(X), ncat - 1), sep = ":")
    }
    colnames(Xextended) <- c(paste0(".nuisance", sprintf(nd, int)),
                             ofInterest)

    ## Set up the extended response
    Yextended <- c(Y[keep] * w[keep])



    fit <- brglmFit(x = Xextended, y = Yextended,
                    start = NULL,
                    family = poisson("log"), control = control, intercept = TRUE, fixed_totals = fixed_totals)

    ## Fitted values
    fitted <- do.call("rbind", tapply(fit$fitted, fixed_totals, function(x) x/sum(x)))
    rownames(fitted) <- rownames(X)
    colnames(fitted) <- lev
    fit$fitted.values <- fitted
    fit$parallel  <- parallel
    fit$call <- call
    class(fit) <- c("bracl", "brmultinom", fit$class, "glm")
    fit$ofInterest <- ofInterest
    fit$ncat <- ncat
    fit$lev <- lev
    fit$ref <- ncat
    if (model) {
        fit$model  <- mf
    }
    if (x) {
        fit$x  <- X
    }
    fit$contrasts <- attr(X, "contrasts")
    fit$xlevels = .getXlevels(Terms, mf)
    fit$terms <- Terms
    fit$coefNames <- colnames(X)
    fit
}

fitted.bracl <- function(object, ...) {
    object$fitted.values
}

coef.bracl <- function(object, ...) {
    if (length(object$ofInterest)) {
        if (object$parallel) {
            with(object, {
                coefs <- coefficients[ofInterest]
                intercept_names <- paste0(lev[-ref], ":", "(Intercept)")
                coefs[intercept_names] <- -diff(c(coefs[intercept_names], 0))
                names(coefs[intercept_names]) <- intercept_names
                coefs
            })
        }
        else {
            with(object, {
                coefs <- matrix(coefficients[ofInterest], nrow = ncat - 1, byrow = TRUE)
                coefs <- -apply(rbind(coefs, 0), 2, diff)
                dimnames(coefs) <- list(lev[-ref], coefNames)
                coefs
            })
        }
    }
    else {
        NULL
    }
}

vcov.bracl <- function(object, ...) {
    vc <- vcov.brglmFit(object, ...)
    vc <- vc[object$ofInterest, object$ofInterest]
    if (object$parallel) {

    }
    else {

    }
}

summary.bracl <- function (object, correlation = FALSE, digits = options()$digits,
                                Wald.ratios = FALSE, ...) {

}
