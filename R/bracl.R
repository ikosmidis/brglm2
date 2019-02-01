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
#' COMPLETE ME
#'
#' @author Ioannis Kosmidis \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso \code{\link[nnet]{multinom}}, \code{\link[nnet]{brmultinom}}
#'
#' @references
#'
#' Kosmidis I., Kenne Pagui E. C. and Sartori N. (2018). Mean and
#' median bias reduction in generalized linear models. *arxiv*,
#' **arxiv:1804.04085**
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
    fit$null.deviance <- NULL
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
                coefs <- matrix(coefficients[ofInterest], ncol = ncat - 1)
                coefs <- -apply(cbind(coefs, 0), 1, diff)
                coefs <- c(coefs)
                names(coefs) <- c(sapply(coefNames, function(x) paste(object$lev[-ref], x, sep = ":")))
                coefs
            })
        }
    }
    else {
        NULL
    }
}

## Incomplete
vcov.bracl <- function(object, ...) {
    vc <- vcov.brglmFit(object, ...)
    coefNames <- names(coefficients(object))
    ofInterest <- object$ofInterest
    vc <- vc[ofInterest, ofInterest]
    levs <- object$lev[-object$ref]
    intercept_names <- paste0(levs, ":", "(Intercept)")
    ddiff <- function(mat) {
        mat <- diff(rbind(mat, 0))
        diff(rbind(t(mat), 0))
    }
    if (object$parallel) {
        beta_names <- ofInterest
        beta_names <- beta_names[!(beta_names %in% intercept_names)]
        vbeta <- vc[beta_names, beta_names]
        vint <- ddiff(vc[intercept_names, intercept_names])
        vintslo <- -diff(rbind(vc[intercept_names, beta_names], 0))
        par_names <- c(intercept_names, beta_names)
        vc[par_names, par_names] <- rbind(cbind(vint, vintslo),
                                          cbind(t(vintslo), vbeta))
    }
    else {
        betas <- object$coefNames
        npar <- length(betas)
        for (j in 1:npar) {
            for (k in 1:npar) {
                par_names1 <- paste(levs, betas[j], sep = ":")
                par_names2 <- paste(levs, betas[k], sep = ":")
                mat <- ddiff(vc[par_names1, par_names2])
                vc[par_names1, par_names2] <- mat
            }
        }
        ## re-order
        vc <- vc[coefNames, coefNames]
    }
    vc
}

summary.bracl <- function(object, correlation = FALSE, digits = 3, ...) {
    object$digits <- digits
    object$logLik <- logLik(object)
    coefs <- coefficients(object)
    aliased <- is.na(coefs)
    vc <- vcov(object)
    var_coef <- diag(vc)
    s.err <- sqrt(var_coef)
    tvalue <- coefs/s.err
    dn <- c("Estimate", "Std. Error")
    pvalue <- 2 * pnorm(-abs(tvalue))
    coef_table <- cbind(coefs, s.err, tvalue, pvalue)
    dimnames(coef_table) <- list(names(coefs), c(dn, "z value", "Pr(>|z|)"))
    object$coefficients <- coef_table
    if (correlation) {
        object$correlation <- vc/outer(s.err, s.err)
    }
    class(object) <- "summary.bracl"
    return(object)
}

print.summary.bracl <- function(x, digits = x$digits, ...) {
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control = NULL)
    }
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, digits = digits)
    cat("\nResidual Deviance:", format(x$deviance), "\n")
    cat("Log-likelihood:", format(x$logLik), "\n")
    cat("AIC:", format(x$aic), "\n")
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

