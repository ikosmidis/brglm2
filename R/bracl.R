# Copyright (C) 2016-2021 Ioannis Kosmidis

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


#' Bias reduction for adjacent category logit models for ordinal
#' responses using the Poisson trick.
#'
#' \code{bracl} is a wrapper of \code{\link{brglmFit}} that fits
#' adjacent category logit models with or without proportional odds
#' using implicit and explicit bias reduction methods. See Kosmidis &
#' Firth (2011) for details.
#'
#' @inheritParams MASS::polr
#' @param control a list of parameters for controlling the fitting
#'     process. See \code{\link{brglmControl}} for details.
#' @param parallel if \code{FALSE} (default), then a non-proportional
#'     odds adjacent category model is fit, assuming different
#'     effects per category; if \code{TRUE} then a proportional odds
#'     adjacent category model is fit. See Details.
#' @param x should the model matrix be included with in the result
#'     (default is \code{TRUE}).
#' @param ... arguments to be used to form the default 'control'
#'     argument if it is not supplied directly.
#'
#' @details
#'
#' The \code{bracl} function fits adjacent category models, which
#' assume multinomial observations with probabilities with
#' proportional odds of the form
#'
#' \deqn{\log\frac{\pi_{ij}}{\pi_{ij + 1}} = \alpha_j + \beta^T x_i}{log(pi[i, j]/pi[i, j+1]) = alpha[j] + sum(beta * x[i, ])}
#'
#' or with non-proportional odds of the form
#'
#' \deqn{\log\frac{\pi_{ij}}{\pi_{ij + 1}} = \alpha_j + \beta_j^T x_i}{log(pi[i, j]/pi[i, j+1]) = alpha[j] + sum(beta[j, ] * x[i, ])}
#'
#' where \eqn{x_i}{x[i, ]} is a vector of covariates and \eqn{\pi_{ij}}{pi[i, j]} is the
#' probability that category \eqn{j} is observed at the covariate setting \eqn{i}.
#'
#' @author Ioannis Kosmidis \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso \code{\link[nnet]{multinom}}, \code{\link{brmultinom}}
#'
#' @references
#'
#' Kosmidis I, Kenne Pagui E C, Sartori N (2020). Mean and median bias
#' reduction in generalized linear models. *Statistics and Computing*,
#' **30**, 43-59 \doi{10.1007/s11222-019-09860-6}
#'
#' Agresti, A (2010). *Analysis of Ordinal Categorical Data* (2nd
#' edition).  Wiley Series in Probability and Statistics. Wiley.
#'
#' Albert A, Anderson J A (1984). On the Existence of Maximum
#' Likelihood Estimates in Logistic Regression Models. *Biometrika*,
#' **71**, 1--10 \doi{10.2307/2336390}
#'
#' Kosmidis I, Firth D (2011). Multinomial logit bias reduction
#' via the Poisson log-linear model. *Biometrika*, **98**,
#' 755-759 \doi{10.1093/biomet/asr026}
#'
#' Palmgren J (1981). The Fisher Information Matrix for Log Linear
#' Models Arguing Conditionally on Observed Explanatory
#' Variables. *Biometrika*, **68**,
#' 563-566 \doi{10.1093/biomet/68.2.563}
#'
#' @examples
#'
#' data("stemcell", package = "brglm2")
#'
#' # Adjacent category logit (non-proportional odds)
#' fit_bracl <- bracl(research ~ as.numeric(religion) + gender, weights = frequency,
#'                    data = stemcell, type = "ML")
#' # Adjacent category logit (proportional odds)
#' fit_bracl_p <- bracl(research ~ as.numeric(religion) + gender, weights = frequency,
#'                     data = stemcell, type = "ML", parallel = TRUE)
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
        cn <- colnames(X)
        X <- cbind(1L, X)
        xint <- 1
        colnames(X) <- c("(Intercept)", cn)
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

    ## cats0 <- apply(Y, 1, function(x) which(x == 1))
    cats <- rep(seq.int(ncat), each = nrow(X))[keep]

    fixed_totals <- rep(seq.int(nkeep), ncat)

    ## Set up the model matrix for the poisson fit
    Xnuisance <- Matrix::Diagonal(nkeep)
    int <- seq.int(nkeep)
    nd <- paste0("%0", nchar(max(int)), "d")
    if (parallel) {
        Xextended <- cbind(Matrix::kronecker(rep(1, ncat), Xnuisance),
                           Matrix::kronecker(Matrix::Diagonal(ncat)[, -ncat, drop = FALSE], X[keep, xint]),
        (ncat - cats) * Matrix::kronecker(c(rep(1, ncat - 1), 0), X[keep, -xint, drop = FALSE]))
        ofInterest <- c(paste(lev[-ncat], rep("(Intercept)", ncat - 1), sep = ":"), colnames(X)[-xint])
    }
    else {
        Xextended <- cbind(Matrix::kronecker(rep(1, ncat), Xnuisance),
                           Matrix::kronecker(Matrix::Diagonal(ncat)[, -ncat, drop = FALSE], X[keep, , drop = FALSE]))
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
    rownames(fitted) <- rownames(X)[keep]
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

#' @method fitted bracl
#' @export
fitted.bracl <- function(object, ...) {
    object$fitted.values
}

#' @method coef bracl
#' @export
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

#' @method vcov bracl
#' @export
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
        if (nrow(vintslo) == 1) {
            vintslo <- drop(vintslo)
        }
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

#' @method summary bracl
#' @export
summary.bracl <- function(object, correlation = FALSE, digits = 3, ...) {
    object$digits <- digits
    object$logLik <- logLik(object)
    object$AIC <- AIC(object)
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

#' @method print summary.bracl
#' @export
print.summary.bracl <- function(x, digits = x$digits, ...) {
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control = NULL)
    }
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, digits = digits)
    cat("\nResidual Deviance:", format(x$deviance), "\n")
    cat("Log-likelihood:", format(x$logLik), "\n")
    cat("AIC:", format(x$AIC), "\n")
    cat("\n\nType of estimator:", x$type, get_type_description(x$type))
    cat("\n", "Number of Fisher Scoring iterations: ", x$iter, "\n", sep = "")
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

#' Predict method for \code{bracl} fits
#'
#' Obtain class and probability predictions from a fitted adjacent
#' category logits model.
#'
#' @param object a fitted object of class inheriting from
#'     \code{"bracl"}.
#' @param newdata optionally, a data frame in which to look for
#'     variables with which to predict.  If omitted, the fitted linear
#'     predictors are used.
#' @param type the type of prediction required. The default is
#'     \code{"class"}, which produces predictions of the response
#'     category at the covariate values supplied in \code{"newdata"},
#'     selecting the category with the largest probability; the
#'     alternative \code{"probs"} returns all category probabilities
#'     at the covariate values supplied in \code{"newdata"}.
#' @param ... further arguments passed to or from other methods.
#'
#'
#' @details
#'
#' If \code{newdata} is omitted the predictions are based on the data
#' used for the fit.
#'
#' @return
#'
#' If \code{type = "class"} a vector with the predicted response
#' categories; if \code{type = "probs"} a matrix of probabilities for
#' all response categories at \code{newdata}.
#'
#' @examples
#'
#' data("stemcell", package = "brglm2")
#'
#' # Adjacent category logit (non-proportional odds)
#' fit_bracl <- bracl(research ~ as.numeric(religion) + gender, weights = frequency,
#'                    data = stemcell, type = "ML")
#' # Adjacent category logit (proportional odds)
#' fit_bracl_p <- bracl(research ~ as.numeric(religion) + gender, weights = frequency,
#'                     data = stemcell, type = "ML", parallel = TRUE)
#'
#' # New data
#' newdata <- expand.grid(gender = c("male", "female"),
#'                        religion = c("liberal", "moderate", "fundamentalist"))
#'
#' # Predictions
#' sapply(c("class", "probs"), function(what) predict(fit_bracl, newdata, what))
#' sapply(c("class", "probs"), function(what) predict(fit_bracl_p, newdata, what))
#'
#' @method predict bracl
#' @export
predict.bracl <- function(object, newdata, type = c("class", "probs"), ...) {
    ## Adapted from nnet:::predict.multinom
    if (!inherits(object, "bracl"))
        stop("not a \"bracl\" fit")
    type <- match.arg(type)
    X <- if (missing(newdata)) model.matrix(object) else model.matrix(object, data = newdata)
    rn <- attr(X, "rn_data")
    keep <- attr(X, "rn_kept")
    cc <- coef(object)
    ## Ignore unidentifiable parameters
    cc[is.na(cc)] <- 0
    nams <- names(cc)
    if (object$parallel) {
        int <- (object$ncat - 1):1
        sl <- nams[-int]
        coefs <- cbind(rev(cumsum(cc[int])),
                       int * matrix(cc[sl], nrow = object$ncat - 1, ncol = length(sl), byrow = TRUE))
        rownames(coefs) <- object$lev[-object$ref]
    }
    else {
        coefs <- matrix(cc, nrow = object$ncat - 1)
        rownames(coefs) <- object$lev[-object$ref]
        coefs <- apply(coefs, 2, function(x) cumsum(rev(x)))
    }

    fits <- matrix(0, nrow = nrow(X), ncol = object$ncat, dimnames = list(rn[keep], object$lev))
    fits1 <- apply(coefs, 1, function(b) X %*% b)
    fits[, rownames(coefs)] <- fits1
    Y1 <- t(apply(fits, 1, function(x) exp(x) / sum(exp(x))))
    Y <- matrix(NA, length(rn), ncol(Y1), dimnames = list(rn, colnames(Y1)))
    Y[keep, ] <- Y1
    switch(type, class = {
        if (length(object$lev) > 2L) Y <- factor(max.col(Y),
            levels = seq_along(object$lev), labels = object$lev)
        if (length(object$lev) == 2L) Y <- factor(1 + (Y > 0.5),
            levels = 1L:2L, labels = object$lev)
        if (length(object$lev) == 0L) Y <- factor(max.col(Y),
            levels = seq_along(object$lab), labels = object$lab)
    }, probs = {
    })
    drop(Y)
}
