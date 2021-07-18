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
#' @param ref the reference category to use for multinomial
#'     regression. Either an integer, in which case
#'     levels(response)[ref] is used as a baseline, or a character
#'     string. Default is 1.
#' @param x should the model matrix be included with in the result
#'     (default is \code{TRUE}).
#' @param ... arguments to be used to form the default 'control'
#'     argument if it is not supplied directly.
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
#' tests). Albert and Andreson (1984) have categorized the possible
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
#' the categorical response.
#'
#' @author Ioannis Kosmidis \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @seealso \code{\link[nnet]{multinom}}, \code{\link{bracl}} for adjacent category logit models for ordinal responses
#'
#' @references
#'
#' Kosmidis I, Kenne Pagui E C, Sartori N (2020). Mean and median bias
#' reduction in generalized linear models. *Statistics and Computing*,
#' **30**, 43-59 \doi{10.1007/s11222-019-09860-6}
#'
#' Agresti A (2002). *Categorical data analysis* (2nd edition). Wiley
#' Series in Probability and Statistics. Wiley.
#'
#' Albert A, Anderson J A (1984). On the Existence of Maximum
#' Likelihood Estimates in Logistic Regression Models. *Biometrika*,
#' **71** 1--10, \doi{10.2307/2336390}
#'
#' Kosmidis I, Firth D (2011). Multinomial logit bias reduction
#' via the Poisson log-linear model. *Biometrika*, **98**, 755-759
#' \doi{10.1093/biomet/asr026}
#'
#' Palmgren, J (1981). The Fisher Information Matrix for Log Linear
#' Models Arguing Conditionally on Observed Explanatory
#' Variables. *Biometrika*, **68**, 563-566
#' \doi{10.1093/biomet/68.2.563}
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
#' all.equal(fitted(houseML3), fitted(houseML1), tolerance = 1e-10)
#'
#' # Bias reduction
#' houseBR3 <- update(houseML3, type = "AS_mean")
#' # Bias correction
#' houseBC3 <- update(houseML3, type = "correction")
#'
#'
#' @export
brmultinom <- function(formula, data, weights, subset, na.action,
                       contrasts = NULL, ref = 1,
                       model = TRUE, x = TRUE,
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
    nkeep <- sum(keep)
    fixed_totals <- rep(seq.int(nkeep), ncat)

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
                    family = poisson("log"), control = control, intercept = TRUE, fixed_totals = fixed_totals)

    ## TODO:
    ## + starting values
    ## + subset
    ## + na.action
    ## + control

    ## Fitted values
    fitted <- do.call("rbind", tapply(fit$fitted, fixed_totals, function(x) x/sum(x)))
    rownames(fitted) <- rownames(X)[keep]
    colnames(fitted) <- lev
    fit$fitted.values <- fitted
    fit$call <- call
    ## fit$fitted.values <- matrix(fit$fitted.values, ncol = ncat)/w[keep]
    ## rownames(fit$fitted.values) <- rownames(X)[keep]
    ## colnames(fit$fitted.values) <- lev

    class(fit) <- c("brmultinom", fit$class, "glm")
    fit$ofInterest <- ofInterest
    fit$ncat <- ncat
    fit$lev <- lev
    fit$ref <- ref
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

#' @method fitted brmultinom
#' @export
fitted.brmultinom <- function(object, ...) {
    object$fitted.values
}

#' Residuals for multinomial logistic regression and adjacent category logit models
#'
#' @param object the object coming out of \code{\link{bracl}} and
#'     \code{\link{brmultinom}}.
#' @param type the type of residuals which should be returned.  The
#'     options are: \code{"pearson"} (default), \code{"response"},
#'     \code{"deviance"}, \code{"working"}. See Details.
#' @param ... Currently not used.
#'
#' @details
#'
#' The residuals computed are the residuals from the equivalent
#' Poisson log-linear model fit, organized in a form that matches the
#' output of \code{fitted(object, type = "probs")}. As a result, the
#' output is residuals defined in terms of the object and expected
#' multinomial counts.
#'
#' @seealso brmultinom bracl
#'
#' @method residuals brmultinom
#' @export
residuals.brmultinom <- function(object, type = c("pearson", "response", "deviance", "working"), ...) {
    type <- match.arg(type)
    ## This is a Poisson log-linear models, so the working weights are
    ## the fitted counts
    fitted <- weights(object, type = "working")
    ## The poisson responses
    y <- object$y
    out <- switch(type,
                  "pearson" = (y - fitted)/sqrt(fitted),
                  "response" = (y - fitted),
                  "working" = object$residuals,
                  "deviance" = object$family$dev.resids(y, fitted, 1))
    matrix(out, ncol = object$ncat, dimnames = dimnames(fitted(object)))
}

#' @method coef brmultinom
#' @export
coef.brmultinom <- function(object, ...) {
    if (length(object$ofInterest)) {
        with(object, {
            coefs <- matrix(coefficients[ofInterest], nrow = ncat - 1, byrow = TRUE)
            dimnames(coefs) <- list(lev[-object$ref], coefNames)
            coefs
        })
    }
    else {
        NULL
    }
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
     cat("\nResidual Deviance:", format(x$deviance), "\n")
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
summary.brmultinom <- function(object, correlation = FALSE, digits = options()$digits, Wald.ratios = FALSE, ...) {
    ncat <- object$ncat
    coefficients <- coef.brmultinom(object)
    object$digits <- digits
    object$AIC <- AIC(object)
    object$logLik <- logLik(object)
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
        ## object$AIC <- AIC(object)
        if (Wald.ratios) {
            object$Wald.ratios <- coefficients/ses
            object$Wald.pvalues <-  2 * pnorm(-abs(object$Wald.ratios))
        }
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
print.summary.brmultinom <- function(x, digits = x$digits, ...)
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control = NULL)
    }
    cat("\nCoefficients:\n")
    print(x$coefficients, digits = digits)
    cat("\nStd. Errors:\n")
    print(x$standard.errors, digits = digits)
    if (!is.null(x$Wald.ratios)) {
        cat("\nValue/SE (Wald statistics):\n")
        print(x$coefficients/x$standard.errors, digits = digits)
    }
    cat("\nResidual Deviance:", format(x$deviance), "\n")
    cat("Log-likelihood:", format(x$logLik), "\n")
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

#' Predict method for \code{brmultinom} fits
#'
#' Obtain class and probability predictions from a fitted baseline
#' category logits model.
#'
#' @param object a fitted object of class inheriting from
#'     \code{"brmultinom"}.
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
#' data("housing", package = "MASS")
#'
#' # Maximum likelihood using brmultinom with baseline category 'Low'
#' houseML1 <- brmultinom(Sat ~ Infl + Type + Cont, weights = Freq,
#'                        data = housing, type = "ML", ref = 1)
#'
#' # New data
#' newdata <- expand.grid(Infl = c("Low", "Medium"),
#'                        Type = c("Tower", "Atrium", "Terrace"),
#'                        Cont = c("Low", NA, "High"))
#'
#' ## Predictions
#' sapply(c("class", "probs"), function(what) predict(houseML1, newdata, what))
#'
#' @method predict brmultinom
#' @export
predict.brmultinom <- function(object, newdata, type = c("class", "probs"), ...)
{
    ## Adapted from nnet:::predict.multinom
    if (!inherits(object, "brmultinom"))
        stop("not a \"brmultinom\" fit")
    type <- match.arg(type)
    X <- if (missing(newdata)) model.matrix(object) else model.matrix(object, data = newdata)
    rn <- attr(X, "rn_data")
    keep <- attr(X, "rn_kept")
    coefs <- coef(object)
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

#' @method model.matrix brmultinom
#' @export
model.matrix.brmultinom <- function(object, data, ...) {
    if (!inherits(object, "brmultinom"))
        stop("not a \"brmultinom\" fit")
    if (missing(data)) {
        data <- model.frame(object)
    }
    else {
        data <- as.data.frame(data)
    }
    Terms <- delete.response(object$terms)
    m <- model.frame(Terms, data, na.action = na.omit,
                     xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses")))
        .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts = object$contrasts)
    rn <- row.names(data)
    attr(X, "rn_data") <- rn
    attr(X, "rn_kept") <-  match(row.names(m), rn)
    X
}


#' Method for computing confidence intervals for one or more
#' regression parameters in a \code{\link{brmultinom}} object
#'
#' @inheritParams stats::confint
#'
#' @export
confint.brmultinom <- function (object, parm, level = 0.95, ...)  {
    ## Apart from formatting changes this function is identical to
    ## nnet:::confint.multinom
    cf <- coef(object)
    pnames <- if (is.matrix(cf)) colnames(cf) else names(cf)
    if (missing(parm)) {
        parm <- seq_along(pnames)
    }
    else {
        if (is.character(parm))  {
            parm <- match(parm, pnames, nomatch = 0L)
        }
    }
    a <- (1 - level) / 2
    a <- c(a, 1 - a)
    pct <- paste(round(100 * a, 1), "%")
    fac <- qnorm(a)
    if (is.matrix(cf)) {
        ses <- matrix(sqrt(diag(vcov(object))), ncol = ncol(cf),
            byrow = TRUE)[, parm, drop = FALSE]
        cf <- cf[, parm, drop = FALSE]
        ci <- array(NA, dim = c(dim(cf), 2L), dimnames = c(dimnames(cf),
            list(pct)))
        ci[, , 1L] <- cf + ses * fac[1L]
        ci[, , 2L] <- cf + ses * fac[2L]
        aperm(ci, c(2L, 3L, 1L))
    }
    else {
        ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(pnames[parm],
            pct))
        ses <- sqrt(diag(vcov(object)))[parm]
        ci[] <- cf[parm] + ses %o% fac
        ci
    }
}


#' Method for simulating a data set from \code{\link{brmultinom}} and
#' \code{\link{bracl}} objects
#'
#' @param object an object of class \code{\link{brmultinom}} or
#'     \code{\link{bracl}}.
#' @param ... currently not used.
#'
#' @return
#'
#' A \code{\link{data.frame}} with \code{object$ncat} times the rows
#' that \code{model.frame(object)} have and the same variables. If
#' \code{weights} has been specified in the call that generated
#' \code{object}, then the simulate frequencies will populate the
#' weights variable. Otherwise, the resulting \code{data.frame} will
#' have a \code{".weights"} variable with the simulated multinomial
#' counts.
#'
#' @examples
#'
#' ## Multinomial logistic regression
#' data("housing", package = "MASS")
#' houseML1 <- brmultinom(Sat ~ Infl + Type + Cont, weights = Freq,
#'                        data = housing, type = "ML", ref = 1)
#' simulate(houseML1)
#'
#' ## Adjacent-category logits
#' data("stemcell", package = "brglm2")
#' stemML1 <- bracl(research ~ religion + gender, weights = frequency,
#'                 data = stemcell, type = "ML")
#'
#' simulate(stemML1)
#'
#' @export
simulate.brmultinom <- function(object, ...) {
    mf <- model.frame(object)
    probs <- predict(object, type = "probs")
    categories <- colnames(probs)
    ncat <- object$ncat
    weights <- model.weights(mf)
    if (is.null(weights)) {
        weights <- rep.int(1L, nrow(mf))
    }
    samples <- sapply(1:nrow(probs), function(j) rmultinom(1, weights[j], probs[j, ]))
    mf <- mf[rep(1:nrow(mf), each = ncat), ]
    mf[, 1] <- factor(colnames(probs),
                      levels = levels(mf[, 1]),
                      ordered = is.ordered(mf[, 1]))
    weights_ind <- grep("(weights)", names(mf))
    if (length(weights_ind)) {
        weights_nam <- as.character(object$call$weights)
        names(mf)[weights_ind] <- weights_nam
    }
    else {
        weights_nam <- ".weights"
    }
    mf[[weights_nam]] <- c(samples)
    mf
}
