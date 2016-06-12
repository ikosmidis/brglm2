## Suggestion by Kurt Hornik to avoid a warning related to the binding
## of n which is evaluated by family$initialize
if(getRversion() >= "2.15.1") globalVariables("n")

#' Bias reduction in Binomial-response GLMs
#'
#' Fits binomial-response GLMs using the bias-reduction method developed in
#' Firth (1993) for the removal of the leading (\eqn{\mathop{\rm
#' O}(n^{-1})}{O(n^{-1})}) term from the asymptotic expansion of the bias
#' of the maximum likelihood estimator. Fitting is performed using
#' pseudo-data representations, as described in Kosmidis (2007, Chapter 5). For
#' estimation in binomial-response GLMs, the bias-reduction method is an
#' improvement over traditional maximum likelihood because:
#' \itemize{
#' \item the bias-reduced estimator is second-order unbiased and has smaller
#'       variance than the maximum likelihood estimator and
#' \item the resultant estimates and their corresponding standard errors
#'       are \bold{always} finite while the maximum likelihood estimates
#'       can be infinite (in situations where complete or quasi separation
#'       occurs).
#' }
#'
#' @aliases brglm0 brglm0.fit print.brglm0 summary.brglm0 print.summary.brglm0
#' @useDynLib brglm2
#' @export
brglm0 <- function (formula, family = binomial, data, weights, subset,
    na.action, start = NULL, etastart, mustart, offset, control.glm = glm.control1(...),
    model = TRUE, method = "brglm.fit", pl = FALSE, x = FALSE,
    y = TRUE, contrasts = NULL, control.brglm = brglm0.control(...),
    ...)
{
    call <- match.call()
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    br <- method == "brglm.fit"
    ####################
    ## More families to be implemented
    if (br & family$family != "binomial")
        stop("families other than 'binomial' are not supported by brglm0. Use the brglmFit fitting method instead")
    ####################
    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    switch(method, model.frame = return(mf), glm.fit = fit.proc <- glm.fit,
        brglm.fit = fit.proc <- brglm0.fit, stop("invalid 'method': ",
            method))
    ####################
    ## Arg control of fit.proc
    if (br) {
        formals(fit.proc)$control.brglm <- control.brglm
    }
    if (pl)
        formals(fit.proc)$pl <- TRUE
    ####################
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    Xor <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0)
    Xmax <- apply(abs(Xor), 2, max)
    Xmax[Xmax==0] <- 1
    X <- sweep(Xor, 2, Xmax, "/")
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")
    if (!is.null(offset)) {
        if (length(offset) == 1)
            offset <- rep(offset, NROW(Y))
        else if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    par.names <- colnames(X)
    fit <- fit.proc(x = X, y = Y, weights = weights, start = start,
        etastart = etastart, mustart = mustart, offset = offset,
        family = family, control = control.glm, intercept = attr(mt,
            "intercept") > 0)
    if (length(offset) && attr(mt, "intercept") > 0) {
        fit$null.deviance <- glm.fit(x = X[, "(Intercept)", drop = FALSE],
            y = Y, weights = weights, offset = offset, family = family,
            control = control.glm, intercept = TRUE)$deviance
    }
    if (model)
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    ## Move back to the original scale
    if (nPars <- ncol(X)) {
        redundant <- if (br)
            fit$redundant
        else rep.int(0, nPars)
        fit$coefficients <- fit$coefficients/Xmax[!redundant]
        #fit$qr <- qr(sqrt(fit$weights) * Xor[, !redundant])
        fit$qr <- qr(sqrt(fit$weights) * Xor)
        if (br) {
            fit$FisherInfo <- fit$FisherInfo * tcrossprod(Xmax[!redundant])
            fit$control.brglm <- control.brglm
        }
        ####################
        ## Aliasing
        coefs <- rep(NA, ncol(X))
        names(coefs) <- par.names
        coefs[!redundant] <- fit$coefficients
        fit$coefficients <- coefs
        ####################
    }
    fit$control.glm <- control.glm
    if (x)
        fit$x <- Xor
    if (!y)
        fit$y <- NULL
    if (br)
        fit$penalized.deviance <- if (all(family$link == "logit") |
            pl)
            fit$deviance - log(det(fit$FisherInfo))
        else NULL
    fit <- c(fit, list(call = call, formula = formula, terms = mt,
        data = data, offset = offset, method = method, pl = pl,
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt,
            mf)))
    class(fit) <- c("brglm0", "glm", "lm")
    fit
}

#' @export
brglm0.fit <- function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
    mustart = NULL, offset = rep(0, nobs), family = binomial(),
    control = glm.control(), control.brglm = brglm0.control(),
    intercept = TRUE, pl = FALSE)
{
    x <- as.matrix(x)
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    linkfun <- family$linkfun
    dmu.deta <- family$mu.eta
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object")
    if (EMPTY) {
        return(glm.fit(x = x, y = y, weights = weights, start = start,
            etastart = etastart, mustart = mustart, offset = offset,
            family = family, control = control, intercept = intercept))
    }
    valideta <- family$valideta
    if (is.null(valideta))
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu))
        validmu <- function(mu) TRUE
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    ## Suggestion by Kurt Hornik to reset the "warn" option value to
    ## what the user has set prior to the execution of brglm.fit
    warnValue <- options(warn = -1)
    cur.repr <- modifications(family, pl = pl)
    ### Find rows with zero weight
    nonzero.w <- which(weights!=0)
    y.count <- y * weights
    if (is.null(control.brglm$br.consts)) {
      wt <- weights + nvars/nobs
      y.adj <- (y.count + 0.5*nvars/nobs)/wt
    }
    else {
       wt <- weights + 2*control.brglm$br.consts
       y.adj <- (y.count + control.brglm$br.consts)/wt
    }
    # Find any aliased out parameters after the removal of any zero weight observations
    temp.fit <- glm.fit(x = x[nonzero.w,], y = y.adj[nonzero.w],
                        weights = wt[nonzero.w], start = start,
                        etastart = etastart[nonzero.w],
                        mustart = mustart[nonzero.w],
                        offset = offset[nonzero.w], family = family,
                        control = control, intercept = intercept)
    redundant <- is.na(temp.fit$coefficients)
    # Remove columns corresponding to aliased out parameters
    if (any(redundant)) {
      x <- x[, -which(redundant), drop = FALSE]
      nvars <- nvars - sum(redundant)
    }
    # Refit to match the dimension of the original data set
    temp.fit <- glm.fit(x = x, y = y.adj, weights = wt, start = start,
                        etastart = etastart, mustart = mustart,
                        offset = offset, family = family,
                        control = control, intercept = intercept)
    nIter <- 0
    test <- TRUE
    x.t <- t(x)
    while (test & (nIter < control.brglm$br.maxit)) {
        nIter <- nIter + 1
        ps <- temp.fit$fitted.values
        etas <- linkfun(ps)
        ww <- rep(0, nobs)
        ww[nonzero.w] <- temp.fit$weights[nonzero.w]/wt[nonzero.w] * weights[nonzero.w]
        W.X <- sqrt(ww[nonzero.w]) * x[nonzero.w, ]
        XWXinv <- chol2inv(chol(crossprod(W.X)))
        hats <- gethats(nobs, nvars, x.t, XWXinv, ww)
        #hats <- diag(x%*%XWXinv%*%t(ww * x))
        cur.model <- cur.repr(ps)
        wt <- weights + hats * cur.model$at
        y.adj <- rep(0, nobs)
        y.adj[nonzero.w] <- (y.count[nonzero.w] + hats[nonzero.w] * cur.model$ar[nonzero.w])/wt[nonzero.w]
        temp.fit <- glm.fit(x = x, y = y.adj, weights = wt, etastart = etas,
            offset = offset, family = family, control = control,
            intercept = intercept)
        modscore <- t(dmu.deta(etas)/variance(ps) * x) %*% ((y.adj -
            ps) * wt)
        if (control.brglm$br.trace) {
            cat("Iteration:", nIter, "\n")
            cat("Modified scores:", modscore, "\n")
        }
        test <- sum(abs(modscore)) > control.brglm$br.epsilon
    }
    options(warnValue)
    temp.fit$converged <- nIter < control.brglm$br.maxit
    if (!temp.fit$converged)
        warning("Iteration limit reached")
    temp.fit$ModifiedScores <- c(modscore)
    ww <- rep(0, nobs)
    ww[nonzero.w] <- temp.fit$weights[nonzero.w]/wt[nonzero.w] * weights[nonzero.w]
    temp.fit$weights <- ww
    W.X <- sqrt(ww[nonzero.w]) * x[nonzero.w, ]
    temp.fit$FisherInfo <- crossprod(W.X)
    XWXinv <- chol2inv(chol(temp.fit$FisherInfo))
    temp.fit$hats <- gethats(nobs, nvars, x.t, XWXinv, ww)
    temp.fit$qr <- qr(W.X)
    temp.fit$nIter <- nIter
    temp.fit$prior.weights <- weights
    temp.fit$y <- y
    temp.fit$deviance <- sum(dev.resids(temp.fit$y, temp.fit$fitted.values,
        temp.fit$prior.weights))
    temp.fit$cur.model <- cur.model
    temp.fit$redundant <- redundant
    aic <- family$aic
    aic.model <- aic(y, n, ps, weights, temp.fit$deviance)
    temp.fit$aic <- aic.model + 2 * temp.fit$rank
    temp.fit
}

#' @export
print.brglm0 <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    if (x$method == "glm.fit" | !(nPars <- length(coef(x)))) {
      class(x) <- class(x)[-match("brglm", class(x))]
      return(print(x, digits, ...))
    }
    cat("\nCall: ", deparse(x$call), "\n\n")
    if (nPars) {
        cat("Coefficients")
        if (is.character(co <- x$contrasts))
            cat("  [contrasts: ", apply(cbind(names(co), co),
                1, paste, collapse = "="), "]")
        cat(":\n")
        print.default(format(x$coefficients, digits = digits),
            print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n\n")
    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
        x$df.residual, "Residual\n")
    if (nchar(mess <- naprint(x$na.action)))
        cat("  (", mess, ")\n", sep = "")
    if (!is.null(x$penalized.deviance))
        cat("Deviance:\t   ", format(round(x$deviance, digits)),
            "\nPenalized Deviance:", format(round(x$penalized.deviance,
                digits)), "\tAIC:", format(round(x$aic, digits)),
            "\n")
    else cat("Deviance:\t   ", format(round(x$deviance, digits)),
        "\tAIC:", format(round(x$aic, digits)), "\n")
    invisible(x)
}

#' @export
print.summary.brglm0 <-
function (x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor,
    signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    if (length(x$aliased) == 0) {
        cat("\nNo Coefficients\n")
    }
    else {
        df <- if ("df" %in% names(x))
            x[["df"]]
        else NULL
        if (!is.null(df) && (nsingular <- df[3] - df[1]))
            cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
                sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients
        if (!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn,
                colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
            na.print = "NA", ...)
    }
    cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ",
        format(x$dispersion), ")\n\n", apply(cbind(paste(format(c("Null",
            "Residual"), justify = "right"), "deviance:"), format(unlist(x[c("null.deviance",
            "deviance")]), digits = max(5, digits + 1)), " on",
            format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"),
            1, paste, collapse = " "), sep = "")
    if (!is.null(x$penalized.deviance))
        cat("Penalized deviance:", format(round(x$penalized.deviance,
            digits = max(5, digits + 1))), "\n")
    if (nchar(mess <- naprint(x$na.action)))
        cat("  (", mess, ")\n", sep = "")
    cat("AIC: ", format(x$aic, digits = max(4, digits + 1)),
        "\n")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2,
                  digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    cat("\n")
    invisible(x)
}

#' @export
summary.brglm0 <-
function (object, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE,
    ...)
{
    if (object$method == "glm.fit")
        return(summary.glm(object, dispersion = NULL, correlation = FALSE,
            symbolic.cor = FALSE, ...))
    df.r <- object$df.residual
    if (is.null(dispersion))
        dispersion <- 1
    aliased <- is.na(coef(object))
    p <- object$rank
    if (p > 0) {
        p1 <- 1:p
        Qr <- object$qr
        coef.p <- object$coefficients[Qr$pivot[p1]]
        covmat.unscaled <- chol2inv(chol(object$FisherInfo))
        covmat <- dispersion * covmat.unscaled
        var.cf <- diag(covmat)
        s.err <- sqrt(var.cf)
        tvalue <- coef.p/s.err
        dn <- c("Estimate", "Std. Error")
        pvalue <- 2 * pnorm(-abs(tvalue))
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p), c(dn, "z value",
            "Pr(>|z|)"))
        df.f <- NCOL(Qr$qr)
    }
    else {
        coef.table <- matrix(, 0, 4)
        dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
            "t value", "Pr(>|t|)"))
        covmat.unscaled <- covmat <- matrix(, 0, 0)
        df.f <- length(aliased)
    }
    keep <- match(c("call", "terms", "family", "deviance", "aic",
        "contrasts", "df.residual", "null.deviance", "df.null",
        "iter", "na.action", "penalized.deviance"), names(object),
        0)
    ans <- c(object[keep], list(deviance.resid = residuals(object,
        type = "deviance"), coefficients = coef.table, aliased = aliased,
        dispersion = dispersion, df = c(object$rank, df.r, df.f),
        cov.unscaled = covmat.unscaled, cov.scaled = covmat))
    if (correlation && p > 0) {
        dd <- sqrt(diag(covmat.unscaled))
        ans$correlation <- covmat.unscaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.brglm0"
    return(ans)
}

brglm0.control <- function (br.epsilon = 1e-08, br.maxit = 100, br.trace = FALSE,
                                  br.consts = NULL, ...)
{
    if (!is.numeric(br.epsilon) || br.epsilon <= 0)
        stop("value of 'epsilon' must be > 0")
    if (!is.numeric(br.maxit) || br.maxit <= 0)
        stop("maximum number of iterations must be > 0")
    list(br.epsilon = br.epsilon, br.maxit = br.maxit, br.trace = br.trace,
         br.consts = br.consts)
}

gethats <- function (nobs, nvars, x.t, XWXinv, ww)
{
    .C("hatsc", n = as.integer(nobs), p = as.integer(nvars),
        x = as.double(x.t), invfisherinf = as.double(XWXinv),
        w = as.double(ww), hat = double(nobs), PACKAGE = "brglm2")$hat
}

## 'glm.control1' is a minor modification of 'glm.conrol'
## The only different is the addition of a ... argument
## Ioannis Kosmidis <I.Kosmidis@warwick.ac.uk> [15/02/2008]
glm.control1 <- function (epsilon = 1e-08, maxit = 25, trace = FALSE, ...)
{
    glm.control(epsilon, maxit, trace)
}

checkModifications <-
function (fun, Length = 100)
{
    p <- seq(.Machine$double.neg.eps, 1 - 1e-10, length = Length)
    te <- fun(p)
    if (!is.list(te))
        stop("The result should be a list of length two.")
    if (length(te) != 2)
        stop("The result should be a list of length two.")
    if (any(is.na(match(names(te), c("ar", "at")))))
        stop("The result should be a list with elements 'ar' and'at'.")
    if (length(te$ar) != Length)
        stop("'ar' should be of the same length as 'p'")
    if (length(te$at) != Length)
        stop("'at' should be of the same length as 'p'")
    if (any(te$ar >= te$at))
        stop("'ar' cannot take larger values than 'at'")
    if (any(te$ar < 0))
        stop("'ar' cannot be negative")
    if (any(te$ar < 0))
        stop("'at' cannot be negative")
    plot(p, te$at, ylim = c(0, 10), type = "l")
    points(p, te$ar, type = "l", col = "grey")
    drop(TRUE)
}

modifications <-
function (family, pl = FALSE)
{
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    distr.link <- paste(family$family, family$link[1], sep = ".")
    distr.link <- gsub(pattern = "[(]", x = distr.link, replacement = ".")
    distr.link <- gsub(pattern = "[)]", x = distr.link, replacement = ".")
    if (pl) {
        out <- switch(distr.link, binomial.logit = function(p) {
            etas <- family$linkfun(p)
            list(ar = 0.5 * p/p, at = 1 * p/p)
        }, binomial.probit = function(p) {
            etas <- family$linkfun(p)
            list(ar = p * (1 - etas * (etas < 0)/dnorm(etas)),
                at = etas * ((etas >= 0) * (1 - p) - (etas <
                  0) * p)/dnorm(etas) + 0.5/p)
        }, binomial.cloglog = function(p) {
            etas <- family$linkfun(p)
            list(ar = -p/log(1 - p), at = 0.5/p)
        }, binomial.cauchit = function(p) {
            etas <- family$linkfun(p)
            list(ar = -2 * pi * etas * p * (etas < 0) + (p -
                0.5) * (etas >= 0) + p, at = 2 * pi * etas *
                ((etas >= 0) - p) - (p - 0.5)/p * (etas < 0) +
                1)
        }, NULL)
        if (is.null(out))
            out <- match.fun("mpl.custom.family")
    }
    else {
        out <- switch(distr.link, binomial.logit = function(p) {
            etas <- family$linkfun(p)
            list(ar = 0.5 * p/p, at = 1 * p/p)
        }, binomial.probit = function(p) {
            etas <- family$linkfun(p)
            list(ar = -0.5 * p * etas * (etas < 0)/dnorm(etas) +
                p, at = 0.5 * etas * ((etas >= 0) - p)/dnorm(etas) +
                1)
        }, binomial.cloglog = function(p) {
            etas <- family$linkfun(p)
            list(ar = -0.5 * p/log(1 - p) + p, at = 0.5 * p/p +
                1)
        }, binomial.cauchit = function(p) {
            etas <- family$linkfun(p)
            list(ar = -pi * etas * p * (etas < 0) + p, at = pi *
                etas * ((etas >= 0) - p) + 1)
        }, NULL)
        if (is.null(out))
            out <- match.fun("br.custom.family")
    }
    out
}
