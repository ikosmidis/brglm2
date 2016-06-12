brglmFitStable <- function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
                      mustart = NULL, offset = rep(0, nobs), family = gaussian(),
                      control = list(), intercept = TRUE)
{
    control <- do.call("brglmControl", control)

    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y))
                  rownames(y)
              else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    ## Enrich the family object with the required derivatives There is
    ## scope for improvement here by adding an argument to enrich*
    ## functions that controls what you get (e.g. only d2mu.deta and
    ## d1afun-d3afun are needed for bias reduction)
    family <- enrichFamily(family)

    ## Extract functions from the enriched family object
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object",
             call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    d2mu.deta <- family$d2mu.deta
    d1afun <- family$d1afun
    d2afun <- family$d2afun
    d3afun <- family$d3afun
    d1TransDisp <- DD(control$Trans, "disp", order = 1)
    d2TransDisp <- DD(control$Trans, "disp", order = 2)

    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta))
            stop("invalid linear predictor values in empty model",
                 call. = FALSE)
        mu <- linkinv(eta)
        if (!validmu(mu))
            stop("invalid fitted means in empty model", call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric()
        iter <- 0L
    }
    else {
        ## Get Starting values
        coefold <- NULL
        eta <- if (!is.null(etastart)) {
                   etastart
               }
               else {
                   if (!is.null(start)) {
                       if (length(start) != nvars) {
                           stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                                         nvars, paste(deparse(xnames), collapse = ", ")),
                                domain = NA)
                       }
                       else {
                           coefold <- start
                           offset + as.vector(if (NCOL(x) == 1L)
                                                  x * start
                                              else x %*% start)
                       }
                   }
                   else {
                       family$linkfun(mustart)
                   }
               }
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta))) {
            stop("cannot find valid starting values: please specify some", call. = FALSE)
        }
        boundary <- conv <- FALSE
        ## Adjust counts if binomial or Poisson in order to avoid infinite estimates
        if (family$family == "binomial") {
            weights.adj <- weights + (!control$correction) * nvars/nobs
            y.adj <- (weights*y + (!control$correction) * 0.5 * nvars/nobs)/weights.adj
        }
        else {
            weights.adj <- weights
            y.adj <- y + if (family$family == "poisson") (!control$correction) * 0.5 * nvars/nobs else 0
        }
        warn <- getOption("warn")
        ## Get startng valuyes and kill warnings whilst doing that
        options(warn = -1)
        temp.fit <- glm.fit(x = x, y = y.adj, weights = weights.adj, start = start,
                            etastart = etastart, mustart = mustart,
                            offset = offset, family = family,
                            control = list(epsilon = control$epsilon,
                                           maxit = 10000, trace = FALSE),
                            intercept = intercept)
        ## Set warn to its original value
        options(warn = warn)
        df.r <- temp.fit$df.residual



        ## Type 1 is regression coefficients, type 2 is dispersion
        fitfun <- function(pars, type = 1) {
        }

        loglik <- function(pars, what = c("mean", "dispersion"), fit = NULL) {

        }

        grad <- function(pars, what = c("mean", "dispersion"), fit = NULL) {

        }

        info <- function(pars, what = c("mean", "dispersion"), fit = NULL) {

        }

        bias <- function(pars, what = c("mean", "dispersion"), fit = NULL) {

        }


        ## Set the dispersion to 1 if poisson or binomial
        if (family$family %in% c("poisson", "binomial")) {
            disp <- 1
            dispML <- 1
            transdisp <- eval(control$Trans)
        }
        ## Estimate the ML of the dispersion parameter for gaussian, gamma and inverse Gaussian
        else {
            if (df.r > 0) {
                ## Remove zero-weight observations from the estimation of the dispersion
                if (any(temp.fit$weights == 0))
                    warning("observations with zero weight not used for estimating dispersion")
                disp <- sum((temp.fit$weights * temp.fit$residuals^2)[temp.fit$weights >
                                                              0])/df.r
                prec <- 1/disp
                mu <- temp.fit$fitted
                devresids <- dev.resids(y, mu, weights)
                for (i in 1L:100) {
                    zetas <- -weights*prec
                    expDevresids <- weights*d1afun(zetas)
                    precScores <- -1/2*sum(devresids - expDevresids)
                    precInfo <- 1/2*sum(weights^2*d2afun(zetas))
                    disp <- disp - disp^2 * precScores/precInfo
                    prec <- 1/disp
                    transdisp <- eval(control$Trans)
                    if (critDispML <- sum(abs(precScores)) < 1e-08) {
                        dispML <- disp
                        break
                    }
                }
                if (!critDispML) {
                    warning("the ML estimate of the dispersion could not be calculated. An alternative estimate had been used as starting value.")
                }
            }
            else { ## if the model is saturated do not dispML is NA
                disp <- 1 ## A random value
                dispML <- NA
                transdisp <- eval(control$Trans)
            }
        }
        nacoefs <- is.na(temp.fit$coefficients)
        coefnames <- names(temp.fit$coefficients)
        X <- x[, coefnames[!nacoefs], drop = FALSE]
        coefs <- temp.fit$coefficients[!nacoefs]
        nvarsRed <- length(coefs)
        adjScores <- rep(NA, nvars + 1)
        names(adjScores) <- c("Transformed Dispersion", coefnames)
        if (control$correction) {
            control$maxit <- 1
            control$slowit <- 1
        }
        objCur <- .Machine$integer.max

        ## Main iterations
        for (iter in 1L:control$maxit) {
            halfstep <- 0
            testhalf <- TRUE
            coefsPrev <- coefs
            objPrev <- objCur

            ## Bias-reduced estimation of mean effects
            while (testhalf & halfstep < 10) {
                ## Do not use observations with zero weight
                good <- weights > 0
                eta <- drop(X %*% coefs + offset)
                mu <- linkinv(eta)
                d1mu <- mu.eta(eta)
                d2mu <- d2mu.deta(eta)
                varmu <- variance(mu)[good]
                if (any(is.na(varmu)))
                    stop("NAs in V(mu)")
                if (any(varmu == 0))
                    stop("0s in V(mu)")
                if (any(is.na(d1mu[good])))
                    stop("NAs in d(mu)/d(eta)")
                good <- (weights > 0) & (d1mu != 0)
                if (all(!good)) {
                    conv <- FALSE
                    warning("no observations informative at iteration ",
                            iter)
                    break
                }
                ## Estimation of coefficients
                w <- rep.int(0, nobs)
                w[good] <- weights[good]*d1mu[good]^2/varmu
                Wx <- sqrt(w[good]) * X[good, ]
                ## XWXinv <- chol2inv(chol(crossprod(Wx)))
                ## hats <- diag(X%*%XWXinv%*%t(w * X))
                qrWx <- qr(Wx)
                Q <- qr.Q(qrWx)
                R <- qr.R(qrWx)
                XWXinv <- solve(crossprod(R))
                nvalid <- sum(good)
                hats <- rowSums(Q * Q)
                adjExpComp <- 0.5 * hats[good] * d2mu[good]/d1mu[good] * X[good, ]
                adjExp <- .Internal(colSums(adjExpComp, nvalid, nvarsRed, FALSE))
                ## adjExp <- try(colSums(adjExpComp))
                scoresComp <- weights[good] * d1mu[good]/varmu * (y - mu)[good] * X[good,]
                scores <- .Internal(colSums(scoresComp, nvalid, nvarsRed, FALSE))
                ## scores <- colSums(scoresComp)
                firstOrderBiasBeta <- - disp * XWXinv %*% adjExp
                adjScoresBeta <- scores/disp + adjExp
                objCur <- sqrt(disp * adjScoresBeta %*% XWXinv %*% adjScoresBeta)
                coefs <- coefsPrev + 2^(-halfstep) * (drop(XWXinv%*%scores) -
                                                      firstOrderBiasBeta)
                halfstep <- halfstep + 1
                testhalf <- objCur > objPrev
            }

            ## Bias-reduced estimation of dispersion parameters
            devresids <- dev.resids(y, mu, weights)
            if (family$family %in% c("poisson", "binomial")) {
                disp <- 1
                transdisp <- eval(control$Trans)
                adjScoresTransDisp <- NA
            }
            else {
                if (df.r > 0) {
                    zetas <- -weights/disp
                    expDevresids <- weights*d1afun(zetas)
                    dispScores <- sum(devresids - expDevresids)/(2*disp^2)
                    dispInfo <- (denom1 <- sum(weights^2*d2afun(zetas)))/(2*disp^4)
                    adjExpDisp <- (nvarsRed - 2)/(2*disp) +
                        sum(weights^3 * d3afun(zetas))/(2*disp^2*denom1)
                    firstOrderBiasDisp <- - adjExpDisp/dispInfo
                    transDispScores <- dispScores/(firstDeriv <- eval(d1TransDisp))
                    transDispInfo <- dispInfo/firstDeriv^2
                    firstOrderBiasTransDisp <- firstDeriv*firstOrderBiasDisp +
                        eval(d2TransDisp)/(2*dispInfo)
                    adjScoresTransDisp <- transDispScores -
                        transDispInfo * firstOrderBiasTransDisp
                    transdisp <- transdisp + control$slowit *
                        (transDispScores/transDispInfo - firstOrderBiasTransDisp)
                    disp <- eval(control$inverseTrans)
                }
                else {
                    disp <- 1 # No effect to the adjusted scores
                    transdisp <- eval(control$Trans)
                    adjScoresTransDisp <- NA
                }
            }
            adjScores[c(TRUE, !nacoefs)] <- c(adjScoresTransDisp, adjScoresBeta)
            if (control$trace) {
                cat("Iteration:", iter, "\n")
                cat("Adjusted scores:", adjScores, "\n")
            }
            if (conv <- sum(abs(adjScores), na.rm = TRUE) < control$epsilon) {
                break
            }
        }

        ## Wrap things up!
        if (df.r == 0) disp <- NaN
        coef <- rep(NA, length(temp.fit$coefficients))
        coef[!nacoefs] <- coefs
        if ((!conv) & (!control$correction)) {
            warning("brglmFit: algorithm did not converge", call. = FALSE)
        }
        if (boundary) {
            warning("brglmFit: algorithm stopped at boundary value", call. = FALSE)
        }
        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps))
                warning("brglmFit: fitted probabilities numerically 0 or 1 occurred",
                        call. = FALSE)
        }
        if (family$family == "poisson") {
            if (any(mu < eps))
                warning("brglmFit: fitted rates numerically 0 occurred",
                        call. = FALSE)
        }
        qr.Wx <- qr(sqrt(w)*x)
        temp.fit$qr <- as.matrix(qr.Wx$qr)
        temp.fit$pivot <- qr.Wx$pivot
        temp.fit$rank <- qr.Wx$rank
        temp.fit$qraux <- qr.Wx$qraux
        xxnames <- xnames[temp.fit$pivot]
        residuals <- (y - mu)/mu.eta(eta)
        nr <- min(nvalid, nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1L:nr, 1L:nvars] <- temp.fit$qr[1L:nr, 1L:nvars]
        }
        else Rmat <- temp.fit$qr[1L:nvars, 1L:nvars]
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat) > col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(temp.fit$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    if (!EMPTY) {
        names(temp.fit$effects) <- c(xxnames[seq_len(temp.fit$rank)], rep.int("", sum(good) - temp.fit$rank))
    }
    ## For the null deviance:
    ## If there is an intercept and an offset then glm will make a call to the fitter to fit the glm with intercept and the offset
    ## If there is an intercept but not an offset then the fitted value is the weighted average and is calculated easily below
    ## If there is an offset but not an intercept then the fitted value is the inverse link evaluated at the offset
    ## If there is neither an offset nor an intercept then the fitted values is the inverse link at zero (and hence covered by linkinv(zero)
    wtdmu <- if (intercept) { ## fitted value
                 sum(weights * y)/sum(weights)
             }
             else {
                 linkinv(offset)
             }
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY)
                0
            else temp.fit$rank
    resdf <- n.ok - rank
    temp.fit$deviance <- sum(devresids)
    aic.model <- aic(y, n, mu, weights, temp.fit$deviance) + 2 * rank
    list(coefficients = coef,
         residuals = residuals,
         fitted.values = mu,
         effects = if (!EMPTY) temp.fit$effects,
         R = if (!EMPTY) Rmat,
         rank = rank,
         qr = if (!EMPTY) structure(temp.fit[c("qr", "rank", "qraux", "pivot", "tol")], class = "qr"),
         family = family,
         linear.predictors = eta,
         deviance = temp.fit$deviance,
         aic = aic.model,
         null.deviance = nulldev,
         iter = iter,
         weights = wt,
         prior.weights = weights,
         df.residual = resdf,
         df.null = nulldf,
         y = y,
         converged = conv,
         boundary = boundary,
         dispersion = disp,
         dispersionML = dispML,
         transDispersion = transdisp,
         biases = firstOrderBiasBeta,
         adjustedScores = adjScores,
         dispTrans = control$dispTrans,
         cov.unscaled = XWXinv,
         class = "brglmFit")
}

summary.brglmFitStable <- function (object, dispersion = object$dispersion,
                              correlation = FALSE, symbolic.cor = FALSE,
                              ...) {
    summary.glm(object, dispersion = dispersion,
                correlation = correlation,
                symbolic.cor = symbolic.cor, ...)
}

customTransStable <- list(Trans = expression(disp),
                    inverseTrans = expression(transdisp))


DDStable <- function(expr,name, order = 1) {
    if(order < 1) stop("'order' must be >= 1")
    if(order == 1) D(expr,name)
    else DD(D(expr, name), name, order - 1)
}

## Suggestion by Kurt Hornik to avoid a warning related to the binding
## of n which is evaluated by family$initialize
if(getRversion() >= "2.15.1") globalVariables("n")

