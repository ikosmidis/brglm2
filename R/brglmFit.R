#' Fitting function for \code{\link{glm}} for reduced-bias
#' estimation and inference
#'
#' \code{\link{brglmFit}} is a fitting function for \code{\link{glm}}
#' that fits generalized linear models using implicit and explicit
#' bias reduction methods. Currently supported methods include the
#' implicit adjusted scores approach in Firth (1993) and Kosmidis \&
#' Firth (2009), the correction of the asymptotic bias in Cordeiro &
#' McCullagh (1991), and maximum likelihood.  Estimation is performed
#' using a quasi-Fisher scoring iteration based on the iterative
#' correction of the asymptotic bias of the Fisher scoring iterates.
#'
#' @inheritParams stats::glm.fit
#' @param x \code{x} is a design matrix of dimension \code{n * p},
#' @param y \code{y} is a vector of observations of length \code{n}
#' @param control a list of parameters controlling the fitting
#'     process. See \code{\link{brglmControl}} for details.
#' @param start starting values for the parameters in the linear
#'     predictor. If \code{NULL} (default) then the maximum likelihood
#'     estimates are caluclated and used as starting values
#' @param mustart applied only when start is not \code{NULL}. Starting
#'     values for the vector of means to be passed to
#'     \code{\link{glm.fit}} when computing starting values using
#'     maximum likelihood.
#' @param etastart applied only when start is not
#'     \code{NULL}. Starting values for the linear predictor to be
#'     passed to \code{\link{glm.fit}} when computing starting values
#'     using maximum likelihood.
#' @param fixedTotals Effective only when \code{family} is
#'     \code{poisson}. Either \code{NULL} (no effect) or a vector that
#'     indicates which counts must be treated as a group. See
#'     Details for more information and \code{\link{brmultinom}}.
#' @param ... arguments to be used to form the default 'control'
#'     argument if it is not supplied directly.
#'
#' @details
#'
#' Implicit and explicit bias reduction methods are described in
#' detail in Kosmidis (2014). The quasi (or modified) Fisher scoring
#' iteration is described in Kosmidis (2010) and is based on the
#' iterative correction of the asymptotic bias of the Fisher scoring
#' iterates. A quick description of the quasi Fisher scoring iteration
#' is also given in one of the vignettes of the *enrichwith* R package
#' (see,
#' \url{https://cran.r-project.org/web/packages/enrichwith/vignettes/bias.html}).
#'
#'
#' The null deviance is evaluated based on the fitted values using the
#'     bias reduction method specified by the \code{type} argument
#'     (see \code{\link{brglmControl}}).
#'
#' The description of \code{method} argument nd the \code{Fitting
#' functions} section in \code{\link{glm}} gives information on
#' supplying fitting methods to \code{\link{glm}}.
#'
#' If \code{type == "correction"} (see \code{\link{brglmControl}}),
#' then \code{coefficients} and \code{transDispersion} carry the
#' estimated biases as attributes.
#'
#'
#'
#' @seealso \code{\link{glm.fit}} and \code{\link{glm}}
#'
#' @references
#'
#' Cordeiro G. M. & McCullagh, P. (1991). Bias correction in generalized
#' linear models. *Journal of the Royal Statistical Society. Series B
#' (Methodological)*, **53**, 629-643
#'
#' Firth D. (1993). Bias reduction of maximum likelihood estimates,
#' Biometrika, **80**, 27-38
#'
#' Kosmidis I and Firth D (2009). Bias reduction in exponential family
#' nonlinear models. *Biometrika*, **96**, 793-804
#'
#' Kosmidis I and Firth D (2010). A generic algorithm for reducing
#' bias in parametric estimation. *Electronic Journal of Statistics*,
#' **4**, 1097-1112
#'
#' Kosmidis I (2014). Bias in parametric estimation: reduction and
#' useful side-effects. *WIRE Computational Statistics*, **6**,
#' 185-196
#'
#'
#'
#' @examples
#' ## The lizards example from ?brglm::brglm
#' data("lizards")
#' # Fit the model using maximum likelihood
#' lizardsML <- glm(cbind(grahami, opalinus) ~ height + diameter +
#'                      light + time, family = binomial(logit), data = lizards,
#'                  method = "glm.fit")
#' # Now the bias-reduced fit:
#' lizardsBR <- glm(cbind(grahami, opalinus) ~ height + diameter +
#'                      light + time, family = binomial(logit), data = lizards,
#'                  method = "brglmFit")
#' summary(lizardsML)
#' summary(lizardsBR)
#'
#'
#' ## Another example from
#' ## King, Gary, James E. Alt, Nancy Elizabeth Burns and Michael Laver
#' ## (1990).  "A Unified Model of Cabinet Dissolution in Parliamentary
#' ## Democracies", _American Journal of Political Science_, vol. 34,
#' ## no. 3, pp. 846-870.
#'
#' \dontrun{
#' data("coalition", package = "Zelig")
#' # The maximum likelihood fit with log link
#' coalitionML <- glm(duration ~ fract + numst2, family = Gamma, data = coalition)
#' # The bias-reduced fit
#' coalitionBR <- update(coalitionML, method = "brglmFit")
#' # The bias-corrected fit
#' coalitionBC <- update(coalitionML, method = "brglmFit", type = "correction")
#' }
#'
#' \dontrun{
#' ## An example with offsets from Venables & Ripley (2002, p.189)
#' data("anorexia", package = "MASS")
#'
#' anorexML <- glm(Postwt ~ Prewt + Treat + offset(Prewt),
#'                 family = gaussian, data = anorexia)
#' anorexBR <- update(anorexML, method = "brglmFit")
#' anorexBC <- update(anorexML, method = "brglmFit", type = "correction")
#'
#' ## The outputs are identical, because the maximum likelihood
#' ## estimators of the regression parameters are unbiased when family
#' ## is Gaussian, and the bias-reduced estimator of the dispersion is
#' ## the unbiased, by degree of freedom adjustment, estimator of the
#' ## residual variance.
#' summary(anorexML)
#' summary(anorexBR)
#' summary(anorexBC, start = coef(anorexML))
#' }
#'
#' ## endometrial data from Heinze \& Schemper (2002) (see ?endometrial)
#' data("endometrial", package = "brglm2")
#' endometrialML <- glm(HG ~ NV + PI + EH, data = endometrial,
#'                      family = binomial("probit"))
#' endometrialBR <- update(endometrialML, method = "brglmFit",
#'                         type = "adjusted_scores")
#' endometrialBC <- update(endometrialML, method = "brglmFit",
#'                         type = "correction")
#' summary(endometrialML)
#' summary(endometrialBC)
#' summary(endometrialBR)
#'
#' @export
brglmFit <- function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
                      mustart = NULL, offset = rep(0, nobs), family = gaussian(),
                      control = list(), intercept = TRUE,
                      ## Arguments that glm will not use in its call to brglmFit (be wise with defaults!)
                      fixedTotals = NULL)
{

    traceFun <- function(what = "coefficient") {
        if (iter %% control$trace == 0) {
            if (what == "coefficient") {
                st <-  max(abs(stepBeta), na.rm = TRUE)
                gr <- max(abs(adjustedGradBeta), na.rm = TRUE)
                cat("Coefficients update:\t")
                cat("Outer/Inner iteration:\t", sprintf("%03d", iter), "   |", sprintf("%03d", stepFactor), "\n")
            }
            else {
                cat("\n")
                st <- abs(stepZeta)
                gr <- abs(adjustedGradZeta)
                cat("Dispersion update:\t")
                cat("Outer iteration:\t", sprintf("%03d", iter), "\n")
            }
            cat("max |step|:", format(round(st, 6), nsmall = 6, scientific = FALSE), "\t",
                "max |gradient|:", format(round(gr, 6), nsmall = 6, scientific = FALSE), "\n")
            if (what == "coefficient") {

            }
        }
    }

    ## fitFun, grad, info and bias are ALWAYS in beta, dispersion parameterization
    fitFun <- function(pars, y, what = "mean", scaleTotals = FALSE, qr = TRUE) {
        betas <- pars[seq.int(nvars)]
        dispersion <- pars[nvars + 1]
        prec <- 1/dispersion
        etas <- drop(x %*% betas + offset)
        mus <- linkinv(etas)
        if (scaleTotals) {
            ## Rescale mus
            musTotals <-  as.vector(tapply(mus, fixedTotals, sum))[fixedTotals]
            mus <- mus * rowTotals / musTotals
            etas <- linkfun(mus)
        }
        out <- list(precision = prec,
                    betas = betas,
                    dispersion = dispersion,
                    etas = etas,
                    mus = mus,
                    scaleTotals = scaleTotals)
        if (what == "mean") {
                d1mus <- mu.eta(etas)
                d2mus <- d2mu.deta(etas)
                varmus <- variance(mus)
                workingWeights <- weights * d1mus^2 / varmus
                wx <- sqrt(workingWeights) * x
                out$d1mus <- d1mus
                out$d2mus <- d2mus
                out$varmus <- varmus
                out$workingWeights <- workingWeights
                if (qr) out$qrDecomposition <- qr(wx)
        }
        if (what == "dispersion") {
            zetas <- -weights * prec
            out$zetas <- zetas
            ## Evaluate the derivatives of the a function only for
            ## objervations with non-zero weight
            d1afuns <- d2afuns <- d3afuns <- rep(NA, nobs)
            d1afuns[keep] <- d1afun(zetas[keep])
            d2afuns[keep] <- d2afun(zetas[keep])
            d3afuns[keep] <- d3afun(zetas[keep])
            out$d2afuns <- d2afuns
            out$d3afuns <- d3afuns
            out$devianceResiduals <- dev.resids(y, mus, weights)
            out$EdevianceResiduals <- weights * d1afuns
        }
        out
    }

    gradFun <- function(pars, what = "mean", fit = NULL) {
        if (is.null(fit)) {
            fit <- fitFun(pars, y = y, what = what, qr = FALSE)
        }
        with(fit, {
            if (what == "mean") {
                scoreComponents <- weights * d1mus  * (y - mus) / varmus * x
                return(precision * .colSums(scoreComponents, nobs, nvars, TRUE))
                ## return(precision * colSums(scoreComponents))
            }
            if (what == "dispersion") {
                return(1/2 * precision^2 * sum(devianceResiduals - EdevianceResiduals, na.rm = TRUE))
            }
        })
    }

    infoFun <- function(pars, what = "mean", fit = NULL, inverse = FALSE) {
        if (is.null(fit)) {
            fit <- fitFun(pars, y = y, what = what, qr = TRUE)
        }
        with(fit, {
            if (what == "mean") {
                Rmat <- qr.R(qrDecomposition)
                if (inverse) {
                    return(dispersion * tcrossprod(solve(Rmat)))
                }
                else {
                    return(precision * crossprod(Rmat))
                }
            }
            if (what == "dispersion") {
                info <- 0.5 * sum(weights^2 * d2afuns, na.rm = TRUE)/disp^4
                if (inverse) {
                    return(1/info)
                }
                else {
                    return(info)
                }
            }
        })
    }

    hats <- function(pars, fit = NULL) {
        if (is.null(fit)) {
            fit <- fitFun(pars, y = y, what = "mean", qr = TRUE)
        }
        with(fit, {
            Qmat <- qr.Q(qrDecomposition)
            .rowSums(Qmat * Qmat, nobs, nvars, TRUE)
        })
    }

    adjustmentFun <- function(pars, what = "mean", fit = NULL) {
        if (is.null(fit)) {
            fit <- fitFun(pars, y = y, what = what, qr = TRUE)
        }
        with(fit, {
            if (what == "mean") {
                hatvalues <- hats(pars, fit = fit)
                ## User only observations with keep = TRUE to ensure that no division with zero takes place
                return(.colSums(0.5 * hatvalues * d2mus/d1mus * x, nobs, nvars, TRUE))
                ## return(colSums(0.5 * hatvalues * d2mus/d1mus * x))
            }
            if (what == "dispersion") {
                s1 <- sum(weights^3 * d3afuns, na.rm = TRUE)
                s2 <- sum(weights^2 * d2afuns, na.rm = TRUE)
                return((nvars - 2)/(2 * dispersion) + s1/(2 * dispersion^2 * s2))
            }
        })
    }

    ## Estimate the ML of the dispersion parameter for gaussian, gamma and inverse Gaussian
    ## Set the dispersion to 1 if poisson or binomial
    ## coefs is only the regression parameters
    estimateDispersion <- function(coefs, y) {
        if (noDispersion) {
            disp <- 1
            dispML <- 1
        }
        else {
            if (dfResidual > 0) {
                dispFit <- try(uniroot(f = function(phi) {
                    theta <- c(coefs, phi)
                    cfit <- fitFun(theta, y = y, what = "dispersion", qr = FALSE)
                    gradFun(theta, what = "dispersion", fit = cfit)
                }, lower = .Machine$double.eps, upper = 10000, tol = control$epsilon), silent = FALSE)
                if (inherits(dispFit, "try-error")) {
                    warning("the ML estimate of the dispersion could not be calculated. An alternative estimate had been used as starting value.")
                    dispML <- NA
                    disp <- NA
                }
                else {
                    disp <- dispML <- dispFit$root
                }
            }
            else { ## if the model is saturated dispML is NA
                disp <- 1 ## A convenient value
                dispML <- NA
            }
        }
        list(disp = disp, dispML = dispML)
    }

    refit <- function(y, coefs_start = NULL) {
        ## Estimate Beta
        coefs <- coef(glm.fit(x = x, y = y, weights = weights,
                              start = coefs_start,
                              offset = offset,
                              family = family,
                              control = list(epsilon = control$epsilon,
                                             maxit = 10000, trace = FALSE),
                              intercept = intercept))
        dispList <- estimateDispersion(coefs, y)
        disp <- dispList$disp
        dispML <- dispList$dispML
        transdisp <- eval(control$Trans)
        c(coefs, transdisp)
    }

    ## TODO: implement adjustment for IBLA

    control0 <- control
    control <- do.call("brglmControl", control)

    ## Get type
    isML <- control$type == "maximum_likelihood"
    isCor <- control$type == "correction"
    noDispersion <- family$family %in% c("poisson", "binomial")
    justEvaluate <- control$maxit == 0

    ## If fixedTotals is specified the compute rowTotals
    if (is.null(fixedTotals)) {
        hasFixedTotals <- FALSE
    }
    else {
        if (family$family == "poisson") {
            rowTotals <-  as.vector(tapply(y, fixedTotals, sum))[fixedTotals]
            hasFixedTotals <- TRUE
        }
        else {
            hasFixedTotals <- FALSE
        }
    }

    ## Ensure x is a matrix, extract variable names, observation
    ## names, nobs, nvars, and initialize weights and offsets if
    ## needed
    x <- as.matrix(x)
    coefNames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y)) rownames(y) else names(y)
    converged <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) {
        weights <- rep.int(1, nobs)
    }
    if (missingOffset <- is.null(offset)) {
        offset <- rep.int(0, nobs)
    }

    ## Enrich the family object with the required derivatives There is
    ## scope for improvement here by adding an argument to enrich*
    ## functions that controls what you get (e.g. only d2mu.deta and
    ## d1afun-d3afun are needed for bias reduction)
    linkglm <- make.link(family$link)
    family <- enrichwith::enrich(family, with = "function a derivatives")
    linkglm <- enrichwith::enrich(linkglm, with = "d2mu.deta")

    ## Extract functions from the enriched family object
    variance <- family$variance
    linkinv <- linkglm$linkinv
    linkfun <- linkglm$linkfun
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object",
             call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- linkglm$mu.eta
    d2mu.deta <- linkglm$d2mu.deta
    d1afun <- family$d1afun
    d2afun <- family$d2afun
    d3afun <- family$d3afun
    simulate <- family$simulate
    d1TransDisp <- DD(control$Trans, "disp", order = 1)
    d2TransDisp <- DD(control$Trans, "disp", order = 2)

    ## Check for invalid etas and mus
    valideta <- unless_null(family$valideta, function(eta) TRUE)
    validmu <- unless_null(family$validmu, function(mu) TRUE)

    ## FIXME: mustart and etastart set to NULL by default
    mustart <- NULL
    etastart <- NULL

    eval(family$initialize)

    ## If there are no covariates in the model then evaluate only the offset
    if (EMPTY) {
        etas <- rep.int(0, nobs) + offset
        if (!valideta(etas))
            stop("invalid linear predictor values in empty model", call. = FALSE)
        mus <- linkinv(etas)
        if (!validmu(mus))
            stop("invalid fitted means in empty model", call. = FALSE)
        ## deviance <- sum(dev.resids(y, mus, weights))
        workingWeights <- ((weights * mu.eta(etas)^2)/variance(mus))^0.5
        residuals <- (y - mus)/mu.eta(etas)
        keep <- rep(TRUE, length(residuals))
        boundary <- converged <- TRUE
        coefsAll <- numeric()
        rank <- 0
        iter <- 0L
    }
    else {
        boundary <- converged <- FALSE

        ## Detect aliasing
        qrx <- qr(x)
        rank <- qrx$rank
        isFullRank <- all.equal(rank, nvars, tolerance = 1e-06)

        if (!isTRUE(isFullRank)) {
            aliased <- qrx$pivot[seq.int(qrx$rank + 1, nvars)]
            Xall <- x
            x <- x[, -aliased]
            nvarsAll <- nvars
            nvars <- ncol(x)
            coefNamesAll <- coefNames
            coefNames <- coefNames[-aliased]
        }
        else {
            nvarsAll <- nvars
            coefNamesAll <- coefNames
        }

        coefsAll <- structure(rep(NA, nvarsAll), .Names = coefNamesAll)
        keep <- weights > 0

        ## Check for zero weights
        ## if (any(!keep)) {
        ##     warning("Observations with non-positive weights have been omited from the computations")
        ## }
        nkeep <- sum(keep)
        dfResidual <- nkeep - rank

        ## Handle starting values

        ## If start is NULL then start at the ML estimator else use start
        if (is.null(start)) {
            ## Adjust counts if binomial or Poisson in order to avoid infinite estimates
            if (family$family == "binomial") {
                weights.adj <- weights + (!(isCor)) * nvars/nobs
                y.adj <- (weights * y + (!(isCor)) * 0.5 * nvars/nobs)/weights.adj
            }
            else {
                weights.adj <- weights
                y.adj <- y + if (family$family == "poisson") (!(isCor)) * 0.5 * nvars/nobs else 0
            }
            ## ML fit to get starting values
            warn <- getOption("warn")
            ## Get startng values and kill warnings whilst doing that
            options(warn = -1)

            tempFit <- glm.fit(x = x, y = y.adj, weights = weights.adj,
                               etastart = etastart, mustart = mustart,
                               offset = offset, family = family,
                               control = list(epsilon = control$epsilon,
                                              maxit = 10000, trace = FALSE),
                               intercept = intercept)
            ## Set warn to its original value
            options(warn = warn)
            coefs <- coef(tempFit)
            names(coefs) <- coefNames
            dispList <- estimateDispersion(coefs, y = y)
            disp <- dispList$disp
            if (is.na(disp)) disp <- var(y)/variance(sum(weights * y)/sum(weights))
            dispML <- dispList$dispML
            transdisp <- eval(control$Trans)
        }
        else {
            if ((length(start) == nvarsAll) & is.numeric(start)) {
                coefsAll <- start
                names(coefsAll) <- coefNamesAll
                if (!isTRUE(isFullRank)) {
                    coefsAll[aliased] <- NA
                    coefs <- coefsAll[-aliased]
                }
                else {
                    coefs <- coefsAll
                }

                ## Estimate dispersion based on current value for coefs
                dispList <- estimateDispersion(coefs, y = y)
                disp <- dispList$disp
                if (is.na(disp)) disp <- var(y)/variance(sum(weights * y)/sum(weights))
                dispML <- dispList$dispML
                transdisp <- eval(control$Trans)
            }

            if ((length(start) == nvarsAll + 1) & is.numeric(start)) {
                coefsAll <- start[seq.int(nvarsAll)]
                names(coefsAll) <- coefNamesAll
                if (!isTRUE(isFullRank)) {
                    coefsAll[aliased] <- NA
                    coefs <- coefsAll[-aliased]
                }
                else {
                    coefs <- coefsAll
                }
                transdisp <- start[nvarsAll + 1]
                dispML <- NA
                disp <- eval(control$inverseTrans)
            }

            if (length(start) > nvarsAll + 1 | length(start) < nvarsAll) {
                stop(paste(paste(gettextf("length of 'start' should be equal to %d and correspond to initial coefs for %s", nvarsAll, paste(deparse(coefNamesAll), collapse = ", "), "or", gettextf("to %d and also include a starting value for the transformed dispersion", nvarsAll + 1)))), domain = NA)
            }
        }

        adjustedGradAll <- rep(NA, nvarsAll + 1)
        names(adjustedGradAll) <- c(coefNamesAll, "Transformed dispersion")
        if (isCor) {
            control$maxit <- 1
            control$slowit <- 1
        }
        objCur <- .Machine$integer.max


        ## Main iterations

        ## If maxit == 0
        if (justEvaluate) {
            iter <- 0
            theta <- c(coefs, disp)
            ## If fixedTotals is provided (i.e. multinomial
            ## regression via the Poisson trick) then everything
            ## except the score function needs to be evaluated at
            ## the scaled fitted means
            fitBeta <- fitFun(theta, y = y, what = "mean", scaleTotals = hasFixedTotals, qr = TRUE)
            gradBeta <-  gradFun(theta, fit = if (hasFixedTotals) NULL else fitBeta, what = "mean")
                                cInverseInfoBeta <- try(infoFun(theta, inverse = TRUE, fit = fitBeta, what = "mean"))
            if (failedInv <- inherits(cInverseInfoBeta, "try-error")) {
                warning("failed to invert the information matrix: iteration stopped prematurely")
                break
            }
            else {
                inverseInfoBeta <- cInverseInfoBeta
            }
            if (isML) {
                adjustmentBeta <- 0
                failedAdj <- FALSE
            }
            else {
                adjustmentBeta <- adjustmentFun(theta, fit = fitBeta, what = "mean")
                if (failedAdj <- any(is.na(adjustmentBeta))) {
                    warning("failed to calculate the bias-reducing score adjustment: iteration stopped prematurely")
                    break
                }
            }

            adjustedGradBeta <- gradBeta + adjustmentBeta

            ## ADD dispersion + dispersion starting values

            if (noDispersion) {
                    disp <- 1
                    transdisp <- eval(control$Trans)
                    adjustedGradZeta <- NA
                    adjustmentZeta <- NA
                    infoDispersion <- NA
                    d1zeta <- NA
            }
            else {
                if (dfResidual > 0) {
                    theta <- c(coefs, disp)
                    fitDispersion <- fitFun(theta, y = y, what = "dispersion", qr = TRUE)
                    gradDispersion <-  gradFun(theta, fit = fitDispersion, what = "dispersion")
                    infoDispersion <- infoFun(theta, inverse = FALSE, fit = fitDispersion, what = "dispersion")
                    d1zeta <- eval(d1TransDisp)
                    if (isML) {
                        adjustedGradZeta <- 0
                    }
                    else {
                        adjustmentDispersion <- adjustmentFun(theta, fit = fitDispersion, what = "dispersion")
                        ## The adjustment for transDisp (use Kosmidis & Firth, 2010, Remark 3 for derivation)
                        ## FIX: some redundancy below...
                        adjustmentZeta <- adjustmentDispersion/d1zeta - 0.5 * eval(d2TransDisp) / d1zeta^2
                        adjustedGradZeta <- gradDispersion/d1zeta + adjustmentZeta
                    }
                }
                else {
                    disp <- 1 # No effect to the adjusted scores
                    transdisp <- eval(control$Trans)
                    adjustedGradZeta <- NA
                }
            }
        }
        else {
            for (iter in seq.int(control$maxit)) {
                stepFactor <- 0
                testhalf <- TRUE
                coefsPrev <- coefs
                objPrev <- objCur
                ## Bias-reduced estimation of mean effects

                while (testhalf & stepFactor < control$maxStepFactor) {
                    theta <- c(coefs, disp)

                    ## If fixedTotals is provided (i.e. multinomial
                    ## regression via the Poisson trick) then everything
                    ## except the score function needs to be evaluated at
                    ## the scaled fitted means
                    fitBeta <- fitFun(theta, y = y, what = "mean", scaleTotals = hasFixedTotals, qr = TRUE)

                    gradBeta <-  gradFun(theta, fit = if (hasFixedTotals) NULL else fitBeta, what = "mean")

                    cInverseInfoBeta <- try(infoFun(theta, inverse = TRUE, fit = fitBeta, what = "mean"))


                    if (failedInv <- inherits(cInverseInfoBeta, "try-error")) {
                        warning("failed to invert the information matrix: iteration stopped prematurely")
                        break
                    }
                    else {
                        inverseInfoBeta <- cInverseInfoBeta
                    }

                    if (isML) {
                        adjustmentBeta <- 0
                        failedAdj <- FALSE
                    }
                    else {
                        adjustmentBeta <- adjustmentFun(theta, fit = fitBeta, what = "mean")
                        if (failedAdj <- any(is.na(adjustmentBeta))) {
                            warning("failed to calculate the bias-reducing score adjustment: iteration stopped prematurely")
                            break
                        }
                    }

                    adjustedGradBeta <- gradBeta + adjustmentBeta
                    stepBeta <- drop(inverseInfoBeta %*% adjustedGradBeta)
                    coefs <- coefsPrev + 2^(-stepFactor) * stepBeta
                    objCur <- sum(stepBeta^2)
                    stepFactor <- stepFactor + 1
                    testhalf <- objCur > objPrev
                    if (control$trace) {
                        traceFun(what = "coefficient")
                    }
                }

                ## Bias-reduced estimation of (transformed) dispersion
                if (noDispersion) {
                    disp <- 1
                    transdisp <- eval(control$Trans)
                    adjustedGradZeta <- NA
                    adjustmentZeta <- NA
                    infoDispersion <- NA
                    d1zeta <- NA
                    stepZeta <- 0
                }
                else {
                    if (dfResidual > 0) {
                        theta <- c(coefs, disp)
                        fitDispersion <- fitFun(theta, y = y, what = "dispersion", qr = TRUE)
                        gradDispersion <-  gradFun(theta, fit = fitDispersion, what = "dispersion")
                        infoDispersion <- infoFun(theta, inverse = FALSE, fit = fitDispersion, what = "dispersion")
                        d1zeta <- eval(d1TransDisp)
                        if (isML) {
                            adjustedGradZeta <- 0
                        }
                        else {
                            adjustmentDispersion <- adjustmentFun(theta, fit = fitDispersion, what = "dispersion")
                            ## The adjustment for transDisp (use Kosmidis & Firth, 2010, Remark 3 for derivation)
                            ## FIX: some redundancy below...
                            adjustmentZeta <- adjustmentDispersion/d1zeta - 0.5 * eval(d2TransDisp) / d1zeta^2
                            adjustedGradZeta <- gradDispersion/d1zeta + adjustmentZeta
                        }
                        stepZeta <- as.vector(d1zeta^2 * adjustedGradZeta/infoDispersion)
                        transdisp <- transdisp + control$slowit * stepZeta
                        disp <- eval(control$inverseTrans)
                    }
                    else {
                        disp <- 1 # No effect to the adjusted scores
                        transdisp <- eval(control$Trans)
                        adjustedGradZeta <- NA
                        stepZeta <- 0
                    }
                }

                if (control$trace) {
                    traceFun(what = "dispersion")
                }

                if (failedAdj | failedInv | all(abs(c(abs(stepBeta), abs(stepZeta))) < control$epsilon, na.rm = TRUE)) {
                    break
                }
            }
        }

        ## Objects used from above iteration
        ## adjustedGradBeta
        ## adjustedGradZeta
        ## coefs
        ## disp
        ##

        adjustedGradAll[coefNames] <- adjustedGradBeta
        adjustedGradAll["Transformed dispersion"] <- adjustedGradZeta

        coefsAll[coefNames] <- coefs

        ## Convergence analysis
        if ((failedInv | failedAdj | iter >= control$maxit) & !(isCor)) {
            warning("brglmFit: algorithm did not converge", call. = FALSE)
            converged <- FALSE
        }
        else {
            converged <- TRUE
        }

        if (boundary) {
            warning("brglmFit: algorithm stopped at boundary value", call. = FALSE)
        }

        ## QR decomposition and fitted values are AT the final value
        ## for the coefficients

        ## QR decomposition for cov.unscaled
        if (!isTRUE(isFullRank)) {
            x <- Xall
            coefs <- coefsAll
            coefs[is.na(coefs)] <- 0
            nvars <- nvarsAll
        }

        ## If hasFixedTotals = TRUE, then scale fitted values before
        ## calculating QR decompositions, fitted values, etas,
        ## residuals and workingWeights
        fitBeta <- fitFun(c(coefs, disp), y = y, what = "mean", scaleTotals = hasFixedTotals, qr = TRUE)
        qr.Wx <- fitBeta$qrDecomposition

        mus <- fitBeta$mus
        etas <- fitBeta$etas
        ## Residuals
        residuals <- with(fitBeta, (y - mus)/d1mus)
        workingWeights <- fitBeta$workingWeights

        ## Fisher information for the transformed dispersion
        d1zeta <- eval(d1TransDisp)
        if (!noDispersion) {
            fitDispersion <- fitFun(c(coefs, disp), y = y, what = "dispersion", qr = TRUE)
            infoTransDisp <- infoFun(c(coefs, disp), inverse = FALSE, fit = fitDispersion, what = "dispersion")/d1zeta^2
        }

        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if (any(mus > 1 - eps) || any(mus < eps)) {
                warning("brglmFit: fitted probabilities numerically 0 or 1 occurred", call. = FALSE)
                boundary <- TRUE
            }
        }
        if (family$family == "poisson") {
            if (any(mus < eps)) {
                warning("brglmFit: fitted rates numerically 0 occurred", call. = FALSE)
                boundary <- TRUE
            }
        }
        if (dfResidual == 0) disp <- NaN

        ## Estimate of first-order bias FROM the last iteration (so
        ## not at the final value for the coefficients)
        if (isML) {
            ## For now... To be set to calculate biases at a later version
            biasBeta <- biasZeta <- NULL
        }
        else {
            biasBeta <- -drop(inverseInfoBeta %*% adjustmentBeta)
            biasZeta <- -adjustmentZeta/infoDispersion/d1zeta^2
            biasesCoefsAll <- coefsAll
            biasesCoefsAll[coefNames] <- biasBeta
            ## If correction has been requested then add estimated biases an attribute to the coefficients
            if (isCor) {
                attr(coefsAll, "biases") <- biasesCoefsAll
                attr(transdisp, "biases") <- biasZeta
            }
        }


    }

    ## Working weights
    wt <- rep.int(0, nobs)
    wt[keep] <- workingWeights[keep]
    names(wt) <- names(residuals) <- names(mus) <- names(etas) <- names(weights) <- names(y) <- ynames
    ## For the null deviance:
    ##
    ##
    ## If there is an intercept but not an offset then the ML fitted
    ## value is the weighted average and is calculated easily below if
    ## ML is used
    ##

    control0$epsilon <- control$epsilon
    control0$maxit <- control$maxit
    control0$type <- control$type
    control0$slowit <- control$slowit
    if (intercept & missingOffset) {
        nullFit <- brglmFit(x = x[, "(Intercept)"], y = y, weights = weights,
                            offset = rep(0, nobs), family = family, intercept = TRUE,
                            control = control0[c("epsilon", "maxit", "type", "dispTrans", "slowit")])
        nullmus <- nullFit$fitted
    }
    ## If there is an offset but not an intercept then the fitted
    ## value is the inverse link evaluated at the offset
    ##
    ## If there is neither an offset nor an intercept then the fitted
    ## values is the inverse link at zero (and hence covered by
    ## linkinv(offset) because offset is zero
    if (!intercept) {
        nullmus <- linkinv(offset)
    }
    ## If there is an intercept and an offset then, for calculating
    ## the null deviance glm will make a call to the fitter to fit the
    ## glm with intercept and the offset
    if (intercept & !missingOffset) {
        nullmus <- mus
        ## doen't really matter what nullmus is set to. glm will make
        ## a new call to brglmFit and use the deviance from that call
        ## as null
    }
    nulldev <- sum(dev.resids(y, nullmus, weights))
    nulldf <- nkeep - as.integer(intercept)
    deviance <- sum(dev.resids(y, mus, weights))
    aic.model <- aic(y, n, mus, weights, deviance) + 2 * rank

    list(coefficients = coefsAll,
         residuals = residuals,
         fitted.values = mus,
         ## TODO: see effects?
         ## effects = if (!EMPTY) effects,
         R = if (!EMPTY) qr.R(qr.Wx),
         rank = rank,
         qr = if (!EMPTY) structure(qr.Wx[c("qr", "rank", "qraux", "pivot", "tol")], class = "qr"),
         family = family,
         linear.predictors = etas,
         deviance = deviance,
         aic = aic.model,
         null.deviance = nulldev,
         iter = iter,
         weights = wt,
         prior.weights = weights,
         df.residual = dfResidual,
         df.null = nulldf,
         y = y,
         converged = converged,
         boundary = boundary,
         dispersion = disp,
         dispersionML = dispML,
         transDispersion = transdisp,
         infoTransDispersion = if (noDispersion) NA else infoTransDisp,
         grad =  adjustedGradAll,
         dispTrans = control$dispTrans,
         ## cov.unscaled = tcrossprod(Rmat),
         type = control$type,
         class = "brglmFit")
}

#' @export
coef.brglmFit <- function(object, model = c("mean", "full", "dispersion"), ...) {
    model <- match.arg(model)
    switch(model,
           mean = {
               object$coefficients
           },
           dispersion = {
               transDisp <- object$transDispersion
               names(transDisp) <- paste0(object$dispTrans, "(dispersion)")
               transDisp
               ## This will ALWAYS be on the scale of the TRANSFORMED dispersion
           },
           full = {
               transDisp <- object$transDispersion
               ntd <- paste0(object$dispTrans, "(dispersion)")
               names(transDisp) <- ntd
               coefs <- object$coefficients
               thetaTrans <- c(coefs, transDisp)
               if (object$type == "correction") {
                   bcf <- attr(coefs, "biases")
                   btd <- attr(transDisp, "biases")
                   names(btd) <- ntd
                   attr(thetaTrans, "biases") <- c(bcf, btd)
               }
               thetaTrans
           })
}

#' \code{summary} method for \code{\link{brglmFit}} objects
#'
#' @inheritParams stats::summary.glm
#'
#' @details The interface of the summary method for
#'     \code{\link{brglmFit}} objects is identical to that of
#'     \code{\link{glm}} objects. The summary method for
#'     \code{\link{brglmFit}} objects computes the p-values of the
#'     individual Wald statistics based on the standard normal
#'     distribution, unless the family is Gaussian, in which case a t
#'     distribution with appropriate degrees of freedom is used.
#'
#' @seealso \code{\link{summary.glm}} and \code{\link{glm}}
#'
#' @examples
#' ## For examples see examples(brglmFit)
#'
#' @method summary brglmFit
#' @export
summary.brglmFit <- function(object, dispersion = NULL,
                             correlation = FALSE, symbolic.cor = FALSE,
                             ...) {
    if (is.null(dispersion)) {
        if (object$family$family == "Gaussian") {
            dispersion <- NULL
        }
        else {
            dispersion <- object$dispersion
        }
    }
    summary.glm(object, dispersion = dispersion,
                correlation = correlation,
                symbolic.cor = symbolic.cor, ...)
}

#' Method for computing confidence intervals for one or more
#' regression parameters in a \code{\link{brglmFit}} object
#'
#' @inheritParams stats::confint
#'
#' @method confint brglmFit
#' @export
confint.brglmFit <- function(object, parm, level = 0.95, ...) {
    confint.default(object, parm, level, ...)
}

#' Return the variance-covariance matrix for the regression parameters
#' in a \code{\link{brglmFit}} object
#'
#' @inheritParams stats::vcov
#' @param model character specyfying for which component of the model coefficients shoould be extracted
#'
#' @method vcov brglmFit
#' @export
vcov.brglmFit <- function(object, model = c("mean", "full", "dispersion"), ...) {
    model <- match.arg(model)
    switch(model,
           mean = {
               summary.brglmFit(object, ...)$cov.scaled
           },
           dispersion = {
               vtd <- 1/object$infoTransDispersion
               ntd <- paste0(object$dispTrans, "(dispersion)")
               names(vtd) <- ntd
               vtd
           },
           full = {
               vcoefs <- summary.brglmFit(object, ...)$cov.scaled
               vtd <- 1/object$infoTransDispersion
               nCoefsAll <- c(rownames(vcoefs), paste0(object$dispTrans, "(dispersion)"))
               vCoefsAll <- cbind(rbind(vcoefs, 0),
                                  c(numeric(nrow(vcoefs)), vtd))
               dimnames(vCoefsAll) <- list(nCoefsAll, nCoefsAll)
               vCoefsAll
           })
}


DD <- function(expr,name, order = 1) {
    if(order < 1) stop("'order' must be >= 1")
    if(order == 1) D(expr,name)
    else DD(D(expr, name), name, order - 1)
}




