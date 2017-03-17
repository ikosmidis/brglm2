#' Fitting function for \code{\link{glm}} for reduced-bias
#' estimation and inference
#'
#' \code{\link{brglmFit}} is a fitting function for \code{\link{glm}}
#' that fits generalized linear models using implicit and explicit
#' bias reduction methods. Currently supported methods include the
#' implicit adjusted scores approach in Firth (1993) and Kosmidis \&
#' Firth (2009), the correction of the asymptotic bias in Cordeiro &
#' McCullagh (1991), and maximum likelihood.  Estimation is performed
#' using a quasi Fisher scoring iteration based on the iterative
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
#' @param fixed_totals effective only when \code{family} is
#'     \code{poisson}. Either \code{NULL} (no effect) or a vector that
#'     indicates which counts must be treated as a group. See Details
#'     for more information and \code{\link{brmultinom}}.
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
#' \url{https://cran.r-project.org/package=enrichwith/vignettes/bias.html}).
#'
#'
#' The null deviance is evaluated based on the fitted values using the
#'     method specified by the \code{type} argument (see
#'     \code{\link{brglmControl}}).
#'
#' The description of \code{method} argument and the \code{Fitting
#' functions} section in \code{\link{glm}} gives information on
#' supplying fitting methods to \code{\link{glm}}.
#'
#' If \code{type == "correction"} (see \code{\link{brglmControl}}),
#' then \code{coefficients} and \code{transformed_dispersion} carry the
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
#' ## Democracies", _American Journal of Political Science_, **34**, 846-870
#'
#' \dontrun{
#' data("coalition", package = "brglm2")
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
#'                         type = "AS-mean")
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
                      fixed_totals = NULL)
{

    trace_iteration <- function() {
        if (iter %% control$trace == 0) {
            st <-  max(abs(step_beta), na.rm = TRUE)
            gr <- max(abs(adjusted_grad_beta), na.rm = TRUE)
            cat("Coefficients update:\t")
            cat("Outer/Inner iteration:\t", sprintf("%03d", iter), "   |", sprintf("%03d", step_factor), "\n")
            if (!no_dispersion) {
                st <- abs(step_zeta)
                gr <- abs(adjusted_grad_zeta)
                cat("Dispersion update:\t")
                cat("Outer iteration:\t", sprintf("%03d", iter), "\n")
            }
            cat("max |step|:", format(round(st, 6), nsmall = 6, scientific = FALSE), "\t",
                "max |gradient|:", format(round(gr, 6), nsmall = 6, scientific = FALSE), "\n")
        }
    }

    ## key_quantities, grad, info and bias are ALWAYS in beta, dispersion parameterization
    key_quantities <- function(pars, y, what = "mean", scale_totals = FALSE, qr = TRUE) {
        betas <- pars[seq.int(nvars)]
        dispersion <- pars[nvars + 1]
        prec <- 1/dispersion
        etas <- drop(x %*% betas + offset)
        mus <- linkinv(etas)
        if (scale_totals) {
            ## Rescale mus
            mus_totals <-  as.vector(tapply(mus, fixed_totals, sum))[fixed_totals]
            mus <- mus * row_totals / mus_totals
            etas <- linkfun(mus)
        }
        out <- list(precision = prec,
                    betas = betas,
                    dispersion = dispersion,
                    etas = etas,
                    mus = mus,
                    scale_totals = scale_totals)
        mean_quantities <- function(out) {
            d1mus <- mu.eta(etas)
            d2mus <- d2mu.deta(etas)
            varmus <- variance(mus)
            working_weights <- weights * d1mus^2 / varmus
            wx <- sqrt(working_weights) * x
            out$d1mus <- d1mus
            out$d2mus <- d2mus
            out$varmus <- varmus
            out$working_weights <- working_weights
            if (qr) out$qr_decomposition <- qr(wx)
            out
        }
        dispersion_quantities <- function(out) {
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
            out$deviance_residuals <- dev.resids(y, mus, weights)
            out$Edeviance_residuals <- weights * d1afuns
            out
        }
        if (what == "mean") {
            out <- mean_quantities(out)

        }
        if (what == "dispersion") {
            out <- dispersion_quantities(out)
        }
        if (what == "all") {
            out <- mean_quantities(out)
            out <- dispersion_quantities(out)
        }
        out
    }

    gradient <- function(pars, what = "mean", fit = NULL) {
        if (is.null(fit)) {
            fit <- key_quantities(pars, y = y, what = what, qr = FALSE)
        }
        with(fit, {
            if (what == "mean") {
                score_components <- weights * d1mus  * (y - mus) / varmus * x
                return(precision * .colSums(score_components, nobs, nvars, TRUE))
            }
            if (what == "dispersion") {
                return(1/2 * precision^2 * sum(deviance_residuals - Edeviance_residuals, na.rm = TRUE))
            }
        })
    }

    information <- function(pars, what = "mean", fit = NULL, inverse = FALSE) {
        if (is.null(fit)) {
            fit <- key_quantities(pars, y = y, what = what, qr = TRUE)
        }
        with(fit, {
            if (what == "mean") {
                R_matrix <- qr.R(qr_decomposition)
                if (inverse) {
                    return(dispersion * tcrossprod(solve(R_matrix)))
                }
                else {
                    return(precision * crossprod(R_matrix))
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

    hat_values <- function(pars, fit = NULL) {
        if (is.null(fit)) {
            fit <- key_quantities(pars, y = y, what = "mean", qr = TRUE)
        }
        with(fit, {
            Qmat <- qr.Q(qr_decomposition)
            .rowSums(Qmat * Qmat, nobs, nvars, TRUE)
        })
    }

    refit <- function(y, betas_start = NULL) {
        ## Estimate Beta
        betas <- coef(glm.fit(x = x, y = y, weights = weights,
                              start = betas_start,
                              offset = offset,
                              family = family,
                              control = list(epsilon = control$epsilon,
                                             maxit = 2, trace = FALSE),
                              intercept = intercept))
        betas
    }

    ## Estimate the ML of the dispersion parameter for gaussian, gamma and inverse Gaussian
    ## Set the dispersion to 1 if poisson or binomial
    ## betas is only the regression parameters
    estimate_dispersion <- function(betas, y) {
        if (no_dispersion) {
            disp <- 1
            dispML <- 1
        }
        else {
            if (df_residual > 0) {
                dispFit <- try(uniroot(f = function(phi) {
                    theta <- c(betas, phi)
                    cfit <- key_quantities(theta, y = y, what = "dispersion", qr = FALSE)
                    gradient(theta, what = "dispersion", fit = cfit)
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

    ASmean_adjustment <- function(pars, what = "mean", fit = NULL) {
        if (is.null(fit)) {
            fit <- key_quantities(pars, y = y, what = what, qr = TRUE)
        }
        with(fit, {
            if (what == "mean") {
                hatvalues <- hat_values(pars, fit = fit)
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

    customTransformation <- is.list(control$transformation) & length(control$transformation == 2)
    if (customTransformation) {
        transformation0 <- control$transformation
    }
    control <- do.call("brglmControl", control)

    ## FIXME: Add IBLA
    adjustment_function <- switch(control$type,
                            "correction" = ASmean_adjustment,
                            "AS-mean" = ASmean_adjustment)

    ## Some useful quantities
    is_ML <- control$type == "ML"
    is_correction <- control$type == "correction"
    no_dispersion <- family$family %in% c("poisson", "binomial")
    just_evaluate <- control$maxit == 0
    ## If fixed_totals is specified the compute row_totals
    if (is.null(fixed_totals)) {
        has_fixed_totals <- FALSE
    }
    else {
        if (family$family == "poisson") {
            row_totals <-  as.vector(tapply(y, fixed_totals, sum))[fixed_totals]
            has_fixed_totals <- TRUE
        }
        else {
            has_fixed_totals <- FALSE
        }
    }

    ## Ensure x is a matrix, extract variable names, observation
    ## names, nobs, nvars, and initialize weights and offsets if
    ## needed
    x <- as.matrix(x)
    betas_names <- dimnames(x)[[2L]]
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

    ## Enrich the family object with the required derivatives
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
    d1_transformed_dispersion <- DD(control$Trans, "disp", order = 1)
    d2_transformed_dispersion <- DD(control$Trans, "disp", order = 2)

    ## Check for invalid etas and mus
    valid_eta <- unless_null(family$valideta, function(eta) TRUE)
    valid_mu <- unless_null(family$validmu, function(mu) TRUE)

    ## FIXME: mustart and etastart set to NULL by default
    mustart <- NULL
    etastart <- NULL

    ## Initialize as prescribed in family
    eval(family$initialize)

    ## If there are no covariates in the model then evaluate only the offset
    if (EMPTY) {
        etas <- rep.int(0, nobs) + offset
        if (!valid_eta(etas))
            stop("invalid linear predictor values in empty model", call. = FALSE)
        mus <- linkinv(etas)
        if (!valid_mu(mus))
            stop("invalid fitted means in empty model", call. = FALSE)
        ## deviance <- sum(dev.resids(y, mus, weights))
        working_weights <- ((weights * mu.eta(etas)^2)/variance(mus))^0.5
        residuals <- (y - mus)/mu.eta(etas)
        keep <- rep(TRUE, length(residuals))
        boundary <- converged <- TRUE
        betas_all <- numeric()
        rank <- 0
        iter <- 0L
    }
    else {
        boundary <- converged <- FALSE
        ## Detect aliasing
        qrx <- qr(x)
        rank <- qrx$rank
        is_full_rank <- all.equal(rank, nvars, tolerance = 1e-06)
        if (!isTRUE(is_full_rank)) {
            aliased <- qrx$pivot[seq.int(qrx$rank + 1, nvars)]
            X_all <- x
            x <- x[, -aliased]
            nvars_all <- nvars
            nvars <- ncol(x)
            betas_names_all <- betas_names
            betas_names <- betas_names[-aliased]
        }
        else {
            nvars_all <- nvars
            betas_names_all <- betas_names
        }
        betas_all <- structure(rep(NA, nvars_all), .Names = betas_names_all)
        keep <- weights > 0
        ## Check for zero weights
        ## if (any(!keep)) {
        ##     warning("Observations with non-positive weights have been omited from the computations")
        ## }
        nkeep <- sum(keep)
        df_residual <- nkeep - rank
        ## Handle starting values
        ## If start is NULL then start at the ML estimator else use start
        if (is.null(start)) {
            ## Adjust counts if binomial or Poisson in order to avoid infinite estimates
            if (family$family == "binomial") {
                weights.adj <- weights + (!(is_correction)) * nvars/nobs
                y.adj <- (weights * y + (!(is_correction)) * 0.5 * nvars/nobs)/weights.adj
            }
            else {
                weights.adj <- weights
                y.adj <- y + if (family$family == "poisson") (!(is_correction)) * 0.5 * nvars/nobs else 0
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
            betas <- coef(tempFit)
            names(betas) <- betas_names
            dispList <- estimate_dispersion(betas, y = y)
            disp <- dispList$disp
            if (is.na(disp)) disp <- var(y)/variance(sum(weights * y)/sum(weights))
            dispML <- dispList$dispML
            transformed_dispersion <- eval(control$Trans)
        }
        else {
            if ((length(start) == nvars_all) & is.numeric(start)) {
                betas_all <- start
                names(betas_all) <- betas_names_all
                if (!isTRUE(is_full_rank)) {
                    betas_all[aliased] <- NA
                    betas <- betas_all[-aliased]
                }
                else {
                    betas <- betas_all
                }
                ## Estimate dispersion based on current value for betas
                dispList <- estimate_dispersion(betas, y = y)
                disp <- dispList$disp
                if (is.na(disp)) disp <- var(y)/variance(sum(weights * y)/sum(weights))
                dispML <- dispList$dispML
                transformed_dispersion <- eval(control$Trans)
            }
            if ((length(start) == nvars_all + 1) & is.numeric(start)) {
                betas_all <- start[seq.int(nvars_all)]
                names(betas_all) <- betas_names_all
                if (!isTRUE(is_full_rank)) {
                    betas_all[aliased] <- NA
                    betas <- betas_all[-aliased]
                }
                else {
                    betas <- betas_all
                }
                transformed_dispersion <- start[nvars_all + 1]
                dispML <- NA
                disp <- eval(control$inverseTrans)
            }
            if (length(start) > nvars_all + 1 | length(start) < nvars_all) {
                stop(paste(paste(gettextf("length of 'start' should be equal to %d and correspond to initial betas for %s", nvars_all, paste(deparse(betas_names_all), collapse = ", "), "or", gettextf("to %d and also include a starting value for the transformed dispersion", nvars_all + 1)))), domain = NA)
            }
        }
        adjusted_grad_all <- rep(NA, nvars_all + 1)
        names(adjusted_grad_all) <- c(betas_names_all, "Transformed dispersion")
        if (is_correction) {
            control$maxit <- 1
            control$slowit <- 1
        }
        objCur <- .Machine$integer.max
        ## Main iterations
        ## If maxit == 0
        what <- if (no_dispersion) "mean" else "all"



        if (just_evaluate) {
            iter <- 0
            theta <- c(betas, disp)
            ## If fixed_totals is provided (i.e. multinomial regression
            ## via the Poisson trick) then evaluate everything expect
            ## the score function at the scaled fitted means
            quantities <- key_quantities(theta, y = y, what = what, scale_totals = has_fixed_totals, qr = TRUE)
            grad_beta <-  gradient(theta, fit = if (has_fixed_totals) NULL else quantities, what = "mean")

            inverse_info_beta <- try(information(theta, inverse = TRUE, fit = quantities, what = "mean"))
            if (failed_inversion <- inherits(inverse_info_beta, "try-error")) {
                warning("failed to invert the information matrix: iteration stopped prematurely")
                break
            }
            if (is_ML) {
                adjustment_beta <- 0
                failed_adjustment <- FALSE
            }
            else {
                adjustment_beta <- adjustment_function(theta, fit = quantities, what = "mean")
                if (failed_adjustment <- any(is.na(adjustment_beta))) {
                    warning("failed to calculate the bias-reducing score adjustment: iteration stopped prematurely")
                    break
                }
            }
            adjusted_grad_beta <- grad_beta + adjustment_beta
            if (no_dispersion) {
                    disp <- 1
                    transformed_dispersion <- eval(control$Trans)
                    adjusted_grad_zeta <- NA
                    adjustment_zeta <- NA
                    info_dispersion <- NA
                    d1zeta <- NA
            }
            else {
                if (df_residual > 0) {
                    gradDispersion <-  gradient(theta, fit = quantities, what = "dispersion")
                    info_dispersion <- information(theta, inverse = FALSE, fit = quantities, what = "dispersion")
                    d1zeta <- eval(d1_transformed_dispersion)
                    if (is_ML) {
                        adjusted_grad_zeta <- 0
                    }
                    else {
                        adjustmentDispersion <- adjustment_function(theta, fit = quantities, what = "dispersion")
                        ## The adjustment for transDisp (use Kosmidis & Firth, 2010, Remark 3 for derivation)
                        ## FIXME: some redundancy below...
                        adjustment_zeta <- adjustmentDispersion/d1zeta - 0.5 * eval(d2_transformed_dispersion) / d1zeta^2
                        adjusted_grad_zeta <- gradDispersion/d1zeta + adjustment_zeta
                    }
                }
                else {
                    disp <- 1 # No effect to the adjusted scores
                    transformed_dispersion <- eval(control$Trans)
                    adjusted_grad_zeta <- NA
                }
            }


            ## iter <- 0
            ## theta <- c(betas, disp)
            ## ## If fixed_totals is provided (i.e. multinomial regression
            ## ## via the Poisson trick) then evaluate everything expect
            ## ## the score function at the scaled fitted means
            ## fit_beta <- key_quantities(theta, y = y, what = "mean", scale_totals = has_fixed_totals, qr = TRUE)
            ## grad_beta <-  gradient(theta, fit = if (has_fixed_totals) NULL else fit_beta, what = "mean")
            ##                     inverse_info_beta <- try(information(theta, inverse = TRUE, fit = fit_beta, what = "mean"))
            ## if (failed_inversion <- inherits(inverse_info_beta, "try-error")) {
            ##     warning("failed to invert the information matrix: iteration stopped prematurely")
            ##     break
            ## }
            ## if (is_ML) {
            ##     adjustment_beta <- 0
            ##     failed_adjustment <- FALSE
            ## }
            ## else {
            ##     adjustment_beta <- adjustment_function(theta, fit = fit_beta, what = "mean")
            ##     if (failed_adjustment <- any(is.na(adjustment_beta))) {
            ##         warning("failed to calculate the bias-reducing score adjustment: iteration stopped prematurely")
            ##         break
            ##     }
            ## }
            ## adjusted_grad_beta <- grad_beta + adjustment_beta
            ## if (no_dispersion) {
            ##         disp <- 1
            ##         transformed_dispersion <- eval(control$Trans)
            ##         adjusted_grad_zeta <- NA
            ##         adjustment_zeta <- NA
            ##         info_dispersion <- NA
            ##         d1zeta <- NA
            ## }
            ## else {
            ##     if (df_residual > 0) {
            ##         theta <- c(betas, disp)
            ##         fit_dispersion <- key_quantities(theta, y = y, what = "dispersion", qr = TRUE)
            ##         gradDispersion <-  gradient(theta, fit = fit_dispersion, what = "dispersion")
            ##         info_dispersion <- information(theta, inverse = FALSE, fit = fit_dispersion, what = "dispersion")
            ##         d1zeta <- eval(d1_transformed_dispersion)
            ##         if (is_ML) {
            ##             adjusted_grad_zeta <- 0
            ##         }
            ##         else {
            ##             adjustmentDispersion <- adjustment_function(theta, fit = fit_dispersion, what = "dispersion")
            ##             ## The adjustment for transDisp (use Kosmidis & Firth, 2010, Remark 3 for derivation)
            ##             ## FIXME: some redundancy below...
            ##             adjustment_zeta <- adjustmentDispersion/d1zeta - 0.5 * eval(d2_transformed_dispersion) / d1zeta^2
            ##             adjusted_grad_zeta <- gradDispersion/d1zeta + adjustment_zeta
            ##         }
            ##     }
            ##     else {
            ##         disp <- 1 # No effect to the adjusted scores
            ##         transformed_dispersion <- eval(control$Trans)
            ##         adjusted_grad_zeta <- NA
            ##     }
            ## }

        }
        else {
            step_beta <- 0
            step_zeta <- 0
            for (iter in seq.int(control$maxit)) {
                step_factor <- 0
                testhalf <- TRUE
                while (testhalf & step_factor < control$maxStepFactor) {
                    step_beta_previous <- step_beta
                    step_zeta_previous <- step_zeta
                    ## Update betas
                    betas <- betas + 2^(-step_factor) * step_beta
                    ## Update transformed_dispersion, if any
                    if (no_dispersion) {
                        disp <- 1
                        transformed_dispersion <- eval(control$Trans)
                        adjusted_grad_zeta <- NA
                        adjustment_zeta <- NA
                        info_dispersion <- NA
                        d1zeta <- NA
                        step_zeta <- 0
                    }
                    else {
                        if (df_residual > 0) {
                            transformed_dispersion <- transformed_dispersion + 2^(-step_factor) * step_zeta
                            disp <- eval(control$inverseTrans)
                        }
                        else {
                            disp <- 1 # No effect to the adjusted scores
                            transformed_dispersion <- eval(control$Trans)
                            adjusted_grad_zeta <- NA
                            step_zeta <- 0
                        }
                    }
                    # Compute key quantities to update direction
                    theta <- c(betas, disp)
                    quantities <- key_quantities(theta, y = y, what = what,
                                                 scale_totals = has_fixed_totals, qr = TRUE)
                    grad_beta <-  gradient(theta, fit = if (has_fixed_totals) NULL else quantities,
                                           what = "mean")
                    inverse_info_beta <- try(information(theta, inverse = TRUE, fit = quantities,
                                                         what = "mean"))
                    if (failed_inversion <- inherits(inverse_info_beta, "try-error")) {
                        warning("failed to invert the information matrix: iteration stopped prematurely")
                        break
                    }
                    ## Step for betas
                    if (is_ML) {
                        adjustment_beta <- 0
                        failed_adjustment <- FALSE
                    }
                    else {
                        adjustment_beta <- adjustment_function(theta, fit = quantities, what = "mean")
                        if (failed_adjustment <- any(is.na(adjustment_beta))) {
                            warning("failed to calculate the adjusted score equations: iteration stopped prematurely")
                            break
                        }
                    }
                    adjusted_grad_beta <- grad_beta + adjustment_beta
                    step_beta <- drop(inverse_info_beta %*% adjusted_grad_beta)
                    ## Step for zetas
                    if (!no_dispersion & df_residual > 0) {
                        gradDispersion <-  gradient(theta, fit = quantities, what = "dispersion")
                        info_dispersion <- information(theta, inverse = FALSE, fit = quantities, what = "dispersion")
                        d1zeta <- eval(d1_transformed_dispersion)
                        if (is_ML) {
                            adjusted_grad_zeta <- 0
                        }
                        else {
                            adjustmentDispersion <- adjustment_function(theta, fit = quantities, what = "dispersion")
                            ## The adjustment for transformed_dispersion (use Kosmidis & Firth, 2010, Remark 3 for derivation)
                            ## FIXME: some redundancy below...
                            adjustment_zeta <- adjustmentDispersion/d1zeta - 0.5 * eval(d2_transformed_dispersion) / d1zeta^2
                            adjusted_grad_zeta <- gradDispersion/d1zeta + adjustment_zeta
                        }
                        step_zeta <- as.vector(d1zeta^2 * adjusted_grad_zeta/info_dispersion)
                    }
                    ## Only test if iteration is > 1 and step_beta has left its starting value of zero
                    if (step_factor == 0 & iter == 1)  {
                        testhalf <- TRUE
                    }
                    else {
                        s2 <- sum(step_beta^2) + step_zeta^2
                        s1 <- sum(step_beta_previous^2) + step_zeta_previous^2
                        testhalf <- s2 > s1
                    }
                    step_factor <- step_factor + 1
                    ##  Trace here
                    if (control$trace) {
                        trace_iteration()
                    }
                }
                if (failed_adjustment | failed_inversion | all(abs(c(abs(step_beta), abs(step_zeta))) < control$epsilon, na.rm = TRUE)) {
                    break
                }
            }

            ## NEW ITERATION
            ## step_beta <- 0
            ## step_zeta <- 0
            ## for (iter in seq.int(control$maxit)) {
            ##     step_factor <- 0
            ##     testhalf <- TRUE
            ##     while (testhalf & step_factor < control$maxStepFactor) {
            ##         step_beta_previous <- step_beta
            ##         ## Update betas
            ##         previous_betas <- betas
            ##         betas <- betas + 2^(-step_factor) * step_beta
            ##         theta <- c(betas, disp)
            ##         fit_beta <- key_quantities(theta, y = y, what = "mean",
            ##                                    scale_totals = has_fixed_totals, qr = TRUE)
            ##         grad_beta <-  gradient(theta, fit = if (has_fixed_totals) NULL else fit_beta,
            ##                                what = "mean")
            ##         inverse_info_beta <- try(information(theta, inverse = TRUE, fit = fit_beta,
            ##                                              what = "mean"))
            ##         if (failed_inversion <- inherits(inverse_info_beta, "try-error")) {
            ##             warning("failed to invert the information matrix: iteration stopped prematurely")
            ##             ## Revert to the previous value for which inversion was possible
            ##             betas <- previous_betas
            ##             break
            ##         }
            ##         if (is_ML) {
            ##             adjustment_beta <- 0
            ##             failed_adjustment <- FALSE
            ##         }
            ##         else {
            ##             adjustment_beta <- adjustment_function(theta, fit = fit_beta, what = "mean")
            ##             if (failed_adjustment <- any(is.na(adjustment_beta))) {
            ##                 warning("failed to calculate the adjusted score equations: iteration stopped prematurely")
            ##                 ## Revert to the previous value for which inversion was possible
            ##                 betas <- previous_betas
            ##                 break
            ##             }
            ##         }
            ##         adjusted_grad_beta <- grad_beta + adjustment_beta
            ##         step_beta <- drop(inverse_info_beta %*% adjusted_grad_beta)
            ##         ## Only test if iteration is > 1 and step_beta has left its starting value of zero
            ##         testhalf <- if (step_factor == 0 & iter == 1)  TRUE else sum(step_beta^2) > sum(step_beta_previous^2)
            ##         step_factor <- step_factor + 1
            ##         if (control$trace) {
            ##             trace_iteration(what = "coefficient")
            ##         }
            ##     }
            ##     ## Bias-reduced estimation of (transformed) dispersion
            ##     if (no_dispersion) {
            ##         disp <- 1
            ##         transformed_dispersion <- eval(control$Trans)
            ##         adjusted_grad_zeta <- NA
            ##         adjustment_zeta <- NA
            ##         info_dispersion <- NA
            ##         d1zeta <- NA
            ##         step_zeta <- 0
            ##     }
            ##     else {
            ##         if (df_residual > 0) {
            ##             transformed_dispersion <- transformed_dispersion + control$slowit * step_zeta
            ##             disp <- eval(control$inverseTrans)
            ##             theta <- c(betas, disp)
            ##             fit_dispersion <- key_quantities(theta, y = y, what = "dispersion", qr = TRUE)
            ##             gradDispersion <-  gradient(theta, fit = fit_dispersion, what = "dispersion")
            ##             info_dispersion <- information(theta, inverse = FALSE, fit = fit_dispersion, what = "dispersion")
            ##             d1zeta <- eval(d1_transformed_dispersion)
            ##             if (is_ML) {
            ##                 adjusted_grad_zeta <- 0
            ##             }
            ##             else {
            ##                 adjustmentDispersion <- adjustment_function(theta, fit = fit_dispersion, what = "dispersion")
            ##                 ## The adjustment for transformed_dispersion (use Kosmidis & Firth, 2010, Remark 3 for derivation)
            ##                 ## FIXME: some redundancy below...
            ##                 adjustment_zeta <- adjustmentDispersion/d1zeta - 0.5 * eval(d2_transformed_dispersion) / d1zeta^2
            ##                 adjusted_grad_zeta <- gradDispersion/d1zeta + adjustment_zeta
            ##             }
            ##             step_zeta <- as.vector(d1zeta^2 * adjusted_grad_zeta/info_dispersion)
            ##         }
            ##         else {
            ##             disp <- 1 # No effect to the adjusted scores
            ##             transformed_dispersion <- eval(control$Trans)
            ##             adjusted_grad_zeta <- NA
            ##             step_zeta <- 0
            ##         }
            ##     }
            ##     if (control$trace) {
            ##         trace_iteration(what = "dispersion")
            ##     }
            ##     if (failed_adjustment | failed_inversion | all(abs(c(abs(step_beta), abs(step_zeta))) < control$epsilon, na.rm = TRUE)) {
            ##         break
            ##     }
            ## }

            ## ## CAREFUL: Old iteration below
            ## for (iter in seq.int(control$maxit)) {
            ##     step_factor <- 0
            ##     testhalf <- TRUE
            ##     betasPrev <- betas
            ##     objPrev <- objCur
            ##     ## Bias-reduced estimation of mean effects
            ##     while (testhalf & step_factor < control$maxStepFactor) {
            ##         theta <- c(betas, disp)
            ##         fit_beta <- key_quantities(theta, y = y, what = "mean",
            ##                                        scale_totals = has_fixed_totals, qr = TRUE)
            ##         grad_beta <-  gradient(theta, fit = if (has_fixed_totals) NULL else fit_beta,
            ##                                    what = "mean")
            ##         inverse_info_beta <- try(information(theta, inverse = TRUE, fit = fit_beta,
            ##                                               what = "mean"))
            ##         if (failed_inversion <- inherits(inverse_info_beta, "try-error")) {
            ##             warning("failed to invert the information matrix: iteration stopped prematurely")
            ##             break
            ##         }
            ##         if (is_ML) {
            ##             adjustment_beta <- 0
            ##             failed_adjustment <- FALSE
            ##         }
            ##         else {
            ##             adjustment_beta <- adjustment_function(theta, fit = fit_beta, what = "mean")
            ##             if (failed_adjustment <- any(is.na(adjustment_beta))) {
            ##                 warning("failed to calculate the bias-reducing score adjustment: iteration stopped prematurely")
            ##                 break
            ##             }
            ##         }
            ##         adjusted_grad_beta <- grad_beta + adjustment_beta
            ##         step_beta <- drop(inverse_info_beta %*% adjusted_grad_beta)
            ##         betas <- betasPrev + 2^(-step_factor) * step_beta
            ##         objCur <- sum(step_beta^2)
            ##         step_factor <- step_factor + 1
            ##         testhalf <- objCur > objPrev
            ##         if (control$trace) {
            ##             trace_iteration(what = "coefficient")
            ##         }
            ##     }
            ##     ## Bias-reduced estimation of (transformed) dispersion
            ##     if (no_dispersion) {
            ##         disp <- 1
            ##         transformed_dispersion <- eval(control$Trans)
            ##         adjusted_grad_zeta <- NA
            ##         adjustment_zeta <- NA
            ##         info_dispersion <- NA
            ##         d1zeta <- NA
            ##         step_zeta <- 0
            ##     }
            ##     else {
            ##         if (df_residual > 0) {
            ##             theta <- c(betas, disp)
            ##             fit_dispersion <- key_quantities(theta, y = y, what = "dispersion", qr = TRUE)
            ##             gradDispersion <-  gradient(theta, fit = fit_dispersion, what = "dispersion")
            ##             info_dispersion <- information(theta, inverse = FALSE, fit = fit_dispersion, what = "dispersion")
            ##             d1zeta <- eval(d1_transformed_dispersion)
            ##             if (is_ML) {
            ##                 adjusted_grad_zeta <- 0
            ##             }
            ##             else {
            ##                 adjustmentDispersion <- adjustment_function(theta, fit = fit_dispersion, what = "dispersion")
            ##                 ## The adjustment for transformed_dispersion (use Kosmidis & Firth, 2010, Remark 3 for derivation)
            ##                 ## FIXME: some redundancy below...
            ##                 adjustment_zeta <- adjustmentDispersion/d1zeta - 0.5 * eval(d2_transformed_dispersion) / d1zeta^2
            ##                 adjusted_grad_zeta <- gradDispersion/d1zeta + adjustment_zeta
            ##             }
            ##             step_zeta <- as.vector(d1zeta^2 * adjusted_grad_zeta/info_dispersion)
            ##             transformed_dispersion <- transformed_dispersion + control$slowit * step_zeta
            ##             disp <- eval(control$inverseTrans)
            ##         }
            ##         else {
            ##             disp <- 1 # No effect to the adjusted scores
            ##             transformed_dispersion <- eval(control$Trans)
            ##             adjusted_grad_zeta <- NA
            ##             step_zeta <- 0
            ##         }
            ##     }
            ##     if (control$trace) {
            ##         trace_iteration(what = "dispersion")
            ##     }
            ##     if (failed_adjustment | failed_inversion | all(abs(c(abs(step_beta), abs(step_zeta))) < control$epsilon, na.rm = TRUE)) {
            ##         break
            ##     }
            ## }
            ## ## END old iterations

        }

        ## Objects used from above iteration
        ## adjusted_grad_beta
        ## adjusted_grad_zeta
        ## betas
        ## disp

        adjusted_grad_all[betas_names] <- adjusted_grad_beta
        adjusted_grad_all["Transformed dispersion"] <- adjusted_grad_zeta
        betas_all[betas_names] <- betas

        ## Convergence analysis
        if ((failed_inversion | failed_adjustment | iter >= control$maxit) & !(is_correction)) {
            warning("brglmFit: algorithm did not converge", call. = FALSE)
            converged <- FALSE
        }
        else {
            converged <- TRUE
        }

        if (boundary) {
            warning("brglmFit: algorithm stopped at boundary value", call. = FALSE)
        }

        ## QR decomposition and fitted values are at the final value
        ## for the coefficients

        ## QR decomposition for cov.unscaled
        if (!isTRUE(is_full_rank)) {
            x <- X_all
            betas <- betas_all
            betas[is.na(betas)] <- 0
            nvars <- nvars_all
        }

        ## If has_fixed_totals = TRUE, then scale fitted values before
        ## calculating QR decompositions, fitted values, etas,
        ## residuals and working_weights

        qr.Wx <- quantities$qr_decomposition

        mus <- quantities$mus
        etas <- quantities$etas
        ## Residuals
        residuals <- with(quantities, (y - mus)/d1mus)
        working_weights <- quantities$working_weights

        ## Fisher information for the transformed dispersion
        d1zeta <- eval(d1_transformed_dispersion)
        if (!no_dispersion) {
            info_transformed_dispersion <- info_dispersion/d1zeta^2
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
        if (df_residual == 0) disp <- NaN

        ## Estimate of first-order bias from the last iteration (so
        ## not at the final value for the coefficients)
        if (is_ML) {
            ## For now... To be set to calculate biases at a later version
            bias_betas <- bias_zeta <- NULL
        }
        else {
            bias_betas <- -drop(inverse_info_beta %*% adjustment_beta)
            bias_zeta <- -adjustment_zeta/info_dispersion/d1zeta^2
            bias_betas_all <- betas_all
            bias_betas_all[betas_names] <- bias_betas
            ## If correction has been requested then add estimated biases an attribute to the coefficients
            if (is_correction) {
                attr(betas_all, "biases") <- bias_betas_all
                attr(transformed_dispersion, "biases") <- bias_zeta
            }
        }
    }

    ## Working weights
    wt <- rep.int(0, nobs)
    wt[keep] <- working_weights[keep]
    names(wt) <- names(residuals) <- names(mus) <- names(etas) <- names(weights) <- names(y) <- ynames
    ## For the null deviance:
    ##
    ## If there is an intercept but not an offset then the ML fitted
    ## value is the weighted average and is calculated easily below if
    ## ML is used
    ##
    control0 <- control
    if (customTransformation) {
        control0$transformation <- transformation0
    }
    if (intercept & missingOffset) {
        nullFit <- brglmFit(x = x[, "(Intercept)"], y = y, weights = weights,
                            offset = rep(0, nobs), family = family, intercept = TRUE,
                            control = control0[c("epsilon", "maxit", "type", "transformation", "slowit")])
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

    list(coefficients = betas_all,
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
         df.residual = df_residual,
         df.null = nulldf,
         y = y,
         converged = converged,
         boundary = boundary,
         dispersion = disp,
         dispersionML = dispML,
         transformed_dispersion = transformed_dispersion,
         info_transformed_dispersion = if (no_dispersion) NA else info_transformed_dispersion,
         grad = adjusted_grad_all,
         transformation = control$transformation,
         ## cov.unscaled = tcrossprod(R_matrix),
         type = control$type,
         class = "brglmFit")
}

#' @export
coef.brglmFit <- function(object, model = c("mean", "full", "dispersion"), ...) {
    model <- match.arg(model)
    switch(model,
           "mean" = {
               object$coefficients
           },
           "dispersion" = {
               transDisp <- object$transformed_dispersion
               names(transDisp) <- paste0(object$transformation, "(dispersion)")
               transDisp
               ## This will ALWAYS be on the scale of the TRANSFORMED dispersion
           },
           "full" = {
               transDisp <- object$transformed_dispersion
               ntd <- paste0(object$transformation, "(dispersion)")
               names(transDisp) <- ntd
               betas <- object$coefficients
               thetaTrans <- c(betas, transDisp)
               if (object$type == "correction") {
                   bcf <- attr(betas, "biases")
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
               vtd <- 1/object$info_transformed_dispersion
               ntd <- paste0(object$transformation, "(dispersion)")
               names(vtd) <- ntd
               vtd
           },
           full = {
               vbetas <- summary.brglmFit(object, ...)$cov.scaled
               vtd <- 1/object$info_transformed_dispersion
               nBetasAll <- c(rownames(vbetas), paste0(object$transformation, "(dispersion)"))
               vBetasAll <- cbind(rbind(vbetas, 0),
                                  c(numeric(nrow(vbetas)), vtd))
               dimnames(vBetasAll) <- list(nBetasAll, nBetasAll)
               vBetasAll
           })
}


DD <- function(expr,name, order = 1) {
    if(order < 1) stop("'order' must be >= 1")
    if(order == 1) D(expr,name)
    else DD(D(expr, name), name, order - 1)
}



