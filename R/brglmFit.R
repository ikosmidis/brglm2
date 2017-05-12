# Copyright (C) 2016, 2017 Ioannis Kosmidis
# function `AS_median_adjustment`: Copyright (C) 2017, Euloge Clovis Kenne Pagui, Ioannis Kosmidis

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


#' Fitting function for \code{\link{glm}} for reduced-bias
#' estimation and inference
#'
#' \code{\link{brglmFit}} is a fitting function for \code{\link{glm}}
#' that fits generalized linear models using implicit and explicit
#' bias reduction methods. Currently supported methods include the
#' mean bias-reducing adjusted scores approach in Firth (1993) and
#' Kosmidis \& Firth (2009), the median bias-reduction adjusted scores
#' approach in Kenne et al. (2016), the correction of the asymptotic
#' bias in Cordeiro & McCullagh (1991), and maximum likelihood.
#' Estimation is performed using a quasi Fisher scoring iteration,
#' which in the case of mean-bias reduction resembles an iterative
#' correction of the asymptotic bias of the Fisher scoring iterates.
#'
#' @inheritParams stats::glm.fit
#' @aliases brglm_fit
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
#' iterates. A mathematical description of the supported adjustments
#' and the quasi Fisher scoring iteration is provided in the iteration
#' vignette (see,
#' \url{https://cran.r-project.org/package=brglm2/vignettes/iteration.pdf}).
#' A quick description of the quasi Fisher scoring iteration is also
#' given in one of the vignettes of the *enrichwith* R package (see,
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
#' \code{fixed_totals} can be used to constrain the means of a poisson
#' model to add up to the corresponding observed counts according to
#'
#' \code{brglm_fit} is an alias to \code{brglmFit}.
#'
#' @author Ioannis Kosmidis [aut, cre] \email{i.kosmidis@ucl.ac.uk}, Euloge Clovis Kenne Pagui [ctb] \email{kenne@stat.unipd.it}
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
#' Kenne Pagui E C, Salvan A and Sartori N (2016). Median bias
#' reduction of maximum likelihood estimates. *arXiv*,
#' **arXiv:1604.04768**
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
#' @examples
#' ## The lizards example from ?brglm::brglm
#' data("lizards")
#' # Fit the model using maximum likelihood
#' lizardsML <- glm(cbind(grahami, opalinus) ~ height + diameter +
#'                  light + time, family = binomial(logit), data = lizards,
#'                  method = "glm.fit")
#' # Mean bias-reduced fit:
#' lizardsBR_mean <- glm(cbind(grahami, opalinus) ~ height + diameter +
#'                       light + time, family = binomial(logit), data = lizards,
#'                       method = "brglmFit")
#' # Median bias-reduced fit:
#' lizardsBR_median <- glm(cbind(grahami, opalinus) ~ height + diameter +
#'                         light + time, family = binomial(logit), data = lizards,
#'                         method = "brglmFit", type = "AS_median")
#' summary(lizardsML)
#' summary(lizardsBR_median)
#' summary(lizardsBR_mean)
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
#' # The mean bias-reduced fit
#' coalitionBR_mean <- update(coalitionML, method = "brglmFit")
#' # The bias-corrected fit
#' coalitionBC <- update(coalitionML, method = "brglmFit", type = "correction")
#' # The median bias-corrected fit
#' coalitionBR_median <- update(coalitionML, method = "brglmFit", type = "AS_median")
#' }
#'
#' \dontrun{
#' ## An example with offsets from Venables & Ripley (2002, p.189)
#' data("anorexia", package = "MASS")
#'
#' anorexML <- glm(Postwt ~ Prewt + Treat + offset(Prewt),
#'                 family = gaussian, data = anorexia)
#' anorexBC <- update(anorexML, method = "brglmFit", type = "correction")
#' anorexBR_mean <- update(anorexML, method = "brglmFit")
#' anorexBR_median <- update(anorexML, method = "brglmFit", type = "AS_median")
#'
#' ## All methods return the same estimates for the regression
#' ## parameters because the maximum likelihood estimator is normally
#' ## distributed around the `true` value under the model (hence, both
#' ## mean and component-wise median unbiased). The Wald tests for
#' ## anorexBC and anorexBR_mean differ from anorexML
#' ## because the bias-reduced estimator of the dispersion is the
#' ## unbiased, by degree of freedom adjustment (divide by n - p),
#' ## estimator of the residual variance. The Wald tests from
#' ## anorexBR_median are based on the median bias-reduced estimator
#' ## of the dispersion that results from a different adjustment of the
#' ## degrees of freedom (divide by n - p - 2/3)
#' summary(anorexML)
#' summary(anorexBC)
#' summary(anorexBR_mean)
#' summary(anorexBR_median)
#' }
#'
#' ## endometrial data from Heinze \& Schemper (2002) (see ?endometrial)
#' data("endometrial", package = "brglm2")
#' endometrialML <- glm(HG ~ NV + PI + EH, data = endometrial,
#'                      family = binomial("probit"))
#' endometrialBR_mean <- update(endometrialML, method = "brglmFit",
#'                              type = "AS_mean")
#' endometrialBC <- update(endometrialML, method = "brglmFit",
#'                         type = "correction")
#' endometrialBR_median <- update(endometrialML, method = "brglmFit",
#'                                type = "AS_median")
#' summary(endometrialML)
#' summary(endometrialBC)
#' summary(endometrialBR_mean)
#' summary(endometrialBR_median)
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
            cat("Outer/Inner iteration:\t", sprintf("%03d", iter), "/", sprintf("%03d", step_factor), "\n", sep = "")
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
    key_quantities <- function(pars, y, level = 0, scale_totals = FALSE, qr = TRUE) {
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
            out$d1varmus <- d1variance(mus)
            out$working_weights <- working_weights
            if (qr) out$qr_decomposition <- qr(wx)
            out
        }
        dispersion_quantities <- function(out) {
            zetas <- -weights * prec
            out$zetas <- zetas
            ## Evaluate the derivatives of the a function only for
            ## objervations with non-zero weight
            d1afuns <- d2afuns <- d3afuns <- rep(NA_real_, nobs)
            d1afuns[keep] <- d1afun(zetas[keep])
            d2afuns[keep] <- d2afun(zetas[keep])
            d3afuns[keep] <- d3afun(zetas[keep])
            out$d2afuns <- d2afuns
            out$d3afuns <- d3afuns
            out$deviance_residuals <- dev.resids(y, mus, weights)
            out$Edeviance_residuals <- weights * d1afuns
            out
        }
        if (level == 0) {
            out <- mean_quantities(out)
        }
        if (level == 1) {
            out <- dispersion_quantities(out)
        }
        if (level > 1) {
            out <- mean_quantities(out)
            out <- dispersion_quantities(out)
        }
        out
    }

    gradient <- function(pars, level = 0, fit = NULL) {
        if (is.null(fit)) {
            fit <- key_quantities(pars, y = y, level = level, qr = FALSE)
        }
        with(fit, {
            if (level == 0) {
                score_components <- weights * d1mus  * (y - mus) / varmus * x
                return(precision * .colSums(score_components, nobs, nvars, TRUE))
            }
            if (level == 1) {
                return(1/2 * precision^2 * sum(deviance_residuals - Edeviance_residuals, na.rm = TRUE))
            }
        })
    }

    information <- function(pars, level = 0, fit = NULL, inverse = FALSE) {
        if (is.null(fit)) {
            fit <- key_quantities(pars, y = y, level = level, qr = TRUE)
        }
        with(fit, {
            if (level == 0) {
                R_matrix <- qr.R(qr_decomposition)
                if (inverse) {
                    ## return(dispersion * tcrossprod(solve(R_matrix)))
                    return(dispersion * chol2inv(R_matrix))
                }
                else {
                    return(precision * crossprod(R_matrix))
                }
            }
            if (level == 1) {
                info <- 0.5 * sum(weights^2 * d2afuns, na.rm = TRUE)/dispersion^4
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
            fit <- key_quantities(pars, y = y, level = 0, qr = TRUE)
        }
        with(fit, {
            Qmat <- qr.Q(qr_decomposition)
            .rowSums(Qmat * Qmat, nobs, nvars, TRUE)
        })
    }

    ## FIXME: Redundant function for now
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
                    cfit <- key_quantities(theta, y = y, level = 1, qr = FALSE)
                    gradient(theta, level = 1, fit = cfit)
                }, lower = .Machine$double.eps, upper = 10000, tol = control$epsilon), silent = FALSE)
                if (inherits(dispFit, "try-error")) {
                    warning("the ML estimate of the dispersion could not be calculated. An alternative estimate had been used as starting value.")
                    dispML <- NA_real_
                    disp <- NA_real_
                }
                else {
                    disp <- dispML <- dispFit$root
                }
            }
            else { ## if the model is saturated dispML is NA_real_
                disp <- 1 ## A convenient value
                dispML <- NA_real_
            }
        }
        list(dispersion = disp, dispersion_ML = dispML)
    }

    AS_mean_adjustment <- function(pars, level = 0, fit = NULL) {
        if (is.null(fit)) {
            fit <- key_quantities(pars, y = y, level = level, qr = TRUE)
        }
        with(fit, {
            if (level == 0) {
                hatvalues <- hat_values(pars, fit = fit)
                ## Use only observations with keep = TRUE to ensure that no division with zero takes place
                return(.colSums(0.5 * hatvalues * d2mus/d1mus * x, nobs, nvars, TRUE))
            }
            if (level == 1) {
                s1 <- sum(weights^3 * d3afuns, na.rm = TRUE)
                s2 <- sum(weights^2 * d2afuns, na.rm = TRUE)
                return((nvars - 2)/(2 * dispersion) + s1/(2 * dispersion^2 * s2))
            }
        })
    }

    ## Implementation by Euloge Clovis Kenne Pagui, 20 April 2017 (kept here for testing)
    ## AS_median_adjustment <- function(pars, level = 0, fit = NULL) {
    ##     if (is.null(fit)) {
    ##         fit <- key_quantities(pars, y = y, level = level, qr = TRUE)
    ##     }
    ##     with(fit, {
    ##         if (level == 0) {
    ##             R_matrix <- qr.R(qr_decomposition)
    ##             info <- precision * crossprod(R_matrix)
    ##             inverse_info <- try(dispersion * tcrossprod(solve(R_matrix)))
    ##             nu_r_s_t <- nu_r_st <- array(0,c(nvars,nvars,nvars))
    ##             k1 <- k2 <- k3 <- b.beta <- rep(NA,nvars)
    ##             for (r in 1:nvars)
    ##             {
    ##               nu_r_s_t[r,,] <- t(x)%*%((working_weights*d1mus*d1varmus*x[,r]/varmus)*x)
    ##               nu_r_st[r,,] <- -t(x)%*%((working_weights*d1mus*(d1varmus/varmus-d2mus/d1mus^2)*x[,r])*x)
    ##             }
    ##             k2 <- 1/diag(inverse_info)
    ##             for (r in 1:nvars)
    ##             {
    ##               sum_s1 <- rep(0,nvars)
    ##               sum_s3 <- rep(0,nvars)
    ##               nu.tu <- inverse_info-outer(inverse_info[,r]*k2[r],inverse_info[,r])
    ##               for (s in 1:nvars){
    ##                 sum_s1[s] <- sum(diag(nu.tu%*%(nu_r_s_t[s,,]+nu_r_st[s,,])))
    ##                 sum_s3[s] <- sum((inverse_info[r,]%*%nu_r_s_t[s,,])*inverse_info[r,])
    ##               }
    ##               barb1r <- sum(sum_s1*inverse_info[r,])
    ##               barb2r <- k2[r]*sum(sum_s3*inverse_info[r,])
    ##               b.beta[r]<- barb1r/2+barb2r/6
    ##             }
    ##             return(info%*%b.beta)
    ##         }
    ##         if (level == 1) {
    ##               s1 <- sum(weights^3 * d3afuns, na.rm = TRUE)
    ##               s2 <- sum(weights^2 * d2afuns, na.rm = TRUE)
    ##               return(nvars/(2*dispersion) + s1/(6*dispersion^2*s2))
    ##         }
    ##     })
    ## }

    ## Implementation by Ioannis Kosmidis, 02 May 2017
    AS_median_adjustment <- function(pars, level = 0, fit = NULL) {
        if (is.null(fit)) {
            fit <- key_quantities(pars, y = y, level = level, qr = TRUE)
        }
        with(fit, {
            if (level == 0) {
                hatvalues <- hat_values(pars, fit = fit)
                R_matrix <- qr.R(qr_decomposition)
                info_unscaled <- crossprod(R_matrix)
                inverse_info_unscaled <- chol2inv(R_matrix)
                ## FIXME: There is 1) definitely a better way to do this, 2) no time...
                b_vector <- numeric(nvars)
                for (j in seq.int(nvars)) {
                    inverse_info_unscaled_j <- inverse_info_unscaled[j, ]
                    vcov_j <- tcrossprod(inverse_info_unscaled_j) / inverse_info_unscaled_j[j]
                    hats_j <- .rowSums((x %*% vcov_j) * x, nobs, nvars, TRUE) * working_weights
                    b_vector[j] <- inverse_info_unscaled_j %*% .colSums(x * (hats_j * (d1mus * d1varmus / (6 * varmus) - 0.5 * d2mus/d1mus)), nobs, nvars, TRUE)
                }
                return(.colSums(0.5 * hatvalues * d2mus / d1mus * x, nobs, nvars, TRUE) +
                       info_unscaled %*% b_vector)
            }
            if (level == 1) {
                s1 <- sum(weights^3 * d3afuns, na.rm = TRUE)
                s2 <- sum(weights^2 * d2afuns, na.rm = TRUE)
                return(nvars/(2 * dispersion) + s1/(6 * dispersion^2 * s2))
            }
        })
    }

    customTransformation <- is.list(control$transformation) & length(control$transformation) == 2
    if (customTransformation) {
        transformation0 <- control$transformation
    }

    control <- do.call("brglmControl", control)

    ## FIXME: Add IBLA
    adjustment_function <- switch(control$type,
                            "correction" = AS_mean_adjustment,
                            "AS_mean" = AS_mean_adjustment,
                            "AS_median" = AS_median_adjustment,
                            "ML" = function(pars, ...) 0)

    ## compute_step_components does everything on the scale of the /transformed/ dispersion
    compute_step_components <- function(pars, level = 0, fit = NULL) {
        if (is.null(fit)) {
            fit <- key_quantities(pars, y = y, level = level, qr = TRUE)
        }
        if (level == 0) {
            grad <-  gradient(pars, fit = if (has_fixed_totals) NULL else fit, level = 0)
            inverse_info <- try(information(pars, inverse = TRUE, fit = fit, level = 0))
            failed_inversion <- inherits(inverse_info, "try-error")
            adjustment <- adjustment_function(pars, fit = fit, level = 0)
            failed_adjustment <- any(is.na(adjustment))
        }
        if (level == 1) {
            if (no_dispersion | df_residual < 1) {
                grad <- adjustment <- inverse_info <- NA_real_
                failed_adjustment <- failed_inversion <- FALSE
            }
            else {
                d1zeta <- eval(d1_transformed_dispersion)
                d2zeta <- eval(d2_transformed_dispersion)
                grad <-  gradient(theta, fit = fit, level = 1)/d1zeta
                inverse_info <- 1/information(theta, inverse = FALSE, fit = fit, level = 1) * d1zeta^2
                failed_inversion <- !is.finite(inverse_info)
                ## adjustment <- adjustment_function(theta, fit = fit, level = 1)/d1zeta - if (is_ML | is_AS_median) 0 else 0.5 * d2zeta / d1zeta^2
                adjustment <- adjustment_function(theta, fit = fit, level = 1)/d1zeta - 0.5 * d2zeta / d1zeta^2
                failed_adjustment <- is.na(adjustment)
            }
        }
        out <- list(grad = grad,
                    inverse_info = inverse_info,
                    adjustment = adjustment,
                    failed_adjustment = failed_adjustment,
                    failed_inversion = failed_inversion)
        out
    }


    ## Some useful quantities
    is_ML <- control$type == "ML"
    is_AS_median <- control$type == "AS_median"
    is_correction <- control$type == "correction"
    no_dispersion <- family$family %in% c("poisson", "binomial")


    if (is_AS_median) {
        transformation1 <- control$transformation
        Trans1 <- control$Trans
        inverseTrans1 <- control$inverseTrans
        control$transformation <- "identity"
        control$Trans <- expression(dispersion)
        control$inverseTrans <- expression(transformed_dispersion)
    }


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
    if (missing_offset <- is.null(offset)) {
        offset <- rep.int(0, nobs)
    }

    ## Enrich the family object with the required derivatives
    linkglm <- make.link(family$link)
    family <- enrichwith::enrich(family, with = c("d1afun", "d2afun", "d3afun", "d1variance"))
    linkglm <- enrichwith::enrich(linkglm, with = "d2mu.deta")

    ## Extract functions from the enriched family object
    variance <- family$variance
    d1variance <- family$d1variance
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
    d1_transformed_dispersion <- DD(control$Trans, "dispersion", order = 1)
    d2_transformed_dispersion <- DD(control$Trans, "dispersion", order = 2)


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
        betas_all <- structure(rep(NA_real_, nvars_all), .Names = betas_names_all)
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
            dispersion <- dispList$dispersion
            if (is.na(dispersion)) dispersion <- var(y)/variance(sum(weights * y)/sum(weights))
            dispersion_ML <- dispList$dispersion_ML
            transformed_dispersion <- eval(control$Trans)
        }
        else {
            if ((length(start) == nvars_all) & is.numeric(start)) {
                betas_all <- start
                names(betas_all) <- betas_names_all
                if (!isTRUE(is_full_rank)) {
                    betas_all[aliased] <- NA_real_
                    betas <- betas_all[-aliased]
                }
                else {
                    betas <- betas_all
                }
                ## Estimate dispersion based on current value for betas
                dispList <- estimate_dispersion(betas, y = y)
                dispersion <- dispList$dispersion
                if (is.na(dispersion)) dispersion <- var(y)/variance(sum(weights * y)/sum(weights))
                dispersion_ML <- dispList$dispersion_ML
                transformed_dispersion <- eval(control$Trans)
            }
            if ((length(start) == nvars_all + 1) & is.numeric(start)) {
                betas_all <- start[seq.int(nvars_all)]
                names(betas_all) <- betas_names_all
                if (!isTRUE(is_full_rank)) {
                    betas_all[aliased] <- NA_real_
                    betas <- betas_all[-aliased]
                }
                else {
                    betas <- betas_all
                }
                transformed_dispersion <- start[nvars_all + 1]
                dispersion_ML <- NA_real_
                dispersion <- eval(control$inverseTrans)
            }
            if (length(start) > nvars_all + 1 | length(start) < nvars_all) {
                stop(paste(paste(gettextf("length of 'start' should be equal to %d and correspond to initial betas for %s", nvars_all, paste(deparse(betas_names_all), collapse = ", "), "or", gettextf("to %d and also include a starting value for the transformed dispersion", nvars_all + 1)))), domain = NA_real_)
            }
        }

        adjusted_grad_all <- rep(NA_real_, nvars_all + 1)
        names(adjusted_grad_all) <- c(betas_names_all, "Transformed dispersion")
        if (is_correction) {
            if (control$maxit > 0) control$maxit <- 1
            control$slowit <- 1
            control$max_step_factor <- 1
        }

        ## Evaluate at the starting values
        theta <- c(betas, dispersion)
        transformed_dispersion <- eval(control$Trans)
        ## Mean quantities
        ## If fixed_totals is provided (i.e. multinomial regression
        ## via the Poisson trick) then evaluate everything expect
        ## the score function at the scaled fitted means
        quantities <- key_quantities(theta, y = y, level = 2 * !no_dispersion, scale_totals = has_fixed_totals, qr = TRUE)
        step_components_beta <- compute_step_components(theta, level = 0, fit = quantities)
        step_components_zeta <- compute_step_components(theta, level = 1, fit = quantities)
        if (step_components_beta$failed_inversion) {
            warning("failed to invert the information matrix")
        }
        if (step_components_beta$failed_adjustment) {
            warning("failed to calculate score adjustment")
        }
        adjusted_grad_beta <- with(step_components_beta, {
            grad + adjustment
        })
        step_beta <- drop(step_components_beta$inverse_info %*% adjusted_grad_beta)
        ## Dispersion quantities
        if (no_dispersion) {
            adjusted_grad_zeta <- step_zeta <- NA_real_
        }
        else {
            if (step_components_zeta$failed_inversion) {
                warning("failed to invert the information matrix")
            }
            if (step_components_zeta$failed_adjustment) {
                warning("failed to calculate score adjustment")
            }
            adjusted_grad_zeta <- with(step_components_zeta, {
                grad + adjustment
            })
            step_zeta <- as.vector(adjusted_grad_zeta * step_components_zeta$inverse_info)
        }

        ## Main iterations
        slowit <- control$slowit
        if (control$maxit == 0) {
            iter <- 0
            failed <- FALSE
        }
        else {
            ## Outer iteration
            for (iter in seq.int(control$maxit)) {
                step_factor <- 0
                testhalf <- TRUE
                ## Inner iteration
                while (testhalf & step_factor < control$max_step_factor) {
                    step_beta_previous <- step_beta
                    step_zeta_previous <- step_zeta
                    ## Update betas
                    betas <- betas + slowit * 2^(-step_factor) * step_beta
                    ## Update zetas
                    if (!no_dispersion & df_residual > 0) {
                        transformed_dispersion <- eval(control$Trans)
                        transformed_dispersion <- transformed_dispersion + 2^(-step_factor) * step_zeta
                        dispersion <- eval(control$inverseTrans)
                    }
                    ## Compute key quantities
                    theta <- c(betas, dispersion)
                    transformed_dispersion <- eval(control$Trans)
                    ## Mean quantities
                    quantities <- key_quantities(theta, y = y, level = 2 * !no_dispersion, scale_totals = has_fixed_totals, qr = TRUE)
                    step_components_beta <- compute_step_components(theta, level = 0, fit = quantities)
                    step_components_zeta <- compute_step_components(theta, level = 1, fit = quantities)
                    if (failed_inversion_beta <- step_components_beta$failed_inversion) {
                        warning("failed to invert the information matrix")
                        break
                    }
                    if (failed_adjustment_beta <- step_components_beta$failed_adjustment) {
                        warning("failed to calculate score adjustment")
                        break
                    }
                    adjusted_grad_beta <- with(step_components_beta, {
                        grad + adjustment
                    })
                    step_beta <- drop(step_components_beta$inverse_info %*% adjusted_grad_beta)
                    ## Dispersion quantities
                    if (no_dispersion) {
                        adjusted_grad_zeta <- step_zeta <- NA_real_
                        failed_inversion_zeta <- failed_adjustment_zeta <- FALSE
                    }
                    else {
                        if (failed_inversion_zeta <- step_components_zeta$failed_inversion) {
                            warning("failed to invert the information matrix")
                            break
                        }
                        if (failed_adjustment_zeta <- step_components_zeta$failed_adjustment) {
                            warning("failed to calculate score adjustment")
                            break
                        }
                        adjusted_grad_zeta <- with(step_components_zeta, {
                            grad + adjustment
                        })
                        step_zeta <- as.vector(adjusted_grad_zeta * step_components_zeta$inverse_info)
                    }
                    ## Continue inner loop
                    if (step_factor == 0 & iter == 1)  {
                        testhalf <- TRUE
                    }
                    else {
                        s2 <- c(abs(step_beta), abs(step_zeta))
                        s1 <- c(abs(step_beta_previous), abs(step_zeta_previous))
                        testhalf <- sum(s2, na.rm = TRUE) > sum(s1, na.rm = TRUE)
                    }
                    step_factor <- step_factor + 1
                    ##  Trace here
                    if (control$trace) {
                        trace_iteration()
                    }
                }
                failed <- failed_adjustment_beta | failed_inversion_beta | failed_adjustment_zeta | failed_inversion_zeta
                if (failed | sum(abs(c(step_beta, step_zeta)), na.rm = TRUE) < control$epsilon) {
                    break
                }
            }
        }

        adjusted_grad_all[betas_names] <- adjusted_grad_beta
        adjusted_grad_all["Transformed dispersion"] <- adjusted_grad_zeta
        betas_all[betas_names] <- betas

        ## Convergence analysis
        if ((failed | iter >= control$maxit) & !(is_correction)) {
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
        ## d1zeta <- eval(d1_transformed_dispersion)

        if (!no_dispersion) {
            info_transformed_dispersion <- 1/step_components_zeta$inverse_info
            if (is_AS_median) {
                transformed_dispersion <- eval(Trans1)
                d1zeta <- eval(DD(Trans1, "dispersion", order = 1))
                adjusted_grad_all["Transformed dispersion"] <- adjusted_grad_all["Transformed dispersion"] / d1zeta
                info_transformed_dispersion <- info_transformed_dispersion / d1zeta^2
            }
        }
        if (is_AS_median) {
            control$transformation <- transformation1
            control$trans <- Trans1
            control$inverseTrans <- inverseTrans1
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
        if (df_residual == 0) dispersion <- NA_real_

        ## ## Estimate of first-order bias from the last iteration (so
        ## ## not at the final value for the coefficients)
        ## if (is_ML) {
        ##     ## For now... To be set to calculate biases at a later version
        ##     bias_betas <- bias_zeta <- NULL
        ## }
        ## else {
        ##     bias_betas <- with(step_components_beta, -drop(inverse_info %*% adjustment))
        ##     bias_zeta <- with(step_components_zeta, -drop(inverse_info %*% adjustment))
        ##     bias_betas_all <- betas_all
        ##     bias_betas_all[betas_names] <- bias_betas
        ##     ## If correction has been requested then add estimated biases an attribute to the coefficients
        ##     if (is_correction) {
        ##         attr(betas_all, "biases") <- bias_betas_all
        ##         attr(transformed_dispersion, "biases") <- bias_zeta
        ##     }
        ## }
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
    control0$maxit <- 1000
    if (customTransformation) {
        control0$transformation <- transformation0
    }
    if (intercept & missing_offset) {
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
    if (intercept & !missing_offset) {
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
         dispersion = dispersion,
         dispersion_ML = dispersion_ML,
         transformed_dispersion = transformed_dispersion,
         info_transformed_dispersion = if (no_dispersion) NA_real_ else info_transformed_dispersion,
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
               ## if (object$type == "correction") {
               ##     bcf <- attr(betas, "biases")
               ##     btd <- attr(transDisp, "biases")
               ##     names(btd) <- ntd
               ##     attr(thetaTrans, "biases") <- c(bcf, btd)
               ## }
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



