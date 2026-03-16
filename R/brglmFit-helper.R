# Copyright (C) 2016- Ioannis Kosmidis

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

# Helper functions for brglmFit
# These are self-contained functions that don't rely on parent scope

#' Compute all quantities needed for a GLM fit
#'
#' @param pars Vector of parameters (betas and dispersion)
#' @param y Response vector
#' @param x Design matrix
#' @param weights Observation weights
#' @param offset Offset vector
#' @param family GLM family object
#' @param fixed_totals NULL or grouping vector for fixed totals
#' @param row_totals NULL or vector of row totals
#' @param no_dispersion Logical indicating if dispersion is fixed
#' @param nobs Number of observations
#' @param nvars Number of variables
#' @param keep Logical vector indicating which observations to keep
#' @param need_qr Logical indicating if QR decomposition is needed
#' @param need_hatvalues Logical indicating if hat values are needed
#'
#' @return A list of class "brglmFit_quantities" containing all computed quantities
compute_fit <- function(pars, y, x, weights, offset, family, 
                        fixed_totals = NULL, row_totals = NULL, 
                        no_dispersion = FALSE, nobs, nvars, keep,
                        need_qr = TRUE, need_hatvalues = TRUE,
                        need_inverse = TRUE) {
    
    # Extract Parameters
    betas <- pars[seq.int(nvars)]
    dispersion <- pars[nvars + 1]
    precision <- 1 / dispersion

    # Basic Quantities
    etas <- drop(x %*% betas + offset)
    mus <- family$linkinv(etas)
    mus_unscaled <- mus  # Keep original for gradient computation
    if (!is.null(fixed_totals)) {
        mus_totals <- as.vector(tapply(mus, fixed_totals, sum))[fixed_totals]
        mus <- mus * row_totals / mus_totals  # Scaled version
        etas <- family$linkfun(mus)  # Update etas to match scaled mus
    }

    # Mean Quantities
    d1mus <- family$mu.eta(etas)
    d2mus <- family$d2mu.deta(etas)
    varmus <- family$variance(mus)
    d1varmus <- family$d1variance(mus)
    working_weights <- weights * d1mus^2 / varmus

    # QR Decomposition and Hat Values (only if needed)
    qr_decomposition <- R_matrix <- Q_matrix <- hatvalues <- NULL
    info_beta <- inverse_info_beta <- NULL

    if (need_qr) {
        wx <- sqrt(working_weights) * x
        qr_decomposition <- qr(wx)
        R_matrix <- qr.R(qr_decomposition)
        
        # Information Matrices
        info_beta <- precision * crossprod(R_matrix)
        # chol2inv is O(p^3) — skip when solver does not need the explicit inverse
        inverse_info_beta <- if (need_inverse) dispersion * chol2inv(R_matrix) else NULL
        
        # Hat Values (only if needed and we have QR)
        if (need_hatvalues) {
            Q_matrix <- qr.Q(qr_decomposition) 
            hatvalues <- .rowSums(Q_matrix * Q_matrix, nobs, nvars, TRUE)
        }
    }

    # Dispersion Quantities (depends on family)
    if (!no_dispersion) {
        zetas <- -weights * precision
        
        # Derivatives of cumulant function (only for non-zero weights)
        d1afuns <- d2afuns <- d3afuns <- rep(NA_real_, nobs)
        d1afuns[keep] <- family$d1afun(zetas[keep])
        d2afuns[keep] <- family$d2afun(zetas[keep])
        d3afuns[keep] <- family$d3afun(zetas[keep])
        
        # Special case for Gamma family
        if (family$family == "Gamma") {
            d1afuns <- d1afuns - 2
        }
        
        # Deviance residuals
        deviance_residuals <- family$dev.resids(y, mus, weights)
        Edeviance_residuals <- weights * d1afuns
        
        # Information for dispersion parameter
        info_zeta <- 0.5 * sum(weights^2 * d2afuns, na.rm = TRUE) / dispersion^4
        inverse_info_zeta <- 1 / info_zeta
    } else {
        # Fixed dispersion (binomial, Poisson)
        zetas <- d1afuns <- d2afuns <- d3afuns <- NA_real_
        deviance_residuals <- family$dev.resids(y, mus, weights)
        Edeviance_residuals <- NA_real_
        info_zeta <- inverse_info_zeta <- NA_real_
    }

    # Gradient Computation
    # Use unscaled mus for gradient when fixed_totals is used
    mus_for_gradient <- if (!is.null(fixed_totals)) mus_unscaled else mus
    etas_for_gradient <- if (!is.null(fixed_totals)) family$linkfun(mus_unscaled) else etas
    d1mus_for_gradient <- family$mu.eta(etas_for_gradient)
    varmus_for_gradient <- family$variance(mus_for_gradient)
    
    score_components_beta <- weights * d1mus_for_gradient * (y - mus_for_gradient) / 
                            varmus_for_gradient * x
    grad_beta <- precision * .colSums(score_components_beta, nobs, nvars, TRUE)
    
    if (!no_dispersion) {
        grad_zeta <- 0.5 * precision^2 * 
                    sum(deviance_residuals - Edeviance_residuals, na.rm = TRUE)
    } else {
        grad_zeta <- NA_real_
    }

    # Return all computed quantities as class
    structure(list(
        # Parameters
        betas = betas,
        dispersion = dispersion,
        precision = precision,
        
        # Basic quantities
        etas = etas,
        mus = mus,
        mus_unscaled = mus_unscaled,  # For gradient with fixed_totals
        
        # Mean-related quantities
        d1mus = d1mus,
        d2mus = d2mus,
        varmus = varmus,
        d1varmus = d1varmus,
        working_weights = working_weights,
        
        # QR decomposition and related
        qr_decomposition = qr_decomposition,
        R_matrix = R_matrix,
        Q_matrix = Q_matrix, 
        hatvalues = hatvalues,
        
        # Information matrices (pre-computed)
        info_beta = info_beta,
        inverse_info_beta = inverse_info_beta,
        
        # Dispersion-related quantities
        zetas = zetas,
        d1afuns = d1afuns,
        d2afuns = d2afuns,
        d3afuns = d3afuns,
        deviance_residuals = deviance_residuals,
        Edeviance_residuals = Edeviance_residuals,
        info_zeta = info_zeta,
        inverse_info_zeta = inverse_info_zeta,
        
        # Gradients (pre-computed)
        grad_beta = grad_beta,
        grad_zeta = grad_zeta,
        score_components_beta = score_components_beta,
        
        # Metadata
        has_fixed_totals = !is.null(fixed_totals),
        no_dispersion = no_dispersion
    ), class = "brglmFit_quantities")
}

#' Compute mean bias-reducing adjustment
#'
#' @param pars Parameter vector
#' @param fit Object of class "brglmFit_quantities"
#' @param level 0 for beta parameters, 1 for dispersion
#' @param x Design matrix
#' @param nobs Number of observations
#' @param nvars Number of variables
#' @param weights Observation weights
#'
#' @return Adjustment vector
AS_mean_adjustment <- function(pars, fit, level = 0, 
                                       x, nobs, nvars, weights) {
    if (level == 0) {
        adj <- .colSums(0.5 * fit$hatvalues * fit$d2mus / fit$d1mus * x, 
                   nobs, nvars, TRUE)
        return(adj)
    } else {
        s1 <- sum(weights^3 * fit$d3afuns, na.rm = TRUE)
        s2 <- sum(weights^2 * fit$d2afuns, na.rm = TRUE)
        return((nvars - 2) / (2 * fit$dispersion) + 
            s1 / (2 * fit$dispersion^2 * s2))
    }
}

#' Compute Jeffreys prior penalty adjustment
#'
#' @param pars Parameter vector
#' @param fit Object of class "brglmFit_quantities"
#' @param level 0 for beta parameters, 1 for dispersion
#' @param x Design matrix
#' @param nobs Number of observations
#' @param nvars Number of variables
#' @param weights Observation weights
#' @param a Power of Jeffreys prior (default 0.5)
#'
#' @return Adjustment vector
AS_Jeffreys_adjustment <- function(pars, fit, level = 0, 
                                           x, nobs, nvars, weights, a = 0.5) {
    if (level == 0) {
        return(2 * a * .colSums(0.5 * fit$hatvalues * 
            (2 * fit$d2mus/fit$d1mus - fit$d1varmus * fit$d1mus / fit$varmus) * x, 
            nobs, nvars, TRUE))
    } else {
        s1 <- sum(weights^3 * fit$d3afuns, na.rm = TRUE)
        s2 <- sum(weights^2 * fit$d2afuns, na.rm = TRUE)
        return(2 * a * (-(nvars + 4)/(2 * fit$dispersion) + 
            s1/(2 * fit$dispersion^2 * s2)))
    }
}

#' Compute median bias-reducing adjustment
#'
#' @param pars Parameter vector
#' @param fit Object of class "brglmFit_quantities"
#' @param level 0 for beta parameters, 1 for dispersion
#' @param x Design matrix
#' @param nobs Number of observations
#' @param nvars Number of variables
#' @param weights Observation weights
#'
#' @return Adjustment vector
AS_median_adjustment <- function(pars, fit, level = 0, 
                                         x, nobs, nvars, weights) {
    if (level == 0) {
        info_unscaled <- fit$info_beta / fit$precision
        inverse_info_unscaled <- fit$inverse_info_beta / fit$dispersion

        b_vector <- numeric(nvars)
        for (j in seq.int(nvars)) {
            inverse_info_unscaled_j <- inverse_info_unscaled[j, ]
            vcov_j <- tcrossprod(inverse_info_unscaled_j) / inverse_info_unscaled_j[j]
            hats_j <- .rowSums((x %*% vcov_j) * x, nobs, nvars, TRUE) * fit$working_weights
            b_vector[j] <- inverse_info_unscaled_j %*% .colSums(x * (hats_j * 
                (fit$d1mus * fit$d1varmus / (6 * fit$varmus) - 0.5 * fit$d2mus/fit$d1mus)), 
                nobs, nvars, TRUE)
        }
        return(.colSums(0.5 * fit$hatvalues * fit$d2mus / fit$d1mus * x, 
                   nobs, nvars, TRUE) + 
           info_unscaled %*% b_vector)
    } else {
        s1 <- sum(weights^3 * fit$d3afuns, na.rm = TRUE)
        s2 <- sum(weights^2 * fit$d2afuns, na.rm = TRUE)
        return(nvars / (2 * fit$dispersion) + 
            s1 / (6 * fit$dispersion^2 * s2))
    }
}

#' Compute mixed bias-reducing adjustment
#'
#' @param pars Parameter vector
#' @param fit Object of class "brglmFit_quantities"
#' @param level 0 for beta parameters, 1 for dispersion
#' @param x Design matrix
#' @param nobs Number of observations
#' @param nvars Number of variables
#' @param weights Observation weights
#'
#' @return Adjustment vector
AS_mixed_adjustment <- function(pars, fit, level = 0, 
                                        x, nobs, nvars, weights) {
    if (level == 0) {
        return(.colSums(0.5 * fit$hatvalues * fit$d2mus/fit$d1mus * x, 
                   nobs, nvars, TRUE))
    } else {
        s1 <- sum(weights^3 * fit$d3afuns, na.rm = TRUE)
        s2 <- sum(weights^2 * fit$d2afuns, na.rm = TRUE)
        return(nvars / (2 * fit$dispersion) + 
            s1 / (6 * fit$dispersion^2 * s2))
    }
}

#' Estimate the ML of the dispersion parameter for gaussian, gamma and inverse Gaussian
#' Set the dispersion to 1 if Poisson or binomial
#' 
#' @param betas Regression coefficients
#' @param y Response vector
#' @param x Design matrix
#' @param weights Observation weights
#' @param offset Offset vector
#' @param family GLM family object
#' @param fixed_totals NULL or grouping vector for fixed totals
#' @param row_totals NULL or vector of row totals
#' @param no_dispersion Logical indicating if dispersion is fixed
#' @param nobs Number of observations
#' @param nvars Number of variables
#' @param keep Logical vector indicating which observations to keep
#' @param df_residual Residual degrees of freedom
#' @param control Control parameters
#' 
#' @seealso compute_fit
#'
#' @return List with dispersion and dispersion_ML
estimate_dispersion <- function(betas, y, x, weights, offset, family,
                               fixed_totals, row_totals, no_dispersion,
                               nobs, nvars, keep, df_residual, control) {
    if (no_dispersion) {
        disp <- 1
        dispML <- 1
    } else {
        if (df_residual > 0) {
            dispFit <- try(uniroot(f = function(phi) {
                theta <- c(betas, phi)
                cfit <- compute_fit(pars = theta, 
                                            y = y, 
                                            x = x,
                                            weights = weights,
                                            offset = offset,
                                            family = family,
                                            fixed_totals = fixed_totals,
                                            row_totals = row_totals,
                                            no_dispersion = no_dispersion,
                                            nobs = nobs,
                                            nvars = nvars,
                                            keep = keep,
                                            need_qr = FALSE,
                                            need_hatvalues = FALSE)
                cfit$grad_zeta
            }, lower = .Machine$double.eps, upper = 10000, tol = control$epsilon), silent = FALSE)
            if (inherits(dispFit, "try-error")) {
                warning("the ML estimate of the dispersion could not be calculated. An alternative estimate had been used as starting value.")
                dispML <- NA_real_
                disp <- NA_real_
            } else {
                disp <- dispML <- dispFit$root
            }
        } else { ## if the model is saturated dispML is NA_real_
            disp <- 1 ## A convenient value
            dispML <- NA_real_
        }
    }
    list(dispersion = disp, dispersion_ML = dispML)
}

#' Utility function for symbolic differentiation
#'
#' @param expr Expression to differentiate
#' @param name Variable name to differentiate with respect to
#' @param order Order of derivative
#'
#' @return Differentiated expression
DD <- function(expr, name, order = 1) {
    if(order < 1) stop("'order' must be >= 1")
    if(order == 1) D(expr, name)
    else DD(D(expr, name), name, order - 1)
}