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

#' Row-scale a sparse n x p matrix by a length-n weight vector
#'
#' Matrix's S4 `*` does NOT broadcast a length-n vector row-wise across an
#' n x p sparseMatrix (it errors or recycles incorrectly).  The standard fix,
#' `Diagonal(x=w) %*% x`, allocates a full n×n sparse diagonal matrix which is
#' expensive for large n due to S4 dispatch and memory overhead.
#'
#' Operate directly on the sparse column-format (CsparseMatrix)
#' x-slot.  In a CsparseMatrix, x@i gives the 0-based row index of every stored value,
#' so `x@x * w[x@i + 1L]` scales each stored value by the weight of its row - O(nnz),
#' no temporaries, no S4 dispatch overhead.
#' @keywords internal
.sparse_row_scale <- function(x, w) {
    # Coerce to CsparseMatrix once (no-op if already dgCMatrix / CsparseMatrix)
    cx <- methods::as(x, "CsparseMatrix")
    cx@x <- cx@x * w[cx@i + 1L]   # 0-based row indices → 1-based via +1L
    cx
}

#' Sparse-safe crossprod: x^T diag(w) x - returns a p x p dense matrix.
#'
#' Uses .sparse_row_scale so no n x n temporary is ever allocated.
#' Result is always a plain dense matrix (cheap when p << n).
#' @keywords internal
.sparse_wtd_crossprod <- function(x, w) {
    if (inherits(x, "sparseMatrix")) {
        # sqrt-scale rows, then standard crossprod: t(Wx) %*% Wx = X'WX
        # (same numeric result as Diagonal(w) %*% X path but ~10x less overhead)
        Wx <- .sparse_row_scale(x, sqrt(w))
        as.matrix(Matrix::crossprod(Wx))
    } else {
        crossprod(sqrt(w) * x)
    }
}

#' Sparse-safe diagonal of X^T diag(w) X (column squared norms).
#'
#' Reuses .sparse_row_scale to avoid the n x n Diagonal allocation.
#' @keywords internal
.sparse_wtd_colnorms2 <- function(x, w) {
    if (inherits(x, "sparseMatrix")) {
        Wx <- .sparse_row_scale(x, sqrt(w))
        drop(Matrix::colSums(Wx * Wx))   # element-wise ^2 then colSums
    } else {
        colSums(x^2 * w)
    }
}

#' Sparse-safe weighted column sums: x^T v.
#'
#' Replaces `colSums(v * x)` / `t(x) %*% v`.  Safe for both sparse and dense x.
#' @keywords internal
.sparse_crossprod_vec <- function(x, v) {
    as.vector(Matrix::crossprod(x, v))
}

#' Compute all quantities needed for a GLM fit
#'
#' @param pars Vector of parameters (betas and dispersion)
#' @param y Response vector
#' @param x Design matrix (may be a sparseMatrix)
#' @param weights Observation weights
#' @param offset Offset vector
#' @param family GLM family object
#' @param fixed_totals NULL or grouping vector for fixed totals
#' @param row_totals NULL or vector of row totals
#' @param no_dispersion Logical indicating if dispersion is fixed
#' @param nobs Number of observations
#' @param nvars Number of variables
#' @param keep Logical vector indicating which observations to keep
#' @param need_qr Logical indicating if QR/Cholesky decomposition is needed
#' @param need_hatvalues Logical indicating if hat values are needed (ignored when
#'   hatvalues_precomputed is non-NULL)
#' @param hatvalues_precomputed Optional numeric vector of pre-computed hat values.
#'   When non-NULL the expensive hat-value computation is skipped entirely and these
#'   values are injected directly into the returned object.  Used by the trust-region
#'   loop to pass a cached copy into candidate evaluations.
#' @param need_inverse Logical; if FALSE the O(p^3) chol2inv() / solve() call
#'   for inverse_info_beta is skipped (leaving it NULL).  Safe to set
#'   FALSE for families with fixed dispersion when AS_median is not used.
#' @param use_sparse Logical; if TRUE use sparse-optimised paths (no densification).
#'
#' @return A list of class "brglmFit_quantities" containing all computed quantities
compute_fit <- function(pars, y, x, weights, offset, family, 
                        fixed_totals = NULL, row_totals = NULL, 
                        no_dispersion = FALSE, nobs, nvars, keep,
                        need_qr = TRUE, need_hatvalues = TRUE,
                        hatvalues_precomputed = NULL,
                        need_inverse = TRUE,
                        use_sparse = FALSE) {
    
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

    # QR / Cholesky Decomposition and Hat Values (only if needed)
    qr_decomposition <- R_matrix <- Q_matrix <- hatvalues <- NULL
    info_beta <- inverse_info_beta <- NULL

    if (need_qr) {
        if (use_sparse) {
            # ---------------------------------------------------------------
            # Sparse path: avoid ever forming a dense n x p matrix until the
            # very end (hat values) when it's unavoidable.
            #
            # Core identity:  X'WX = (sqrt(W)X)^T (sqrt(W)X)
            # Form sqrt_w_x once (sparse, O(nnz)) and reuse for:
            #   (i) info_beta  = precision * crossprod(sqrt_w_x)       [p×p dense]
            #   (ii) R_chol     = chol(X'WX)                           [p×p dense]
            #   (iii) hat values = w_i ||R_chol^{-T} x_i||^2           [n dense vec]
            # ---------------------------------------------------------------
            sqrt_w_x  <- .sparse_row_scale(x, sqrt(working_weights))  # sparse n×p
            XtWX_dense <- as.matrix(Matrix::crossprod(sqrt_w_x))       # p×p dense

            info_beta <- precision * XtWX_dense

            # Always attempt Cholesky when need_qr=TRUE: R_matrix is needed for
            # hat values AND the periodic hat-refresh in the trust-region loop.
            R_chol <- tryCatch(chol(XtWX_dense), error = function(e) NULL)
            if (!is.null(R_chol)) {
                R_matrix <- R_chol
                if (need_inverse) {
                    inverse_info_beta <- dispersion * chol2inv(R_chol)
                }
            } else {
                # Cholesky failed (rank-deficient / numerically singular)
                if (need_inverse) {
                    inverse_info_beta <- tryCatch(dispersion * solve(XtWX_dense),
                                                  error = function(e) NULL)
                }
            }

            # Hat values: h_i = w_i * x_i^T (X'WX)^{-1} x_i
            #                 = ||R_chol^{-T} x_i||^2  * w_i
            #
            # Since R_chol is upper-triangular with R_chol^T R_chol = X'WX,
            # we need  z_i = R_chol^{-T} x_i  via a *lower*-triangular solve:
            #   forwardsolve(t(R_chol), x_i) <=> backsolve(R_chol, x_i, transpose=TRUE)
            #
            # Equivalently: the rows of sqrt_w_x are already sqrt(w_i)*x_i, and
            #   QR of sqrt_w_x gives Q with h_i = rowSums(Q^2).
            # But to avoid forming sqrt_w_x dense we use the R_chol route:
            #   R_chol^T R_chol = X'WX  ->  R_chol^{-T} x_i  ->  ||.||^2 * w_i
            if (is.null(hatvalues_precomputed) && need_hatvalues) {
                if (!is.null(R_matrix)) {
                    # backsolve with transpose=TRUE solves R_chol^T z = x_i for each col
                    # t(as.matrix(x)) is p×n; result is p×n; t(.) gives n×p
                    Xd <- as.matrix(x)                                    # n×p dense - unavoidable
                    Z  <- backsolve(R_matrix, t(Xd), transpose = TRUE)    # p×n: R^{-T} X^T
                    hatvalues <- working_weights * colSums(Z * Z)               # n-vec: w_i ||z_i||^2
                    Xd <- Z <- NULL
                } else if (!is.null(inverse_info_beta)) {
                    Xd <- as.matrix(x)
                    XV <- Xd %*% inverse_info_beta
                    hatvalues <- working_weights * rowSums(XV * Xd)
                    Xd <- XV <- NULL
                }
            }
            if (!is.null(hatvalues_precomputed)) hatvalues <- hatvalues_precomputed
            sqrt_w_x <- NULL   # free sparse n×p immediately

            qr_decomposition <- list(R = R_matrix, sparse_chol = TRUE)

        } else {
            # ---------------------------------------------------------------
            # Dense path (unchanged)
            # ---------------------------------------------------------------
            wx <- sqrt(working_weights) * x
            qr_decomposition <- qr(wx)
            R_matrix <- qr.R(qr_decomposition)
            
            # Information Matrices
            info_beta <- precision * crossprod(R_matrix)
            # chol2inv is O(p^3) - skip when solver does not need the explicit inverse
            inverse_info_beta <- if (need_inverse) dispersion * chol2inv(R_matrix) else NULL
            
            # qr.Q() call.  Otherwise compute from scratch only if requested.
            if (is.null(hatvalues_precomputed) && need_hatvalues) {
                Q_matrix  <- qr.Q(qr_decomposition)
                Q_matrix  <- as.matrix(Q_matrix)
                hatvalues <- rowSums(Q_matrix * Q_matrix)
            }
            if (!is.null(hatvalues_precomputed)) hatvalues <- hatvalues_precomputed
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

    # Sparse-safe gradient: use crossprod instead of colSums(scalar * x)
    # For dense x, .sparse_crossprod_vec falls through to as.vector(crossprod(x, v))
    score_weights <- weights * d1mus_for_gradient * (y - mus_for_gradient) / varmus_for_gradient
    grad_beta <- precision * .sparse_crossprod_vec(x, score_weights)
    # Keep score_components_beta for backward compat (only materialised if used downstream)
    score_components_beta <- NULL  # computed lazily to avoid densification

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
        # Cached score weights for adjustment functions (avoids re-densification)
        score_weights = score_weights,
        
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
        # Sparse-safe: crossprod(x, v) instead of colSums(v * x)
        v <- 0.5 * fit$hatvalues * fit$d2mus / fit$d1mus
        adj <- .sparse_crossprod_vec(x, v)
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
        v <- 2 * a * 0.5 * fit$hatvalues *
             (2 * fit$d2mus/fit$d1mus - fit$d1varmus * fit$d1mus / fit$varmus)
        return(.sparse_crossprod_vec(x, v))
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
 
        v_hat <- 0.5 * fit$hatvalues * fit$d2mus / fit$d1mus
 
        if (inherits(x, "sparseMatrix")) {
            # Sparse path: materialise XV = X %*% V once (n x p dense),
            # then reuse it for every j — mirrors AS_median_adjustment_new.
            XV   <- as.matrix(x %*% inverse_info_unscaled)  # n x p dense
            d_V  <- diag(inverse_info_unscaled)              # length-p diagonal
 
            c_vec    <- fit$d1mus * fit$d1varmus / (6 * fit$varmus) - 0.5 * fit$d2mus / fit$d1mus
            wc       <- fit$working_weights * c_vec
            b_vector <- colSums(XV^3 * wc) / d_V            # O(np), no loop
 
            # sparse-safe: use Matrix::crossprod instead of base .colSums
            return(.sparse_crossprod_vec(x, v_hat) + as.vector(info_unscaled %*% b_vector))
        } else {
            # Dense path: original per-column loop, unchanged
            b_vector <- numeric(nvars)
            for (j in seq.int(nvars)) {
                inverse_info_unscaled_j <- inverse_info_unscaled[j, ]
                vcov_j <- tcrossprod(inverse_info_unscaled_j) / inverse_info_unscaled_j[j]
                hats_j <- .rowSums((x %*% vcov_j) * x, nobs, nvars, TRUE) * fit$working_weights
                b_vector[j] <- inverse_info_unscaled_j %*% .colSums(x * (hats_j * 
                    (fit$d1mus * fit$d1varmus / (6 * fit$varmus) - 0.5 * fit$d2mus/fit$d1mus)), 
                    nobs, nvars, TRUE)
            }
 
            return(.colSums(v_hat * x, nobs, nvars, TRUE) + 
               as.vector(info_unscaled %*% b_vector))
        }
    } else {
        s1 <- sum(weights^3 * fit$d3afuns, na.rm = TRUE)
        s2 <- sum(weights^2 * fit$d2afuns, na.rm = TRUE)
        return(nvars / (2 * fit$dispersion) + 
            s1 / (6 * fit$dispersion^2 * s2))
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
AS_median_adjustment_new <- function(pars, fit, level = 0,
                                 x, nobs, nvars, weights) {
    if (level == 0) {
        info_unscaled         <- fit$info_beta / fit$precision
        inverse_info_unscaled <- fit$inverse_info_beta / fit$dispersion

        # XV = X V where V = inverse_info_unscaled (p x p)
        # For sparse x: XV is n x p — this is the unavoidable O(np^2) step for AS_median
        # but we keep x sparse through the multiply so only XV materialises as dense
        XV    <- as.matrix(x %*% inverse_info_unscaled)   # n x p dense
        d_V   <- diag(inverse_info_unscaled)               # length-p diagonal

        # Per-observation weight for the cubic term
        c_vec <- fit$d1mus * fit$d1varmus / (6 * fit$varmus) -
                 0.5 * fit$d2mus / fit$d1mus                # length n

        # b_vector[j] = colSums(XV^3 * wc) / d_V[j]
        wc    <- fit$working_weights * c_vec
        b_vector <- colSums(XV^3 * wc) / d_V               # O(np)

        # Final term: sparse-safe crossprod
        v_hat <- 0.5 * fit$hatvalues * fit$d2mus / fit$d1mus
        return(.sparse_crossprod_vec(x, v_hat) + as.vector(info_unscaled %*% b_vector))
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
        v <- 0.5 * fit$hatvalues * fit$d2mus / fit$d1mus
        return(.sparse_crossprod_vec(x, v))
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
                               nobs, nvars, keep, df_residual, control,
                               use_sparse = FALSE) {
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
                                            need_hatvalues = FALSE,
                                            use_sparse = use_sparse)
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