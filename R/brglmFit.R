# Copyright (C) 2026- Ioannis Kosmidis, Oliver Clark 

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

#' Trust Region Fitting Method for Bias-Reduced GLMs
#'
#' Implements trust region optimisation with CG-Steihaug subproblem solver
#' for bias-reduced generalised linear models. This method provides robust
#' global convergence and exploits matrix-free operations for scalability.
#'
#' @inheritParams brglmFit
#' 
#' @details
#' This implementation follows the trust region framework from Nocedal & Wright (2006)
#' with adaptations for bias-reduced GLMs. The CG-Steihaug method (Algorithm 7.2)
#' is used to solve the trust region subproblem without forming the Hessian explicitly.
#'
#' Key features:
#' - Matrix-free Hessian-vector products: O(np) per CG iteration vs O(np^2 + p^3)
#' - Handles semi-definite Hessians (detects negative curvature)
#' - Natural sparse matrix support
#' - Adaptive trust region radius
#'
#' Currently supports binomial family with fixed dispersion only.
#' 
#' @return
#'
#' An object inheriting from [`"brglmFit"`][brglmFit()] object, which
#' is a list having the same elements to the list that
#' [stats::glm.fit()] returns, with a few extra arguments.
#'
#' @references
#' Nocedal, J. and Wright, S.J. (2006). Numerical Optimization (2nd ed.). Springer.
#' Steihaug, T. (1983). The conjugate gradient method and trust regions in large
#' scale optimization. SIAM Journal on Numerical Analysis, 20(3), 626-637.
#' 
#' @examples
#' ## The lizards example from ?brglm::brglm
#' data("lizards", package = "brglm2")
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
#' # Maximum penalized likelihood with Jeffreys prior penatly
#' lizards_Jeffreys <- glm(cbind(grahami, opalinus) ~ height + diameter +
#'                         light + time, family = binomial(logit), data = lizards,
#'                         method = "brglmFit", type = "MPL_Jeffreys")
#' # lizards_Jeffreys is the same fit as lizardsBR_mean (see Firth, 1993)
#' all.equal(coef(lizardsBR_mean), coef(lizards_Jeffreys))
#'
#' # Maximum penalized likelihood with powers of the Jeffreys prior as
#' # penalty. See Kosmidis & Firth (2021) for the finiteness and
#' # shrinkage properties of the maximum penalized likelihood
#' # estimators in binomial response models
#' \donttest{
#' a <- seq(0, 20, 0.5)
#' coefs <- sapply(a, function(a) {
#'       out <- glm(cbind(grahami, opalinus) ~ height + diameter +
#'              light + time, family = binomial(logit), data = lizards,
#'              method = "brglmFit", type = "MPL_Jeffreys", a = a)
#'       coef(out)
#' })
#' # Illustration of shrinkage as a grows
#' matplot(a, t(coefs), type = "l", col = 1, lty = 1)
#' abline(0, 0, col = "grey")
#'}
#' 
#' ## endometrial data from Heinze & Schemper (2002) (see ?endometrial)
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
brglmFit <- function(x, y, weights = rep(1, nobs), 
                                    start = NULL, etastart = NULL,
                                    mustart = NULL, offset = rep(0, nobs), 
                                    family = gaussian(),
                                    control = list(), intercept = TRUE,
                                    fixed_totals = NULL, singular.ok = TRUE) {
    
    control <- do.call("brglmControl", control)
    
    adjustment_function <- switch(control$type,
                                    "correction" = AS_mean_adjustment,
                                    "AS_mean" = AS_mean_adjustment,
                                    "AS_median" = AS_median_adjustment,
                                    "AS_mixed" = AS_mixed_adjustment,
                                    "MPL_Jeffreys" = AS_Jeffreys_adjustment,
                                    "ML" = function(pars, ...) 0)
    
    is_ML <- control$type == "ML"
    no_dispersion <- family$family %in% c("poisson", "binomial")
    
    # Family compatibility check
    if (!no_dispersion) {
        warning("Trust region method currently only supports binomial and Poisson families (fixed dispersion)")
        # Return a minimal valid object so tests can continue
        # Need to set up minimal dimensions to avoid indexing errors
        x <- as.matrix(x)
        nvars <- ncol(x)
        nobs <- NROW(y)
        
        return(list(
            coefficients = structure(rep(NA_real_, nvars), .Names = colnames(x)),
            residuals = rep(NA_real_, nobs),
            fitted.values = rep(NA_real_, nobs),
            R = matrix(NA_real_, nvars, nvars),
            rank = 0,
            qr = list(qr = matrix(NA_real_, nobs, nvars), rank = 0, 
                     qraux = rep(NA_real_, nvars), pivot = seq_len(nvars), tol = 1e-7),
            family = family,
            linear.predictors = rep(NA_real_, nobs),
            deviance = NA,
            aic = NA,
            null.deviance = NA,
            iter = 0,
            weights = rep(0, nobs),
            prior.weights = weights,
            df.residual = 0,
            df.null = 0,
            y = y,
            converged = FALSE,
            boundary = FALSE,
            dispersion = NA,
            dispersion_ML = NA,
            transformed_dispersion = NA,
            info_transformed_dispersion = NA,
            grad = structure(rep(NA_real_, nvars + 1), 
                           .Names = c(colnames(x), "Transformed dispersion")),
            transformation = if (is.null(control$transformation)) "identity" else control$transformation,
            type = control$type,
            control = control,
            class = "brglmFit"
        ))
    }
    
    # Handle fixed_totals for Poisson
    if (is.null(fixed_totals)) {
        has_fixed_totals <- FALSE
        row_totals <- NULL
    } else {
        if (family$family == "poisson") {
            row_totals <- as.vector(tapply(y, fixed_totals, sum))[fixed_totals]
            has_fixed_totals <- TRUE
        } else {
            has_fixed_totals <- FALSE
            row_totals <- NULL
        }
    }
    
    # Matrix setup
    x <- as.matrix(x)
    betas_names <- dimnames(x)[[2L]]
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(betas_names) & !EMPTY) {
        betas_names <- colnames(x) <- paste0("x", seq.int(nvars))
    }
    ynames <- if (is.matrix(y)) rownames(y) else names(y)
    converged <- FALSE
    nobs <- NROW(y)
    if (is.null(weights)) {
        weights <- rep.int(1, nobs)
    }
    if (missing_offset <- is.null(offset)) {
        offset <- rep.int(0, nobs)
    }
    
    # Enrich family (needed for d2mu.deta, etc.)
    family <- enrichwith::enrich(family, with = c("d1afun", "d2afun", "d3afun", "d1variance"))
    ok_links <- c("logit", "probit", "cauchit", "cloglog", "identity", "log", "sqrt", "inverse")

    if (isTRUE(family$family %in% c("quasi", "quasibinomial", "quasipoisson"))) {
        stop("`brglmFit` does not currently support the `quasi`, `quasipoisson` and `quasibinomial` families.")
    }

    if ((family$link %in% ok_links) | (grepl("mu\\^", family$link))) {
        linkglm <- make.link(family$link)
        linkglm <- enrichwith::enrich(linkglm, with = "d2mu.deta")
        family[names(linkglm)] <- linkglm
    }
    if (is.null(family$d2mu.deta)) {
        family$d2mu.deta <- function(eta) numDeriv::grad(family$mu.eta, eta)
    }
    
    # Extract family functions
    variance <- family$variance
    linkinv <- family$linkinv
    linkfun <- family$linkfun
    mu.eta <- family$mu.eta
    dev.resids <- family$dev.resids
    valid_eta <- unless_null(family$valideta, function(eta) TRUE)
    valid_mu <- unless_null(family$validmu, function(mu) TRUE)
    
    # Initialize
    eval(family$initialize)
    
    # Handle empty model
    if (EMPTY) {
        stop("Trust region method requires at least one predictor")
    }
    
    # Detect aliasing
    if (!isTRUE(control$check_aliasing)) {
        is_full_rank <- TRUE
        rank <- nvars_all <- nvars
        betas_names_all <- betas_names
    } else {
        qrx <- qr(x)
        rank <- qrx$rank
        is_full_rank <- rank == nvars
        if (!isTRUE(singular.ok) && !isTRUE(is_full_rank)) {
            stop("singular fit encountered")
        }
        if (!isTRUE(is_full_rank)) {
            aliased <- qrx$pivot[seq.int(qrx$rank + 1, nvars)]
            X_all <- x
            x <- x[, -aliased]
            nvars_all <- nvars
            nvars <- ncol(x)
            betas_names_all <- betas_names
            betas_names <- betas_names[-aliased]
        } else {
            nvars_all <- nvars
            betas_names_all <- betas_names
        }
    }
    betas_all <- structure(rep(NA_real_, nvars_all), .Names = betas_names_all)
    keep <- weights > 0
    nkeep <- sum(keep)
    df_residual <- nkeep - rank
    
    # ===== STARTING VALUES =====
    
    if (is.null(start)) {
        # Use ML fit with adjusted response (same as brglmFit)
        adj <- control$response_adjustment
        if (is.null(adj)) adj <- nvars/nobs
        
        if (family$family == "binomial") {
            weights.adj <- weights + adj
            y.adj <- (weights * y + 0.5 * adj) / weights.adj
        } else {
            weights.adj <- weights
            y.adj <- y + if (family$family == "poisson") 0.5 * adj else 0
        }
        
        suppressWarnings(
            tempFit <- glm.fit(x = x, y = y.adj, weights = weights.adj,
                                etastart = etastart, mustart = mustart,
                                offset = offset, family = family,
                                control = list(epsilon = control$epsilon,
                                                maxit = 10000, trace = FALSE),
                                intercept = intercept)
        )
        betas <- coef(tempFit)
        names(betas) <- betas_names
    } else {
        if (length(start) == nvars_all) {
            betas_all <- start
            names(betas_all) <- betas_names_all
            if (!isTRUE(is_full_rank)) {
                betas_all[aliased] <- NA_real_
                betas <- betas_all[-aliased]
            } else {
                betas <- betas_all
            }
        } else {
            stop("length of 'start' should equal number of parameters")
        }
    }
    
    # Fixed dispersion
    dispersion <- 1
    dispersion_ML <- 1
    
    # ===== TRUST REGION PARAMETERS =====
    # Following Nocedal & Wright defaults 
    Delta <- 1.0                    # Initial trust region radius
    Delta_max <- 1e10               # Maximum radius
    eta_shrink <- 0.25              # Shrink if rho < eta_shrink
    eta_expand <- 0.75              # Expand if rho > eta_expand
    shrink_factor <- 0.25           # Multiply Delta by this when shrinking
    expand_factor <- 2.0            # Multiply Delta by this when expanding
    
    # CG-Steihaug parameters
    cg_tol <- 0.1                   # Relative CG tolerance (More appropriate number?)
    cg_maxiter <- min(nvars, 50)    # Maximum CG iterations (More appropriate number?)
    
    # Hat value recomputation frequency
    hat_recompute_freq <- 5         # Recompute exact hat values every N iterations
    
    # ===== TRUST REGION ITERATION =====
    boundary <- FALSE
    
    #browser()
    # Initial evaluation
    obj <- compute_objective(betas = betas, y = y, x = x, weights = weights,
                            offset = offset, family = family,
                            adjustment_function = adjustment_function,
                            fixed_totals = fixed_totals, row_totals = row_totals,
                            no_dispersion = no_dispersion, nobs = nobs, nvars = nvars,
                            keep = keep)
    f_current <- obj$value
    grad <- obj$gradient 
    J <- obj$jacobian
    H <- obj$hessian
    r <- obj$residual

    if (control$trace) {
        cat("Trust Region Method - Merit Function Formulation\n")
        cat("Minimizing f(beta) = ||r(beta)||²/2 where r = adjusted_score\n")
        cat("Initial f:", f_current, "\n")
        cat("Initial ||r||:", sqrt(2*f_current), "\n")
        cat("Initial ||Nabla f||:", sqrt(sum(grad^2)), "\n\n")
    }
    
    # Main trust region loop
    for (iter in seq_len(control$maxit)) {
        
        # Solve trust region subproblem: min_p { g'p + (1/2)p'Hp } s.t. ||p|| <= Delta
        # Standard formulation - grad is ∇f (gradient of objective)
        step <- cg_steihaug_subproblem_with_H(
            grad = grad,  # Pass gradient of f 
            H = H,
            Delta = Delta,
            tol = cg_tol,
            maxiter = cg_maxiter
        )
        
        # Predicted reduction (standard quadratic model)
        # pred = m(0) - m(p) = -(g'p + (1/2)p'Hp)
        Bp <- drop(H %*% step$p)
        pred_reduction <- -(sum(grad * step$p) + 0.5 * sum(step$p * Bp))
        #Jp <- drop(J %*% step$p)
        #pred_reduction <- -(2*f_current - (sum((r + Jp)^2)))

        # Try the step
        betas_new <- betas + step$p
        
        #browser()
        # Evaluate at new point
        obj_new <- compute_objective(betas = betas_new, y = y, x = x, weights = weights,
                              offset = offset, family = family,
                              adjustment_function = adjustment_function,
                              fixed_totals = fixed_totals, row_totals = row_totals,
                              no_dispersion = no_dispersion, nobs = nobs, nvars = nvars,
                              keep = keep)
        f_new <- obj_new$value
        grad_new <- obj_new$gradient
        H_new <- obj_new$hessian
        J_new <- obj_new$jacobian
        r_new <- obj_new$residual
        
        # Actual reduction
        actual_reduction <- (f_current - f_new)

        #browser()
        
        # Reduction ratio
        rho <- if (abs(pred_reduction) < 1e-16) {
            if (actual_reduction > 0) 1.0 else 0.0
        } else {
            actual_reduction / pred_reduction
        }
        
        # Update trust region radius
        if (rho < eta_shrink) {
            Delta <- shrink_factor * Delta 
        } else if (rho > eta_expand) {
            Delta <- min(expand_factor * Delta, Delta_max)
        }

        # Accept or reject step
        if (rho > eta_shrink) { # Try smaller factor
            # Accept step
            betas <- betas_new
            f_current <- f_new
            grad <- grad_new
            H <- H_new
            J <- J_new
            r <- r_new
            accept_status <- "ACCEPT"
        } else {
            # Reject step
            accept_status <- "REJECT"
        }
        
        # Check for minimum radius
        if (Delta < 1e-16) {
            warning("Trust region radius became too small")
            break
        }
        
        # Trace output
        if (control$trace) {
            grad_norm <- sqrt(sum(grad^2))
            step_norm <- sqrt(sum(step$p^2))
            cat(sprintf("Iter %3d: f = %.6e, ||grad|| = %.6e, ||step|| = %.6e, Delta = %.3e, rho = %.3f, %s%s\n",
                       iter, f_current, grad_norm, step_norm, Delta, rho,
                       if (step$on_boundary) "[BOUND] " else "",
                       accept_status))
        }
        
        # Convergence check on GRADIENT of objective
        # We've converged when ∇f ≈ 0, which implies r ≈ 0
        grad_norm <- sqrt(sum(grad^2))
        if (grad_norm < control$epsilon) {
            converged <- TRUE
            break
        }
        
        # Alternative: converge on residual norm
        # r_norm <- sqrt(sum(r_current^2))
        # if (r_norm < control$epsilon) {
        #     converged <- TRUE
        #     break
        # }
    }
    
    # Final convergence check
    if (!converged && iter >= control$maxit) {
        warning("Trust region algorithm did not converge within maxit iterations")
    }
    
    # ===== FINAL EVALUATION AND OUTPUT FORMATTING =====
    
    # Restore full parameter vector if aliasing occurred
    if (!isTRUE(is_full_rank)) {
        x <- X_all
        betas_all[betas_names] <- betas
        betas <- betas_all
        betas[is.na(betas)] <- 0
        nvars <- nvars_all
    } else {
        betas_all <- betas
    }
    
    # Final fit with full QR decomposition
    theta <- c(betas, dispersion)
    fit <- compute_fit(pars = theta, y = y, x = x, weights = weights,
                      offset = offset, family = family,
                      fixed_totals = fixed_totals, row_totals = row_totals,
                      no_dispersion = no_dispersion, nobs = nobs, nvars = nvars,
                      keep = keep, need_qr = TRUE, need_hatvalues = TRUE)
    
    qr.Wx <- fit$qr_decomposition
    mus <- fit$mus
    etas <- fit$etas
    residuals <- (y - mus) / fit$d1mus
    working_weights <- fit$working_weights
    
    # Boundary check
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
        if (any(mus > 1 - eps) || any(mus < eps)) {
            warning("brglmFit: fitted probabilities numerically 0 or 1 occurred")
            boundary <- TRUE
        }
    }
    if (family$family == "poisson") {
        if (any(mus < eps)) {
            warning("brglmFit: fitted rates numerically 0 occurred")
            boundary <- TRUE
        }
    }
    
    # Working weights and residuals
    wt <- rep.int(0, nobs)
    wt[keep] <- working_weights[keep]
    names(wt) <- names(residuals) <- names(mus) <- names(etas) <- ynames
    names(weights) <- names(y) <- ynames
    
    # Null deviance (same logic as brglmFit)
    control0 <- control
    control0$maxit <- 1000
    if (intercept & missing_offset) {
        nullFit <- brglmFit(
            x = x[, "(Intercept)", drop = FALSE], y = y, weights = weights,
            offset = rep(0, nobs), family = family, intercept = TRUE,
            control = control0[c("epsilon", "maxit", "type", "transformation", "slowit")],
            start = linkfun(mean(y))
        )
        nullmus <- nullFit$fitted.values
    }
    if (!intercept) {
        nullmus <- linkinv(offset)
    }
    if (intercept & !missing_offset) {
        nullmus <- mus
    }
    
    nulldev <- sum(dev.resids(y, nullmus, weights))
    nulldf <- nkeep - as.integer(intercept)
    deviance <- sum(dev.resids(y, mus, weights))
    aic.model <- family$aic(y, nobs, mus, weights, deviance) + 2 * rank
    
    # Gradient for output (residual at final point for compatibility)
    adjustment <- adjustment_function(theta, fit = fit, level = 0,
                                     x = x, nobs = nobs, nvars = nvars, weights = weights)
    r_final <- fit$grad_beta + adjustment
    
    adjusted_grad_all <- rep(NA_real_, nvars_all + 1)
    names(adjusted_grad_all) <- c(betas_names_all, "Transformed dispersion")
    adjusted_grad_all[betas_names] <- r_final  # Store residual for compatibility
    adjusted_grad_all["Transformed dispersion"] <- NA_real_
    
    # Return in brglmFit format
    list(
        coefficients = betas_all,
        residuals = residuals,
        fitted.values = mus,
        R = if (!EMPTY) qr.R(qr.Wx),
        rank = rank,
        qr = if (!EMPTY) structure(qr.Wx[c("qr", "rank", "qraux", "pivot", "tol")], 
                                   class = "qr"),
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
        transformed_dispersion = 1,
        info_transformed_dispersion = NA_real_,
        grad = adjusted_grad_all,
        transformation = control$transformation,
        type = control$type,
        control = control,
        class = "brglmFit"
    )
}

# ===== HELPER FUNCTIONS FOR TRUST REGION METHOD =====

#' Matrix-free Hessian-vector product for GLMs
#'
#' Computes (X'WX)v without forming the matrix explicitly
#' Cost: O(np) instead of O(np^2 + p^3)
#'
#' @param v Vector to multiply
#' @param x Design matrix (n x p)
#' @param weights Working weights (length n)
#' @return (X'WX)v as vector of length p
hessian_vector_product <- function(v, x, weights) {
    # (X'WX)v = X'(W(Xv))
    Xv <- drop(x %*% v)          # O(np)
    WXv <- weights * Xv          # O(n)
    XtWXv <- drop(crossprod(x, WXv))  # O(np)
    return(XtWXv)
}

#' CG-Steihaug Trust Region Subproblem Solver
#'
#' Solves: min_p { g'p + 0.5 p'Bp } subject to ||p|| <= Delta
#' using conjugate gradient with early termination (Steihaug, 1983)
#'
#' @param neg_grad Negative gradient vector (length p)
#' @param x Design matrix (n x p)
#' @param weights_sqrt Square root of working weights (length n)
#' @param Delta Trust region radius
#' @param tol Convergence tolerance for CG
#' @param maxiter Maximum CG iterations
#'
#' @return List with:
#'   - p: Step vector
#'   - on_boundary: Logical indicating if step hit trust region boundary
#'   - cg_iter: Number of CG iterations performed
#'
#' @references
#' Steihaug, T. (1983). The conjugate gradient method and trust regions.
#' SIAM J. Numer. Anal., 20(3), 626-637.
#' Nocedal & Wright (2006), Algorithm 7.2
cg_steihaug_subproblem <- function(neg_grad, x, weights, Delta, 
                                   tol = 0.1, maxiter = 50) {
    
    p <- length(neg_grad)
    
    # Initialize
    z <- numeric(p)              # Current iterate
    r <- neg_grad                # Residual (initially Bz = 0)
    d <- -r                      # Search direction
    
    on_boundary <- FALSE
    
    for (j in seq_len(maxiter)) {
        
        # Compute Bd using matrix-free product
        # B = X'WX, so we need working_weights
        Bd <- hessian_vector_product(d, x, weights)
        
        # Curvature
        dBd <- sum(d * Bd)
        
        # Check for non-positive curvature
        if (dBd <= 0) {
            # Find tau such that ||z + tau*d|| = Delta
            tau <- find_boundary_step(z, d, Delta)
            z <- z + tau * d
            on_boundary <- TRUE
            break
        }
        
        # CG step size
        r_norm_sq <- sum(r * r)
        alpha <- r_norm_sq / dBd
        
        # Check if step would exit trust region
        z_new <- z + alpha * d
        if (sqrt(sum(z_new * z_new)) > Delta) {
            # Find tau such that ||z + tau*d|| = Delta
            tau <- find_boundary_step(z, d, Delta)
            z <- z + tau * d
            on_boundary <- TRUE
            break
        }
        
        # Accept CG step
        z <- z_new
        
        # Update residual
        r <- r + alpha * Bd
        r_norm_sq_new <- sum(r * r)
        
        # Check convergence with relative tolerance (tol% of initial grad norm)
        if (sqrt(r_norm_sq_new) < tol * sqrt(sum(neg_grad * neg_grad))) { 
            break
        }
        
        # Compute new search direction
        beta <- r_norm_sq_new / r_norm_sq
        d <- -r + beta * d
    }
    
    list(
        p = z,
        on_boundary = on_boundary,
        cg_iter = j
    )
}

#' Find step to trust region boundary
#'
#' Finds tau >= 0 such that ||z + tau*d|| = Delta
#' Solves: ||z + tau*d||^2 = Delta^2
#'
#' @param z Current point
#' @param d Direction
#' @param Delta Trust region radius
#' @return Scalar tau
find_boundary_step <- function(z, d, Delta) {
    # ||z + tau*d||^2 = Delta^2
    # (z + tau*d)'(z + tau*d) = Delta^2
    # z'z + 2*tau*z'd + tau^2*d'd = Delta^2
    # tau^2*(d'd) + tau*(2*z'd) + (z'z - Delta^2) = 0
    
    a <- sum(d * d)
    b <- 2 * sum(z * d)
    c <- sum(z * z) - Delta^2
    
    # Quadratic formula: tau = (-b + sqrt(b^2 - 4ac)) / (2a)
    # We want the positive root
    discriminant <- b^2 - 4 * a * c
    
    if (discriminant < 0) {
        # Numerical issue - should not happen in theory
        return(0)
    }
    
    tau1 <- (-b + sqrt(discriminant)) / (2 * a)
    tau2 <- (-b - sqrt(discriminant)) / (2 * a)

    # Choose smallest positive root (first boundary intersection)
    if (tau1 > 0 && tau2 > 0) {
        return(min(tau1, tau2))
    } else if (tau1 > 0) {
        return(tau1)
    } else if (tau2 > 0) {
        return(tau2)
    } else {
        return(0)  # Fallback
    }
}

# ===== S3 METHODS FOR brglmFit  =====

#' Extract model coefficients from [`"brglmFit"`][brglmFit] objects
#'
#' @inheritParams stats::coef
#' @param model one of `"mean"` (default), `"dispersion"`, `"full",
#'     to return the estimates of the parameters in the linear
#'     prediction only, the estimate of the dispersion parameter only,
#'     or both, respectively.
#'
#' @details
#'
#' See [coef()] for more details.
#'
#' @seealso
#'
#' [coef()]
#'
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

#' [summary()] method for [`"brglmFit"`][brglmFit()] objects
#'
#' @inheritParams stats::summary.glm
#'
#' @details The interface of the summary method for [`"brglmFit"`][brglmFit]
#'     objects is identical to that of [`"glm"`][glm] objects. The summary
#'     method for [`"brglmFit"`][brglmFit] objects computes the p-values of the
#'     individual Wald statistics based on the standard normal
#'     distribution, unless the family is Gaussian, in which case a t
#'     distribution with appropriate degrees of freedom is used.
#'
#' @seealso [summary.glm()] and [glm()]
#'
#' @examples
#' ## For examples see `examples(brglmFit)`
#'
#' @method summary brglmFit
#' @export
summary.brglmFit <- function(object, dispersion = NULL,
                             correlation = FALSE, symbolic.cor = FALSE,
                             ...) {
    if (is.null(dispersion)) {
        if (object$family$family == "Gaussian") {
            dispersion <- NULL
        } else {
            dispersion <- object$dispersion
        }
    }
    out <- summary.glm(object, dispersion = dispersion,
                       correlation = correlation,
                       symbolic.cor = symbolic.cor, ...)
    out$type <- object$type
    class(out) <- c("summary.brglmFit", class(out))
    out
}

#' Method for computing confidence intervals for one or more
#' regression parameters in a [`"brglmFit"`][brglmFit()] object
#'
#' @inheritParams stats::confint
#'
#' @method confint brglmFit
#' @export
confint.brglmFit <- function(object, parm, level = 0.95, ...) {
    confint.default(object, parm, level, ...)
}

#' Return the variance-covariance matrix for the regression parameters
#' in a [brglmFit()] object
#'
#' @inheritParams stats::vcov.glm
#' @param model character specifying for which component of the model coefficients should be extracted.
#'
#' @details
#'
#' The options for `model` are `"mean"` for mean regression parameters
#' only (default), `"dispersion"` for the dispersion parameter (or the
#' transformed dispersion; see [brglm_control()]), and `"full"` for
#' both the mean regression and the (transformed) dispersion
#' parameters.
#'
#' @method vcov brglmFit
#' @export
vcov.brglmFit <- function(object, model = c("mean", "full", "dispersion"), complete = TRUE, ...) {
    model <- match.arg(model)
    switch(model,
           mean = {
        vcov(summary.brglmFit(object, ...), complete = complete)
    },
    dispersion = {
        vtd <- 1/object$info_transformed_dispersion
        ntd <- paste0(object$transformation, "(dispersion)")
        names(vtd) <- ntd
        vtd
    },
    full = {
        vbetas <- vcov(summary.brglmFit(object, ...), complete = complete)
        vtd <- 1/object$info_transformed_dispersion
        nBetasAll <- c(rownames(vbetas), paste0(object$transformation, "(dispersion)"))
        vBetasAll <- cbind(rbind(vbetas, 0),
                           c(numeric(nrow(vbetas)), vtd))
        dimnames(vBetasAll) <- list(nBetasAll, nBetasAll)
        vBetasAll
    })
}

## Almost all code in print.summary.brglmFit is from
## stats:::print.summary.glm apart from minor modifications
#' @rdname summary.brglmFit
#' @method print summary.brglmFit
#' @export
print.summary.brglmFit <- function (x, digits = max(3L, getOption("digits") - 3L),
                                    symbolic.cor = x$symbolic.cor,
                                    signif.stars = getOption("show.signif.stars"), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    cat("Deviance Residuals: \n")
    if (x$df.residual > 5) {
        x$deviance.resid <- setNames(quantile(x$deviance.resid,
                                              na.rm = TRUE), c("Min", "1Q", "Median", "3Q", "Max"))
    }
    xx <- zapsmall(x$deviance.resid, digits + 1L)
    print.default(xx, digits = digits, na.print = "", print.gap = 2L)
    if (length(x$aliased) == 0L) {
        cat("\nNo Coefficients\n")
    } else {
        df <- if ("df" %in% names(x))
                  x[["df"]]
              else NULL
        if (!is.null(df) && (nsingular <- df[3L] - df[1L]))
            cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
                sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients
        if (!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn,
                                                                     colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                     na.print = "NA", ...)
    }
    cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ",
        format(x$dispersion), ")\n\n", apply(cbind(paste(format(c("Null",
                                                                  "Residual"), justify = "right"), "deviance:"), format(unlist(x[c("null.deviance",
                                                                                                                                   "deviance")]), digits = max(5L, digits + 1L)), " on",
                                                   format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"),
                                             1L, paste, collapse = " "), sep = "")
    if (nzchar(mess <- naprint(x$na.action)))
        cat("  (", mess, ")\n", sep = "")
    cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)))
    cat("\n\nType of estimator:", x$type, get_type_description(x$type))
    cat("\n", "Number of Fisher Scoring iterations: ", x$iter, "\n", sep = "")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            } else {
                correl <- format(round(correl, 2L), nsmall = 2L,
                                 digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    invisible(x)
}

cg_steihaug_subproblem_with_H <- function(grad, H, Delta, 
                                          tol = 0.1, maxiter = 50) {
    p <- length(grad)
    z <- numeric(p)
    r <- -grad
    d <- r
    grad_norm <- sqrt(sum(grad^2))
    
    for (j in seq_len(maxiter)) {
        Bd <- drop(H %*% d)  # Use H directly
        dBd <- sum(d * Bd)
        
        if (dBd <= 0) {
            tau <- find_boundary_step(z, d, Delta)
            return(list(p = z + tau * d, on_boundary = TRUE, cg_iter = j))
        }
        
        r_norm_sq <- sum(r^2)
        alpha <- r_norm_sq / dBd
        z_new <- z + alpha * d
        
        if (sqrt(sum(z_new^2)) > Delta) {
            tau <- find_boundary_step(z, d, Delta)
            return(list(p = z + tau * d, on_boundary = TRUE, cg_iter = j))
        }
        
        z <- z_new
        r <- r - alpha * Bd
        
        if (sqrt(sum(r^2)) < tol * grad_norm) {
            return(list(p = z, on_boundary = FALSE, cg_iter = j))
        }
        
        r_norm_sq_new <- sum(r^2)
        beta <- r_norm_sq_new / r_norm_sq
        d <- r + beta * d
    }
    
    list(p = z, on_boundary = FALSE, cg_iter = maxiter)
}

compute_objective <- function(betas = betas, y = y, x = x, weights = weights,
                              offset = offset, family = family,
                              adjustment_function = adjustment_function,
                              fixed_totals = fixed_totals, row_totals = row_totals,
                              no_dispersion = no_dispersion, nobs = nobs, nvars = nvars,
                              keep = keep, need_qr = TRUE, need_hatvalues = TRUE) {

    # Define the objective function: f(beta) = ||adjusted_score(beta)||^2 / 2
    objective <- function(beta_vec) {
        theta <- c(beta_vec, 1)  # Fixed dispersion

        fit <- compute_fit(pars = theta, y = y, x = x, weights = weights,
                        offset = offset, family = family,
                        fixed_totals = fixed_totals, row_totals = row_totals,
                        no_dispersion = no_dispersion, nobs = nobs, nvars = nvars,
                        keep = keep, need_qr = TRUE, need_hatvalues = TRUE)
        adj <- adjustment_function(theta, fit = fit, level = 0,
                                    x = x, nobs = nobs, nvars = nvars,
                                    weights = weights)

        r <- fit$grad_beta + adj
        sum(r^2) / 2
    }

    residual <- function(beta_vec) {
        theta <- c(beta_vec, 1)  # Fixed dispersion

        fit <- compute_fit(pars = theta, y = y, x = x, weights = weights,
                    offset = offset, family = family,
                    fixed_totals = fixed_totals, row_totals = row_totals,
                    no_dispersion = no_dispersion, nobs = nobs, nvars = nvars,
                    keep = keep, need_qr = TRUE, need_hatvalues = TRUE)
        adj <- adjustment_function(theta, fit = fit, level = 0,
                                x = x, nobs = nobs, nvars = nvars,
                                weights = weights)

        r <- fit$grad_beta + adj
        r
    }

    theta <- c(betas, 1)
    fit <- compute_fit(pars = theta, y = y, x = x, weights = weights,
                        offset = offset, family = family,
                        fixed_totals = fixed_totals, row_totals = row_totals,
                        no_dispersion = no_dispersion, nobs = nobs, nvars = nvars,
                        keep = keep, need_qr = TRUE, need_hatvalues = TRUE)

    ## Computations using numDeriv
    F <- objective(betas)
    G <- numDeriv::grad(func = objective, x = betas)
    H_actual <- numDeriv::hessian(func = objective, x = betas)
    j <- numDeriv::jacobian(func = objective, x = betas)
    J <- numDeriv::jacobian(func = residual, x = betas)
    
    H_Gauss <- t(J) %*% J  # Gauss-Newton approximation to Hessian

    working_weights <- fit$working_weights

    info <- t(x) %*% (diag(working_weights) %*% x)
    H_approx <- info %*% info

    r <- residual(betas)

    #browser()

    list(value = F, gradient = G, jacobian = info, hessian = H_approx, residual = r)
} 