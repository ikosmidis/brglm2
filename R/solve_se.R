# Copyright (C) 2025- Ioannis Kosmidis, Federico Boiocchi, Philipp Sterzinger

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

#' Solve the MDYPL state evolution equations with or without
#' intercept, with signal strength or contaminated signal strength
#'
#' @param kappa asymptotic ratio of columns/rows of the design
#'     matrix. `kappa` should be in `(0, 1)`.
#' @param ss signal strength or corrupted signal strength, depending
#'     on whether `corrupted = TRUE` or not. See Details.
#' @param alpha the shrinkage parameter of the MDYPL
#'     estimator. `alpha` should be in `(0, 1)`.
#' @param intercept if `NULL` (default) then the MDYPL state evolution
#'     equations for the model with no intercept parameter are
#'     solved. If a real then the equations for the models with
#'     intercept parameter equal to `intercept` are solved. See
#'     Details.
#' @param start a vector with starting values for `mu`, `b`,`sigma`
#'     (and `iota` if `intercept` is numeric).
#' @param corrupted if `FALSE` (default) then `ss` is signal strength
#'     and `intercept`, if numeric, is the oracle intercept value. If
#'     `TRUE`, then `ss` is the corrupted signal strength, and
#'     `intercept`, if numeric, is the limit of the estimator computed
#'     by [mdyplFit()] with shrinkage parameter `alpha`. See Details.
#' @param gh A list with the Gauss-Hermite quadrature nodes and
#'     weights, as returned from `statmod::gauss.quad()` with `kind =
#'     "hermite"`. Default is `NULL`, in which case `gh` is set to
#'     `statmod::gauss.quad(200, kind = "hermite")`.
#' @param prox_tol tolerance for the computation of the proximal
#'     operator; default is `1e-10`.
#' @param transform if `TRUE` (default), the optimization is with
#'     respect to `log(mu)`, `log(b)`,`log(sigma)`, (and `iota` if
#'     `intercept` is numeric). If `FALSE`, then it is over `mu`, `b`,
#'     `sigma` (and `iota` if `intercept` is numeric). The solution is
#'     returned in terms of the latter set, regardless of how
#'     optimization took place.
#' @param init_iter how many iterations of [optim()] should we make to
#'     get starting values for [nleqslv::nleqslv()]? Default is `50`,
#'     but can also be `0` in which case `start` is directly passed to
#'     `nleqslv:nleqslv()`. `init_iter = "only"` results in only
#'     [optim()] being used. See Details.
#' @param init_method The method to be passed to [optim()]. Default is
#'     `"Nelder-Mead"`.
#' @param ... further arguments to be passed to [nleqslv::nleqslv()],
#'     unless `init_iter = "only"`, in which case `...` is further
#'     arguments to be passed to [optim()].
#'
#' @details
#'
#' `init_iter` iterations of [optim()] with `method = init_method` are
#' used towards minimizing `sum(se)^2`, where `se` is a vector of the
#' state evolution functions. The solution is then passed to
#' `nleqslv::nleqslv()` for a more aggressive iteration. The state
#' evolution equations are given in expressions (8) (model without
#' intercept) and expression (15) (model with intercept) in Sterzinger
#' & Kosmidis (2024).
#'
#' If `corrupted = FALSE` (default), then `ss` is the square root of
#' the signal strength, which is the limit \eqn{\gamma^2} of
#' \eqn{var(X \beta)}. If `corrupted = TRUE`, then `ss` is the square
#' root of the corrupted signal strength which is the limit
#' \eqn{\nu^2} of \eqn{var(X hat(beta)(\alpha))}, where
#' \eqn{hat(\beta)(\alpha)} is the maximum Diaconis-Ylvisaker prior
#' penalized likelihood (MDYPL) estimator as computed by [mdyplFit()]
#' with shrinkage parameter \eqn{alpha}.
#'
#' If `intercept = NULL`, then the state evolution equations are
#' solved for the model without intercept. If `intercept` is a real
#' number, then the state evolution equations for the model with
#' intercept are solved (i.e. with predictor \eqn{\eta_i = \theta +
#' x_i^T \beta}). In that case, what `intercept` represents depends on
#' the value of `corrupted`. If `corrupted = FALSE`, `intercept`
#' represents the oracle value of \eqn{\theta}, otherwise it represents
#' the limit `iota` of the MDYPL estimator of \eqn{\theta} as computed by
#' [mdyplFit()] with shrinkage parameter `alpha`.
#'
#' Note that `start` is always for `mu`, `b`,`sigma`, as is the
#' result, regardless whether `transform = TRUE` or
#' not. Transformations during optimization are done internally.
#'
#' @return
#'
#' If `intercept = NULL`, a vector with the values of `mu`,
#' `b`,`sigma`. Otherwise, a vector with the values of `mu`,
#' `b`,`sigma`, and `iota`, if `corrupted = FALSE`, or the value of
#' the intercept otherwise. The vector has attributes the state
#' evolution functions at the solution (`"funcs"`), the number of
#' iterations used by the last optimization method (`"iter"`), any
#' messages from the last optimization method (`"message"`), and
#' information on the optimization methods used
#' (`"optimization-chain"`).
#'
#' @author Ioannis Kosmidis `[aut, cre]` \email{ioannis.kosmidis@warwick.ac.uk}, Federico Boiocchi `[ctb]` \email{federico.boiocchi@gmail.com}, Philipp Sterzinger `[ctb, earlier Julia code by]` \email{P.Sterzinger@lse.ac.uk}
#'
#'
#' @references
#'
#' Zhao Q, Sur P, Cand\`es E J (2022). The asymptotic distribution of
#' the MLE in high-dimensional logistic models: Arbitrary
#' covariance. *Bernoulli*, **28**, 1835–1861. \doi{10.3150/21-BEJ1401}.
#'
#' Sterzinger P, Kosmidis I (2024). Diaconis-Ylvisaker prior
#' penalized likelihood for \eqn{p/n \to \kappa \in (0,1)} logistic
#' regression. *arXiv*:2311.07419v2, \url{https://arxiv.org/abs/2311.07419}.
#'
#'
#' @examples
#'
#' ## Reproducing Table 13 of Zhao et al. (2022, DOI: 10.3150/21-BEJ1401)
#' \dontrun{
#'
#' thetas <- c(0, 0.5, 1, 2, 2.5)
#' gamma0 <- 5
#' pars3 <- matrix(NA, length(thetas), 3)
#' pars4 <- matrix(NA, length(thetas), 4)
#' colnames(pars4) <- c("I_mu", "I_b", "I_sigma", "I_iota")
#' colnames(pars3) <- c("II_mu", "II_b", "II_sigma")
#' for (i in seq_along(thetas)) {
#'     start3 <- c(0.5, 1, 1)
#'     pars3[i, ] <- solve_se(kappa = 0.2, ss = sqrt(5 + thetas[i]^2),
#'                            alpha = 1, start = start3, init_iter = 0)
#'     start4 <- c(pars3[i, ], thetas[i])
#'     pars4[i, ] <- solve_se(kappa = 0.2, ss = sqrt(5), intercept = thetas[i],
#'                            alpha = 1, start = start4, init_iter = 0)
#' }
#'
#' cbind(pars3, pars4)
#'
#' }
#' @export
solve_se <- function(kappa, ss, alpha, intercept = NULL, start, corrupted = FALSE, gh = NULL, prox_tol = 1e-10, transform = TRUE, init_method = "Nelder-Mead", init_iter = 50, ...) {
    is_corrupted <- isTRUE(corrupted)
    init_solver <- if (is_corrupted) optim_se_corrupted else optim_se
    main_solver <- if (is_corrupted) nleqslv_se_corrupted else nleqslv_se
    if (is.null(gh))
        gh <- gauss.quad(200, kind = "hermite")
    if (isTRUE(init_iter == "only")) {
        start <- init_solver(kappa, ss, alpha, intercept, start, gh, prox_tol, method = init_method, ...)
        opt_chain <- paste0("optim(method = ", init_method, ")")
    } else {
        if (init_iter > 0) {
            start <- init_solver(kappa, ss, alpha, intercept, start, gh, prox_tol, method = init_method, control = list(maxit = init_iter));
            opt_chain <- paste0("optim(method = ", init_method, ", maxit = ", init_iter, ") -> ")
        } else {
            opt_chain <- ""
        }
        start <- main_solver(kappa, ss, alpha, intercept, start, gh, prox_tol, transform, ...)
        opt_chain <- paste0(opt_chain, "nleqslv()")
    }
    attr(start, "optimization-chain") <- opt_chain
    start
}

se_funcs <- function(kappa, ss, alpha, intercept = NULL, iota = NULL,
                     gh, prox_tol, corrupted = FALSE, transform = TRUE) {
    ## corrupted signal strength
    if (isTRUE(corrupted)) {
        if (is.null(iota)) {
            if (transform) {
                g <- function(pars) {
                    pars <- exp(pars)
                    mu <- pars[1]
                    b <- pars[2]
                    sigma <- pars[3]
                    gamma <- sqrt(ss^2 - kappa * sigma^2) / mu
                    if (is.na(gamma)) return(rep(NA, 3))
                    se0(mu = mu, b = b, sigma = sigma, kappa = kappa, gamma = gamma, alpha = alpha, gh = gh, prox_tol = prox_tol)
                }
            } else {
                g <- function(pars) {
                    mu <- pars[1]
                    b <- pars[2]
                    sigma <- pars[3]
                    gamma <- sqrt(ss^2 - kappa * sigma^2) / mu
                    if (is.na(gamma)) return(rep(NA, 3))
                    se0(mu = mu, b = b, sigma = sigma, kappa = kappa, gamma = gamma, alpha = alpha, gh = gh, prox_tol = prox_tol)
                }
            }
        } else {
            no_int <- 1:3
            if (transform) {
                g <- function(pars) {
                    pars[no_int] <- exp(pars[no_int])
                    mu <- pars[1]
                    b <- pars[2]
                    sigma <- pars[3]
                    gamma <- sqrt(ss^2 - kappa * sigma^2) / mu
                    if (is.na(gamma)) return(rep(NA, 4))
                    se1(mu = mu, b = b, sigma = sigma, iota = iota, kappa = kappa, gamma = gamma, alpha = alpha, intercept = pars[4], gh = gh, prox_tol = prox_tol)
                }
            } else {
                g <- function(pars) {
                    mu <- pars[1]
                    b <- pars[2]
                    sigma <- pars[3]
                    gamma <- sqrt(ss^2 - kappa * sigma^2) / mu
                    if (is.na(gamma)) return(rep(NA, 4))
                    se1(mu = mu, b = b, sigma = sigma, iota = iota, kappa = kappa, gamma = gamma, alpha = alpha, intercept = pars[4], gh = gh, prox_tol = prox_tol)
                }
            }
        }
    } else {
        ## signal strength
        if (is.null(intercept)) {
            if (transform) {
                g <- function(pars) {
                    pars <- exp(pars)
                    se0(mu = pars[1], b = pars[2], sigma = pars[3], kappa = kappa, gamma = ss, alpha = alpha, gh = gh, prox_tol = prox_tol)
                }
            } else {
                g <- function(pars) {
                    se0(mu = pars[1], b = pars[2], sigma = pars[3], kappa = kappa, gamma = ss, alpha = alpha, gh = gh, prox_tol = prox_tol)
                }
            }
        } else {
            no_int <- 1:3
            if (transform) {
                g <- function(pars) {
                    pars[no_int] <- exp(pars[no_int])
                    se1(mu = pars[1], b = pars[2], sigma = pars[3], iota = pars[4], kappa = kappa, gamma = ss, alpha = alpha, intercept = intercept, gh = gh, prox_tol = prox_tol)
                }
            } else {
                g <- function(pars) {
                    se1(mu = pars[1], b = pars[2], sigma = pars[3], iota = pars[4], kappa = kappa, gamma = ss, alpha = alpha, intercept = intercept, gh = gh, prox_tol = prox_tol)
                }
            }
        }
    }
    g
}

nleqslv_se <- function(kappa, gamma, alpha, intercept = NULL, start, gh = NULL, prox_tol = 1e-10, transform = TRUE, ...) {
    no_intercept <- is.null(intercept)
    g <- se_funcs(kappa, gamma, alpha, intercept, iota = NULL, gh, prox_tol, corrupted = FALSE, transform)
    npar <- 3 + !no_intercept
    stopifnot(length(start) == npar)
    start <- c(if (transform) log(start[1:3]) else start[1:3],
               if (no_intercept) NULL else start[4])
    res <- nleqslv(start, g, ...)
    if (transform) {
        if (no_intercept) {
            soln <- exp(res$x)
        } else {
            soln <- c(exp(res$x[1:3]), res$x[4])
        }
    } else {
        soln <- res$x
    }
    attr(soln, "funcs") <- res$fvec
    attr(soln, "iter") <- res$iter
    attr(soln, "message") <- res$message
    attr(soln, "nleqslv_termination_code") <- res$termcd
    soln
}

nleqslv_se_corrupted <- function(kappa, nu, alpha, iota = NULL, start, gh = NULL, prox_tol = 1e-10, transform = TRUE, ...) {
    no_intercept <- is.null(iota)
    g <- se_funcs(kappa, nu, alpha, intercept = NULL, iota, gh, prox_tol, corrupted = TRUE, transform)
    npar <- 3 + !no_intercept
    stopifnot(length(start) == npar)
    start <- c(if (transform) log(start[1:3]) else start[1:3],
               if (no_intercept) NULL else start[4])
    suppressWarnings(res <- nleqslv(start, g, ...))
    if (transform) {
        if (no_intercept) {
            soln <- exp(res$x)
        } else {
            soln <- c(exp(res$x[1:3]), res$x[4])
        }
    } else {
        soln <- res$x
    }
    attr(soln, "funcs") <- res$fvec
    attr(soln, "iter") <- res$iter
    attr(soln, "message") <- res$message
    attr(soln, "nleqslv_termination_code") <- res$termcd
    soln
}


optim_se <- function(kappa, gamma, alpha, intercept = NULL, start, gh = NULL, prox_tol = 1e-10, transform = TRUE, ...) {
    no_intercept <- is.null(intercept)
    g <- se_funcs(kappa, gamma, alpha, intercept, iota = NULL, gh, prox_tol, corrupted = FALSE, transform = TRUE)
    npar <- 3 + !no_intercept
    stopifnot(length(start) == npar)
    start <- c(log(start[1:3]), if (no_intercept) NULL else start[4])
    obj <- function(pars) {
        sum(g(pars)^2)
    }
    res <- optim(start, obj, ...)
    if (no_intercept) {
        soln <- exp(res$par)
    } else {
        soln <- c(exp(res$par[1:3]), res$par[4])
    }
    attr(soln, "funcs") <- g(soln)
    attr(soln, "iter") <- res$counts
    attr(soln, "convergence") <- res$convergence
    attr(soln, "message") <- res$message
    soln
}

optim_se_corrupted <- function(kappa, nu, alpha, iota = NULL, start, gh = NULL, prox_tol = 1e-10, transform = TRUE, ...) {
    no_intercept <- is.null(iota)
    g <- se_funcs(kappa, nu, alpha, intercept = NULL, iota, gh, prox_tol, corrupted = TRUE, transform = TRUE)
    npar <- 3 + !no_intercept
    stopifnot(length(start) == npar)
    start <- c(log(start[1:3]), if (no_intercept) NULL else start[4])
    obj <- function(pars) {
        sum(g(pars)^2)
    }
    suppressWarnings(res <- optim(start, obj, ...))
    if (no_intercept) {
        soln <- exp(res$par)
    } else {
        soln <- c(exp(res$par[1:3]), res$par[4])
    }
    suppressWarnings(attr(soln, "funcs") <- g(soln))
    attr(soln, "iter") <- res$counts
    attr(soln, "convergence") <- res$convergence
    attr(soln, "message") <- res$message
    soln
}

## solvers based on https://cran.r-project.org/package=trust
## trust_se <- function(kappa, gamma, alpha, intercept = NULL, start, gh = NULL, prox_tol = 1e-10, transform = FALSE, ...) {
##     ssq <- function(x) sum(x^2)
##     no_intercept <- is.null(intercept)
##     if (no_intercept) {
##         npars <- 3
##         g <- function(pars) {
##             pars <- exp(pars)
##             se0(mu = pars[1], b = pars[2], sigma = pars[3], kappa = kappa, gamma = gamma, alpha = alpha, gh = gh, prox_tol = prox_tol) |> ssq()
##         }
##         start <- log(start)
##     } else {
##         npars <- 4
##         stopifnot(length(start) == 4)
##         no_int <- 1:3
##         g <- function(pars) {
##             pars[no_int] <- exp(pars[no_int])
##             se1(mu = pars[1], b = pars[2], sigma = pars[3], iota = pars[4], kappa = kappa, gamma = gamma, alpha = alpha, intercept = intercept, gh = gh, prox_tol = prox_tol) |> ssq()
##         }
##         start[no_int] <- log(start[no_int])
##     }
##     h <- matrix(0, npars, npars)
##     upp_inds <- upper.tri(h, diag = TRUE)
##     low_inds <- lower.tri(h, diag = TRUE)
##     vec2mat <- function(vec, d) {
##         h <- matrix(NA, npars, npars)
##         h[upp_inds] <- vec
##         h[low_inds] <- t(h)[low_inds]
##         h
##     }
##     obj <- function(pars) {
##         v <- numDeriv::genD(g, pars)
##         list(value = v$f0,
##              gradient = v$D[1:npars],
##              hessian = vec2mat(v$D[-c(1:npars)]))
##     }
##     res <- trust(obj, start, rinit = 1, rmax = 5, ...)
##     if (no_intercept) {
##         soln <- exp(res$argument)
##     } else {
##         soln <- c(exp(res$argument[no_int]), res$argument[4])
##     }
##     attr(soln, "objective") <- res$value
##     attr(soln, "iter") <- res$iterations
##     soln
## }

## trust_se_est <- function(kappa, nu, alpha, iota = NULL, start, gh = NULL, prox_tol = 1e-10, transform = FALSE, ...) {
##     ssq <- function(x) sum(x^2)
##     no_intercept <- is.null(iota)
##     if (no_intercept) {
##         npars <- 3
##         g <- function(pars) {
##             pars <- exp(pars)
##             mu <- pars[1]
##             b <- pars[2]
##             sigma <- pars[3]
##             gamma <- sqrt(nu^2 - kappa * sigma^2) / mu
##             if (is.na(gamma)) return(NA)
##             se0(mu = mu, b = b, sigma = sigma, kappa = kappa, gamma = gamma, alpha = alpha, gh = gh, prox_tol = prox_tol) |> ssq()
##         }
##         start <- log(start)
##     } else {
##         npars <- 4
##         stopifnot(length(start) == 4)
##         no_int <- 1:3
##         g <- function(pars) {
##             pars[no_int] <- exp(pars[no_int])
##             mu <- pars[1]
##             b <- pars[2]
##             sigma <- pars[3]
##             gamma <- sqrt(nu^2 - kappa * sigma^2) / mu
##             if (is.na(gamma)) return(NA)
##             se1(mu = mu, b = b, sigma = sigma, iota = iota, kappa = kappa, gamma = gamma, alpha = alpha, intercept = pars[4], gh = gh, prox_tol = prox_tol) |> ssq()
##         }
##         start[no_int] <- log(start[no_int])
##     }
##     h <- matrix(0, npars, npars)
##     upp_inds <- upper.tri(h, diag = TRUE)
##     low_inds <- lower.tri(h, diag = TRUE)
##     vec2mat <- function(vec, d) {
##         h <- matrix(NA, npars, npars)
##         h[upp_inds] <- vec
##         h[low_inds] <- t(h)[low_inds]
##         h
##     }
##     obj <- function(pars) {
##         v <- numDeriv::genD(g, pars)
##         list(value = v$f0,
##              gradient = v$D[1:npars],
##              hessian = vec2mat(v$D[-c(1:npars)]))
##     }
##     res <- trust(obj, start, rinit = 1, rmax = 5, ...)
##     if (no_intercept) {
##         soln <- exp(res$argument)
##     } else {
##         soln <- c(exp(res$argument[no_int]), res$argument[4])
##     }
##     attr(soln, "objective") <- res$value
##     attr(soln, "iter") <- res$iterations
##     soln
## }
