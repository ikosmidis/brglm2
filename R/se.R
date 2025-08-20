#' MDYPL state evolution functions with no intercept
#'
#' @param mu aggregate bias parameter.
#' @param b parameter `b` in the state evolution functions.
#' @param sigma square root of the aggregate variance of the MDYPL
#'     estimator.
#' @param kappa asymptotic ratio of columns/rows of the design
#'     matrix. `kappa` should be in `(0, 1)`.
#' @param gamma the square root of the limit of the variance of the
#'     linear predictor.
#' @param alpha the shrinkage parameter of the MDYPL
#'     estimator. `alpha` should be in `(0, 1)`.
#' @param gh A list with the Gauss-Hermite quadrature nodes and
#'     weights, as returned from `statmod::gauss.quad()` with `kind =
#'     "hermite"`. Default is `NULL`, in which case `gh` is set to
#'     `statmod::gauss.quad(200, kind = "hermite")` is used.
#' @param prox_tol tolerance for the computation of the proximal
#'     operator; default is `1e-10`.
#'
#' @export
se0 <- function(mu, b, sigma, kappa, gamma, alpha, gh = NULL, prox_tol = 1e-10) {
    if (is.null(gh))
        gh <- gauss.quad(200, kind = "hermite")
    a_frac <- 0.5 * (1 + alpha)
    xi <- gh$nodes
    wi <- gh$weights
    n_quad <- length(xi)
    q1 <- rep(sqrt(2) * gamma * xi, times = n_quad)
    q2 <- mu * q1 + rep(sqrt(2 * kappa) * sigma * xi, each = n_quad)
    w2 <- rep(wi, times = n_quad) * rep(wi, each = n_quad)

    w2p <- 2 * w2 * plogis2(q1)  / pi
    p_prox <- plogis2(prox(q2 + a_frac * b, b, prox_tol))
    prox_resid <- a_frac - p_prox
    c(sum(w2p * q1 * prox_resid),
      1 - kappa - sum(w2p / (1 + b * p_prox * (1 - p_prox))),
      kappa^2 * sigma^2 - b^2 * sum(w2p * prox_resid^2))
}

#' MDYPL state evolution functions with intercept
#'
#' @param mu aggregate bias parameter.
#' @param b parameter `b` in the state evolution functions.
#' @param sigma square root of the aggregate variance of the MDYPL
#'     estimator.
#' @param iota limits of the MDYPL estimate for the intercept as the sample size goes to +Inf
#' @param kappa asymptotic ratio of columns/rows of the design
#'     matrix. `kappa` should be in `(0, 1)`.
#' @param gamma the square root of the limit of the variance of the
#'     linear predictor.
#' @param alpha the shrinkage parameter of the MDYPL
#'     estimator. `alpha` should be in `(0, 1)`.
#' @param intercept intercept of the logistic regression model
#' @param gh A list with the Gauss-Hermite quadrature nodes and
#'     weights, as returned from `statmod::gauss.quad()` with `kind =
#'     "hermite"`. Default is `NULL`, in which case `gh` is set to
#'     `statmod::gauss.quad(200, kind = "hermite")` is used.
#' @param prox_tol tolerance for the computation of the proximal
#'     operator; default is `1e-10`. fixed point problem solved via Newton-Raphson
#' @export
se1 <- function(mu, b, sigma, iota, kappa, gamma, alpha, intercept, gh = NULL, prox_tol = 1e-10) {
    if (is.null(gh))
        gh <- gauss.quad(200, kind = "hermite")
    a_frac <- 0.5 * (1 + alpha)
    xi <- gh$nodes
    wi <- gh$weights
    n_quad <- length(xi)
    q1n <- sqrt(2) * gamma * xi
    q1 <- rep(q1n + intercept, times = n_quad)
    q2 <- mu * q1n + rep(sqrt(2 * kappa) * sigma * xi + iota, each = n_quad)
    w2 <- rep(wi, times = n_quad) * rep(wi, each = n_quad)

    p_q1_pos <- plogis2(q1)
    w2pi <- w2 / pi
    w2p_pos <- w2pi * p_q1_pos
    w2p_neg <- w2pi - w2p_pos

    p_prox_pos <- plogis2(prox(a_frac * b + q2, b, prox_tol))
    p_prox_neg <- plogis2(prox(a_frac * b - q2, b, prox_tol))

    prox_pos_resid <- a_frac - p_prox_pos
    prox_neg_resid <- a_frac - p_prox_neg

    c(sum(q1 * (w2p_pos * prox_pos_resid - w2p_neg * prox_neg_resid)),
      1 - kappa - sum(w2p_pos / (1 + b * p_prox_pos * (1 - p_prox_pos)) + w2p_neg / (1 + b * p_prox_neg * (1 - p_prox_neg))),
      kappa^2 * sigma^2 - b^2 * sum(w2p_pos * prox_pos_resid^2 + w2p_neg * prox_neg_resid^2),
      sum(w2p_pos * prox_pos_resid - w2p_neg * prox_neg_resid))
}


