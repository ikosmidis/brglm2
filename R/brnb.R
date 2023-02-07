# Copyright (C) 2020- Euloge Clovis Kenne Pagui, Ioannis Kosmidis

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

#' Bias reduction for negative binomial regression models
#'
#' [brnb()] is a function that fits negative binomial regression
#' models using implicit and explicit bias reduction methods.
#'
#' @inheritParams stats::glm
#' @param link The link function. Currently must be one of `"log"`,
#'     `"sqrt"` or `"identity"`.
#' @param control a list of parameters for controlling the fitting
#'     process. See [brglmControl()] for details.
#' @return A fitted model object of class [`"brnb"`][brnb] inheriting
#'     from [`"negbin"`][negbin] and [`"brglmFit"`][brglmFit]. The
#'     object is similar to the output of [brglmFit()] but contains
#'     four additional components: `theta` for the maximum likelihood
#'     estimate of the dispersion parameter as in [MASS::glm.nb()],
#'     `vcov.mean` for the estimated variance-covariance matrix of the
#'     regression coefficients, `vcov.dispersion` for the estimated
#'     variance of the dispersion parameter in the chosen
#'     parameterization (using the expected information), and
#'     `twologlik` for twice the log-likelihood function.
#'
#' @details
#'
#' A detailed description of the fitting procedure is given in the
#' iteration vignette (see, `vignette("iteration", "brglm2")` and
#' Kosmidis et al, 2020). The number of iterations when estimating
#' parameters are controlled by the `maxit` argument of
#' [brglmControl()].
#'
#' The type of score adjustment to be used is specified through the
#' `type` argument (see [brglmControl()] for details).
#'
#' The available options are:
#'
#' * `type = "AS_mixed"`: the mixed bias-reducing score
#' adjustments in Kosmidis et al (2020) that result in mean bias
#' reduction for the regression parameters and median bias reduction
#' for the dispersion parameter, if any; default.
#'
#' * `type = "AS_mean"`: the mean bias-reducing score
#' adjustments in Firth (1993) and Kosmidis & Firth (2009).
#'
#' * `type = "AS_median"`: the median bias-reducing score
#' adjustments in Kenne Pagui et al. (2017)
#'
#' * `type = "MPL_Jeffreys"`: maximum penalized likelihood
#' with powers of the Jeffreys prior as penalty.
#'
#' * `type = "ML"`: maximum likelihood.
#'
#' * `type = "correction"`: asymptotic bias correction, as in
#' Cordeiro & McCullagh (1991).
#'
#' The choice of the parameterization for the dispersion is controlled
#' by the `transformation` argument (see [brglmControl()] for
#' details).  The default is `"identity"`. Using `transformation =
#' "inverse"` uses the dispersion parameterization that
#' [MASS::glm.nb()] uses.
#'
#' @author Euloge Clovis Kenne Pagui `[aut]` \email{kenne@stat.unipd.it}, Ioannis Kosmidis `[aut, cre]` \email{ioannis.kosmidis@warwick.ac.uk}
#'
#' @references
#'
#' Cordeiro G M, McCullagh P (1991). Bias correction in generalized
#' linear models. *Journal of the Royal Statistical Society. Series B
#' (Methodological)*, **53**, 629-643. \doi{10.1111/j.2517-6161.1991.tb01852.x}.
#'
#' Firth D (1993). Bias reduction of maximum likelihood estimates.
#' *Biometrika*. **80**, 27-38. \doi{10.2307/2336755}.
#'
#' Kenne Pagui E C, Salvan A, Sartori N (2017). Median bias
#' reduction of maximum likelihood estimates. *Biometrika*, **104**,
#' 923–938. \doi{10.1093/biomet/asx046}.
#'
#' Kosmidis I, Kenne Pagui E C, Sartori N (2020). Mean and median bias
#' reduction in generalized linear models. *Statistics and Computing*,
#' **30**, 43-59. \doi{10.1007/s11222-019-09860-6}.
#'
#' Kosmidis I, Firth D (2009). Bias reduction in exponential family
#' nonlinear models. *Biometrika*, **96**, 793-804. \doi{10.1093/biomet/asp055}.
#'
#' @examples
#' ## Example in Saha, K., & Paul, S. (2005). Bias-corrected maximum
#' ## likelihood estimator of the negative binomial dispersion
#' ## parameter.  Biometrics, 61, 179--185.
#' #
#' # Number of revertant colonies of salmonella data
#' salmonella <- data.frame(freq = c(15, 16, 16, 27, 33, 20,
#'                                   21, 18, 26, 41, 38, 27,
#'                                   29, 21, 33, 60, 41, 42),
#'                          dose = rep(c(0, 10, 33, 100, 333, 1000), 3),
#'                          observation = rep(1:3, each = 6))
#'
#' # Maximum likelihood fit with glm.nb of MASS
#' salmonella_fm <- freq ~ dose + log(dose + 10)
#' fitML_glmnb <- MASS::glm.nb(salmonella_fm, data = salmonella)
#'
#' # Maximum likelihood fit with brnb
#' fitML <- brnb(salmonella_fm, data = salmonella,
#'               link = "log", transformation = "inverse", type = "ML")
#'
#' # Mean bias-reduced fit
#' fitBR_mean <- update(fitML, type = "AS_mean")
#'
#' # Median bias-reduced fit
#' fitBR_median <- update(fitML, type = "AS_median")
#'
#' # Mixed bias-reduced fit
#' fitBR_mixed <- update(fitML, type = "AS_mixed")
#'
#' # Mean bias-corrected fit
#' fitBC_mean <- update(fitML, type = "correction")
#'
#' # Penalized likelihood with Jeffreys-prior penalty
#' fit_Jeffreys <- update(fitML, type = "MPL_Jeffreys")
#'
#' # The parameter estimates from glm.nb and brnb with type = "ML" are
#' # numerically the same
#' all.equal(c(coef(fitML_glmnb), fitML_glmnb$theta),
#'             coef(fitML, model = "full"), check.attributes = FALSE)
#'
#' # Because of the invariance properties of the maximum likelihood,
#' # median reduced-bias, and mixed reduced-bias estimators the
#' # estimate of a monotone function of the dispersion should be
#' # (numerically) the same as the function of the estimate of the
#' # dispersion:
#'
#' # ML
#' coef(fitML, model = "dispersion")
#' 1 / coef(update(fitML, transformation = "identity"), model = "dispersion")
#'
#' # Median BR
#' coef(fitBR_median, model = "dispersion")
#' 1 / coef(update(fitBR_median, transformation = "identity"), model = "dispersion")
#'
#' # Mixed BR
#' coef(fitBR_mixed, model = "dispersion")
#' 1 / coef(update(fitBR_mixed, transformation = "identity"), model = "dispersion")
#'
#' ## The same is not true for mean BR
#' coef(fitBR_mean, model = "dispersion")
#' 1 / coef(update(fitBR_mean, transformation = "identity"), model = "dispersion")
#'
#' \donttest{
#' ## An example  from Venables & Ripley (2002, p.169).
#' data("quine", package = "MASS")
#' quineML <- brnb(Days ~ Sex/(Age + Eth*Lrn), link = "sqrt", transformation="inverse",
#'                 data = quine, type="ML")
#' quineBR_mean <- update(quineML, type = "AS_mean")
#' quineBR_median <- update(quineML, type = "AS_median")
#' quineBR_mixed <- update(quineML, type = "AS_mixed")
#' quine_Jeffreys <- update(quineML, type = "MPL_Jeffreys")
#'
#' fits <- list(ML = quineML,
#'              AS_mean = quineBR_mean,
#'              AS_median = quineBR_median,
#'              AS_mixed = quineBR_mixed,
#'              MPL_Jeffreys = quine_Jeffreys)
#' sapply(fits, coef, model = "full")
#' }
#'
#' @export
brnb <- function(formula, data, subset, weights = NULL, offset = NULL,
                 link = "log", start = NULL, etastart = NULL,
                 mustart = NULL,  control = list(...), na.action,
                 model = TRUE, x = FALSE, y = TRUE, contrasts = NULL,
                 intercept = TRUE, singular.ok = TRUE ,...) {

    loglik <- function(n, th, mu, y, w) {
        sum(w * (lgamma(th + y) - lgamma(th) - lgamma(y + 1) +
                 th * log(th) + y * log(mu + (y == 0)) - (th + y) * log(th + mu)))
    }

    trace_iteration <- function() {
        if (control$trace) {
            st <- max(abs(step_beta), na.rm = TRUE)
            gr <- max(abs(adjusted_grad_beta), na.rm = TRUE)
            cat("Coefficients update:\t", sprintf("%03f", betas), "\n", sep = " ")
            cat("Outer/Inner iteration:\t", sprintf("%03d", iter),
                "/", sprintf("%03d", step_factor), "\n", sep = "")
            cat("max |step|:", format(round(st, 6), nsmall = 6,
                                      scientific = FALSE), "\t", "max |gradient|:",
                format(round(gr, 6), nsmall = 6, scientific = FALSE),
                "\n")
            st <- abs(step_dispersion)
            gr <- abs(adjusted_grad_dispersion)
            cat("Dispersion update:\t",sprintf("%03f", dispersion), "\n", sep = "")
            cat("Outer iteration:\t", sprintf("%03d", iter),
                "\n")
            cat("max |step|:", format(round(st, 6), nsmall = 6,
                                      scientific = FALSE), "\t", "max |gradient|:",
                format(round(gr, 6), nsmall = 6, scientific = FALSE),
                "\n")
        }
    }

    mf <- Call <- match.call()
    m <- match(c("formula", "data", "subset", "weights", "na.action",
                 "etastart", "mustart", "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval.parent(mf)
    Terms <- attr(mf, "terms")
    return_x <- x
    return_y <- y

    y <- model.response(mf, "numeric")
    x <- if (!is.empty.model(Terms))
             model.matrix(Terms, mf, contrasts)
         else matrix(, NROW(y), 0)
    weights <- model.weights(mf)
    offset <- model.offset(mf)
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    x <- as.matrix(x)
    nobs <- NROW(x)
    nvars <- NCOL(x)
    if (missing_offset <- is.null(offset)) {
        offset <- rep.int(0, nobs)
    }

    if (is.null(weights)) weights <- rep(1,nobs)
    if (any(weights < 0))
        stop("negative weights not allowed")

    control <- do.call("brglmControl", control)
    ok_links <- c("log", "sqrt", "identity")
    if (!isTRUE(link %in% ok_links))
        stop(link, " is not one of the appropriate link function")
    ok_transformations <- c("log", "sqrt", "identity", "inverse")
    if (!isTRUE(control$transformation %in% ok_transformations))
        stop(control$transformation, " is not one of the implemented dispersion transformations")

    linkobj <- enrichwith::enrich(make.link(link),with="d2mu.deta")
    linkfun <- linkobj$linkfun
    linkinv <- linkobj$linkinv
    dmu.deta <- linkobj$mu.eta
    d2mu.deta <- linkobj$d2mu.deta

    linkobj_disp <- enrichwith::enrich(make.link(control$transformation), with = "d2mu.deta")
    linkfun_disp <- linkobj_disp$linkfun
    linkinv_disp <- linkobj_disp$linkinv
    dkappa.dphi <- linkobj_disp$mu.eta
    d2kappa.dphi <- linkobj_disp$d2mu.deta

    fam0 <- do.call("poisson", list(link = link))

    is_ML <- control$type == "ML"
    is_AS_median <- control$type == "AS_median"
    is_AS_mixed <- control$type == "AS_mixed"
    is_correction <- control$type == "correction"

    ## computation of s1
    s1fun <- function(yy, k) {
        res <- 0
        if (yy > 0) {
            j <- 0:(yy - 1)
            res <- sum(j / (k * j + 1))
        }
        res
    }
    ## Computation of needed expectations almost exactly
    exp_quant <- function(mu, k) {
        n <- length(mu)
        E_s2 <-  E_s2y <- E_s1s2 <- E_s3 <- numeric(n)
        if (k < 0) stop("negative value of k")
        ymax <- max(qnbinom(1 - 10 * .Machine$double.eps, mu = mu, size = 1 / k))
        ## if (ymax > 30000) stop("ymax too much large")  ## need control on largest value?
        if (ymax > 1) {
            out <- .C('expectedValues',
                      as.double(mu),
                      as.double(k),
                      as.integer(ymax),
                      as.integer(n),
                      E_s2 = as.double(E_s2),
                      E_s2y = as.double(E_s2y),
                      E_s1s2 = as.double(E_s1s2),
                      E_s3 = as.double(E_s3))
            E_s2 <- out$E_s2
            E_s2y <- out$E_s2y
            E_s1s2 <- out$E_s1s2
            E_s3 <- out$E_s3
        }
        list(E_s2 = E_s2, E_s2y = E_s2y, E_s1s2 = E_s1s2, E_s3 = E_s3)
    }

    ## ## ## IK, 20 Apr 2021: Optimized R version of exp_quant
    ## exp_quant <- function(mu, k) {
    ##     if (k < 0) {
    ##         stop("encountered negative value of k")
    ##     }
    ##     E_s2 <-  E_s2y <- E_s1s2 <- E_s3 <- numeric(length(mu))
    ##     ymax <- max(qnbinom(1 - 10 * .Machine$double.eps, mu = mu, size = 1 / k), na.rm = TRUE)
    ##     if (ymax > 1) {
    ##         yval <- 0:ymax
    ##         pmat <- sapply(mu, function(x) dnbinom(yval, mu = x, size = 1 / k))
    ##         j <- c(0, 0:(ymax - 1))
    ##         fra <- j / (k * j + 1)
    ##         s1 <- cumsum(fra)
    ##         s2 <- cumsum(fra^2)
    ##         s3 <- cumsum(fra^3)
    ##         temp <- pmat * s2
    ##         E_s2 <- .colSums(temp, ymax + 1, nobs, TRUE)
    ##         E_s2y <- .colSums(temp * yval, ymax + 1, nobs, TRUE)
    ##         E_s1s2 <- .colSums(temp * s1, ymax + 1, nobs, TRUE)
    ##         E_s3 <- .colSums(pmat * s3, ymax + 1, nobs, TRUE)
    ##     }
    ##     list(E_s2 = E_s2, E_s2y = E_s2y, E_s1s2 = E_s1s2, E_s3 = E_s3)
    ## }

    ## required quantities
    key_quantities <- function(pars) {
        betas <- pars[1:nvars]
        phi <- pars[nvars + 1]
        etas <- drop(x %*% betas + offset)
        mus <- linkinv(etas)
        d1 <- dmu.deta(etas)
        d2 <- d2mu.deta(etas)
        k <- linkinv_disp(phi)
        d1phi <- dkappa.dphi(phi)
        d2phi <- d2kappa.dphi(phi)
        varmus <- k * mus^2 + mus  ## variance
        d1varmus <- 2 * k * mus + 1  ## derivative of variance
        working_weights <- weights * d1^2 / varmus
        wx <- sqrt(working_weights) * x
        qr_decomposition <- qr(wx)
        s1 <- sapply(y, function(t) s1fun(t, k))
        expectations <- exp_quant(mus, k)
        out <- c(list(betas = betas,
                      k = k,
                      etas = etas,
                      phi = phi,
                      mus = mus,
                      d1 = d1,
                      d2 = d2,
                      d1phi = d1phi,
                      d2phi = d2phi,
                      working_weights = working_weights,
                      varmus = varmus,
                      d1varmus = d1varmus,
                      s1 = s1,
                      qr_decomposition = qr_decomposition),expectations)
        out
    }

    ## hat values
    hat_values <- function(pars, fit = NULL) {
        if (is.null(fit)) {
            fit <- key_quantities(pars)
        }
        with(fit, {
            Qmat <- qr.Q(qr_decomposition)
            .rowSums(Qmat * Qmat, nobs, nvars, TRUE)
        })
    }

    ## score function
    score <- function(pars, level = 0, fit = NULL) {
        if (is.null(fit)) {
            fit <- key_quantities(pars)
        }
        with(fit, {
            if (level == 0) {
                return(.colSums(x * working_weights * (y - mus) / d1, nobs, nvars, TRUE))
            }
            if (level == 1) {
                return(sum((s1 - mus * y / (k * mus + 1) + ((k * mus + 1) * log(k * mus + 1) - k * mus) / (k^3 * mus + k^2)) * d1phi * weights, na.rm = TRUE))
            }
        })
    }

    ## information matrix
    information <- function(pars, level = 0, inverse = FALSE, fit = NULL) {
        if (is.null(fit)) {
            fit <- key_quantities(pars)
        }
        with(fit, {
            if (level == 0) {
                R_matrix <- qr.R(qr_decomposition)
                if (inverse) {
                    return(chol2inv(R_matrix))
                } else {
                    return(crossprod(R_matrix))
                }
            }
            if (level == 1) {
                iphiphi <- sum( weights * (E_s2 + ((2 * k * mus + 2) * log((k * mus + 1)) - k^2 * mus^2 - 2 * k * mus) / (k^4 * mus + k^3)) * d1phi^2, na.rm = TRUE )
                if (inverse) {
                    return(1/iphiphi)
                } else {
                    return(iphiphi)
                }
            }
        })
    }

    ## Jeffreys adjustment
    AS_Jeffreys_adjustment <- function(pars, level = 0, fit = NULL) {
        if (is.null(fit)) {
            fit <-  key_quantities(pars)
        }
        with(fit, {
            iphiphi <- sum(weights * (E_s2 + ((2 * k * mus + 2) * log((k * mus + 1)) - k^2 * mus^2 - 2 * k * mus) / (k^4 * mus + k^3)) * d1phi^2, na.rm = TRUE )
            hatvalues <- hat_values(pars, fit = fit)
            if (level == 0) {
                ## Use only observations with keep = TRUE to ensure that no division with zero takes place
                return(control$a * ( .colSums(hatvalues * (2 * d2/d1 - d1varmus * d1 / varmus) * x, nobs, nvars, TRUE) +
                                     .colSums(d1phi^2 * x * working_weights * (E_s2y - mus * E_s2 - mus^3/(k * mus + 1)) / d1, nobs, nvars, TRUE) /  iphiphi))
            }
            if (level == 1) {
                return(control$a * (sum(-d1phi * hatvalues * mus^2 / varmus, na.rm = TRUE) +
                                    sum(weights * (-2 * E_s3  + (2 * k^3 * mus^3 + 9 * k^2 * mus^2 + 6 * k * mus - 6 * (k * mus + 1)^2 * log(k * mus + 1)) / (k^4 * (k * mus + 1)^2) + E_s1s2 -
                                                   mus * E_s2y / (k * mus + 1) -  (k * mus-(k * mus + 1) * log(k * mus + 1)) * E_s2 / (k^2 * (k * mus + 1))), na.rm = TRUE) * d1phi^3 +
                                    2 * sum(weights * (E_s2 + ((2 * k * mus + 2) * log((k * mus + 1)) - k^2 * mus^2 - 2 * k * mus)/(k^4 * mus + k^3)), na.rm = TRUE) * d1phi * d2phi / iphiphi))
            }
        })
    }

    ## mean adjustment
    AS_mean_adjustment <- function(pars, level = 0, fit = NULL) {
        if (is.null(fit)) {
            fit <- key_quantities(pars)
        }
        with(fit, {
            hatvalues <- hat_values(pars, fit = fit)
            if (level == 0) {
                return(.colSums(0.5 * hatvalues * d2 / d1 * x, nobs, nvars, TRUE))
            }
            if (level == 1) {
                iphiphi <- sum(weights * (E_s2 + ((2 * k * mus + 2) * log((k * mus + 1)) - k^2 * mus^2 - 2 * k * mus) / (k^4 * mus + k^3)) * d1phi^2, na.rm = TRUE)
                pqphiphi<- sum(weights * (-2 * E_s3 + (2 * k^3 * mus^3 + 9 * k^2 * mus^2 + 6 * k * mus - 6 * (k * mus + 1)^2 * log(k * mus + 1)) / (k^4 * (k * mus + 1)^2) +
                                          2 * E_s1s2 - 2 * mus * E_s2y / (k * mus + 1) - 2 * (k * mus - (k * mus + 1) * log(k * mus + 1)) * E_s2 / (k^2 * (k * mus + 1))), na.rm = TRUE) * d1phi^3 +
                           sum(weights * (E_s2 + ((2 * k * mus + 2) * log((k * mus + 1)) - k^2 * mus^2 - 2 * k * mus) / (k^4 * mus + k^3)) , na.rm = TRUE) * d1phi * d2phi
                return(0.5 * d1phi * sum(weights * hatvalues * d1^2 * mus^2 / (working_weights * varmus^2), na.rm = TRUE) + 0.5 * pqphiphi / iphiphi)
            }
        })
    }

    ## median adjustment
    AS_median_adjustment <- function(pars, level = 0, fit = NULL) {
        if (is.null(fit)) {
            fit <- key_quantities(pars)
        }
        with(fit, {
            hatvalues <- hat_values(pars, fit = fit)
            if (level == 0) {
                invInfo <- information(pars = pars, level = 0, inverse = TRUE, fit = fit)
                info <- information(pars = pars, level = 0, inverse = FALSE, fit = fit)
                F2_tilde <- numeric(nvars)
                for (j in seq.int(nvars)) {
                    invInfo_j <- invInfo[j,]
                    vcov_j <- tcrossprod(invInfo_j) / invInfo_j[j]
                    hats_j <- .rowSums((x %*% vcov_j) * x, nobs, nvars, TRUE) * working_weights
                    F2_tilde[j] <- invInfo_j %*% .colSums(x * (hats_j * (d1 * d1varmus / (6 * varmus) - 0.5 * d2 / d1)), nobs, nvars,  TRUE)
                }
                return( .colSums(0.5 * hatvalues * d2 / d1 * x, nobs, nvars, TRUE) + info %*% F2_tilde)
            }
            if (level == 1) {
                iphiphi <- sum(weights * (E_s2 + ((2 * k * mus + 2) * log((k * mus + 1)) - k^2 * mus^2 - 2 * k * mus) / (k^4 * mus + k^3)) * d1phi^2, na.rm = TRUE)
                pqphiphi <- sum(weights * (-2 * E_s3 + (2 * k^3 * mus^3 + 9 * k^2 * mus^2 + 6 * k * mus - 6 * (k * mus + 1)^2 * log(k * mus + 1)) / (k^4 * (k * mus + 1)^2) +
                                           2 * E_s1s2 - 2 * mus * E_s2y / (k * mus + 1) - 2 * (k * mus - (k * mus + 1) * log(k * mus + 1)) * E_s2 / (k^2 * (k * mus + 1))), na.rm = TRUE) * d1phi^3 +
                            sum(weights *(E_s2 + ((2 * k * mus + 2) * log((k * mus + 1)) - k^2 * mus^2 - 2 * k * mus) / (k^4 * mus + k^3)), na.rm = TRUE) * d1phi * d2phi
                pq2phiphi <- (sum(weights * (-2 * E_s3 / 3 + (2 * k^3 * mus^3 + 9 * k^2 * mus^2 + 6 * k * mus - 6 * (k * mus + 1)^2 * log(k * mus + 1)) / (3 * k^4 * (k * mus + 1)^2) +
                                             E_s1s2 / 2 - 0.5 * mus * E_s2y / (k * mus + 1) - 0.5 * (k * mus-(k * mus + 1) * log(k * mus + 1)) * E_s2 / (k^2 * (k * mus + 1))), na.rm = TRUE) * d1phi^3 +
                              0.5 * sum(weights * (E_s2 + ((2 * k * mus + 2) * log((k * mus + 1)) - k^2 * mus^2 - 2 * k * mus)/(k^4 * mus + k^3)), na.rm = TRUE) * d1phi * d2phi)
                return((0.5 * d1phi * sum(weights * hatvalues * d1^2 * mus^2 / (working_weights * varmus^2), na.rm = TRUE) + 0.5 * pqphiphi / iphiphi - pq2phiphi / iphiphi))
            }
        })
    }

    ## mixed adjustment
    AS_mixed_adjustment <- function(pars, level = 0, fit = NULL) {
        if (is.null(fit)) {
            fit <- key_quantities(pars)
        }
        with(fit, {
            hatvalues <- hat_values(pars, fit = fit)
            if (level == 0) {
                return(.colSums(0.5 * hatvalues * d2 / d1 * x, nobs, nvars, TRUE))
            }
            if (level == 1)
            {
                iphiphi <- sum(weights * (E_s2 + ((2 * k * mus + 2) * log((k * mus + 1)) - k^2 * mus^2 - 2 * k * mus) / (k^4 * mus + k^3)) * d1phi^2, na.rm = TRUE)
                pqphiphi <- (sum(weights *( -2 * E_s3 + (2 * k^3 * mus^3 + 9 * k^2 * mus^2 + 6 * k * mus - 6 * (k * mus + 1)^2 * log(k * mus + 1)) / (k^4 * (k * mus + 1)^2) +
                                            2 * E_s1s2 - 2 * mus * E_s2y / (k * mus + 1) - 2 * (k * mus - (k * mus + 1) * log(k * mus + 1)) * E_s2 / (k^2 * (k * mus + 1))), na.rm = TRUE) * d1phi^3 +
                             sum(weights * (E_s2 + ((2 * k * mus + 2) * log((k * mus + 1)) - k^2 * mus^2 - 2 * k * mus) / (k^4 * mus + k^3)), na.rm = TRUE) * d1phi * d2phi)
                pq2phiphi <- (sum(weights *(-2 * E_s3 / 3 + (2 * k^3 * mus^3 + 9 * k^2 * mus^2 + 6 * k * mus - 6 * (k * mus + 1)^2 * log(k * mus + 1)) / (3 * k^4 * (k * mus + 1)^2) +
                                            E_s1s2 / 2 - 0.5 * mus * E_s2y / (k * mus + 1) - 0.5 * (k * mus - (k * mus + 1) * log(k * mus + 1)) * E_s2 / (k^2 * (k * mus + 1))), na.rm = TRUE) * d1phi^3 +
                              0.5 * sum(weights * (E_s2 + ((2 * k * mus + 2) * log((k * mus + 1)) - k^2 * mus^2 - 2 * k * mus) / (k^4 * mus + k^3)), na.rm = TRUE) * d1phi * d2phi)
                return((0.5 * d1phi * sum(weights * hatvalues * d1^2 * mus^2 / (working_weights * varmus^2), na.rm = TRUE) + 0.5 * pqphiphi / iphiphi - pq2phiphi / iphiphi))
            }
        })
    }

    ## adjustment term function ##
    adjustment_function <- switch(control$type,
                                  AS_mean = AS_mean_adjustment,
                                  AS_median = AS_median_adjustment,
                                  AS_mixed = AS_mixed_adjustment,
                                  MPL_Jeffreys = AS_Jeffreys_adjustment,
                                  correction = function(pars, ...) 0,
                                  ML = function(pars, ...) 0)

    ## required components for computing the adjusted scores

    compute_step_components <- function(pars, level = 0, fit = NULL) {
        if (is.null(fit)) {
            fit <- key_quantities(pars)
        }
        if (level == 0) {
            grad <- score(pars, level = 0, fit = fit)
            inverse_info <- try(information(pars, level = 0, inverse = TRUE,fit = fit))
            failed_inversion <- inherits(inverse_info, "try-error")
            adjustment <- adjustment_function(pars, level = 0, fit = fit)
            failed_adjustment <- any(is.na(adjustment))
        }
        if (level == 1) {
            grad <- score(pars, level = 1, fit = fit)
            inverse_info <- try(information(pars, level = 1, inverse = TRUE,fit = fit))
            failed_inversion <- inherits(inverse_info, "try-error")
            adjustment <- adjustment_function(pars, level = 1, fit = fit)
            failed_adjustment <- any(is.na(adjustment))
        }
        out <- list(grad = grad,
                    inverse_info = inverse_info,
                    adjustment = adjustment,
                    failed_inversion = failed_inversion,
                    failed_adjustment = failed_adjustment)
        out
    }

    betas_names <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y))
                  rownames(y)
              else names(y)
    EMPTY <- nvars == 0

    valid_eta <- fam0$valideta
    valid_mu <- fam0$validmu
    eval(fam0$initialize)
    ## to be completed
    if (EMPTY)
        stop("invalid linear predictor values in empty model")

    qrx <- qr(x)
    rank <- qrx$rank
    is_full_rank <- rank == nvars
    if (!singular.ok && !is_full_rank) {
        stop("singular fit encountered")
    }

    if (!isTRUE(is_full_rank)) {
        aliased <- qrx$pivot[seq.int(qrx$rank + 1,nvars)]
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
    betas_all <- structure(rep(NA_real_, nvars_all), .Names = betas_names_all)
    keep <- weights > 0
    nkeep <- sum(keep)
    df_residual <- nkeep - rank
    warn <- getOption("warn")
    options(warn = -1)

    if (is.null(start)) {
        fit<- glm.fit(x = x, y = y, weights = weights,
                      start = start, offset = offset, family =fam0,
                      control = list(epsilon = control$epsilon, maxit = 100,
                                     trace = FALSE), intercept =  intercept)
        dispersion <- linkfun_disp(1/as.vector(MASS::theta.ml(y = y, mu = fitted(fit), n = nobs, weights = weights,
                                                              trace = control$trace > 2)))
        options(warn = warn)
        betas <- coef(fit)
        names(betas) <- betas_names
    } else {
        if ((length(start) == nvars_all) & is.numeric(start) ) {
            betas_all <- start
            names(betas_all) <- betas_names_all
            if (!isTRUE(is_full_rank)) {
                betas_all[aliased] <- NA_real_
                betas <- betas_all[-aliased]
            } else {
                betas <- betas_all
            }
            etas <- drop(x %*% betas + offset)
            dispersion <- linkfun_disp(1/as.vector(MASS::theta.ml(y = y, mu = linkinv(etas), n = nobs, weights = weights,
                                                                  trace = control$trace > 2)))
        }

        if ((length(start) == nvars_all+1) & is.numeric(start) ) {
            betas_all <- start[seq.int(nvars_all)]
            names(betas_all) <- betas_names_all
            if (!isTRUE(is_full_rank)) {
                betas_all[aliased] <- NA_real_
                betas <- betas_all[-aliased]
            } else {
                betas <- betas_all
            }
            dispersion <- start[nvars_all + 1]
        }

        if (length(start) > nvars_all + 1 | length(start) < nvars_all ) {
            stop(paste(paste(gettextf("length of 'start' should be equal to %d and correspond to initial betas for %s",
                                      nvars_all, paste(deparse(betas_names_all),
                                                       collapse = ", "), "or", gettextf("to %d and also include a starting value for the transformed dispersion",
                                                                                        nvars_all)))), domain = NA_real_)
        }
    }

    adjusted_grad_all <- rep(NA_real_, nvars_all + 1)
    names(adjusted_grad_all) <- c(betas_names_all, "dispersion")
    names(dispersion) <- "dispersion"
    par <- c(betas, dispersion)
    quantities <- key_quantities(par)
    ## keeped to test
                                        # if (control$type == "MPL_Jeffreys"){
                                        #   opt <- try(nlminb(par, loglikMPL,fit = NULL,lower = c(rep(-Inf,nvars),10^{-8}), upper = rep(+Inf, nvars + 1) ),TRUE)
                                        #   #betas <-opt$par[-(nvars + 1)]
                                        #   #dispersion <- opt$par[nvars + 1]
                                        #   cat("\n","Estimates par. with Jeff obtained with nlminb as test: ","\n\n")
                                        #   print(opt)
                                        # }
    step_components_beta <- compute_step_components(par, level = 0, fit = quantities)
    step_components_dispersion <- compute_step_components(par, level = 1, fit = quantities)
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
    if (step_components_dispersion$failed_inversion)  {
        warning("failed to invert the information matrix")
    }
    if (step_components_dispersion$failed_adjustment) {
        warning("failed to calculate score adjustment")
    }
    adjusted_grad_dispersion <- with(step_components_dispersion, {
        grad + adjustment
    })
    step_dispersion <- as.vector(adjusted_grad_dispersion * step_components_dispersion$inverse_info)

    if (control$maxit == 0) {
        iter <- 0
        failed <- FALSE
    } else {
        for (iter in seq.int(control$maxit)) {
            step_factor <- 0
            testhalf <- TRUE
            while (testhalf & step_factor < control$max_step_factor) {
                step_beta_previous <- step_beta
                step_dispersion_previous <- step_dispersion
                betas <- betas + control$slowit * 2^(-step_factor) * step_beta
                dispersion <- dispersion + 2^(-step_factor) * step_dispersion
                par <- c(betas, dispersion)
                quantities <- key_quantities(par)
                step_components_beta <- compute_step_components(par, level = 0, fit = quantities)
                step_components_dispersion <- compute_step_components(par,  level = 1, fit = quantities)
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
                step_beta <- drop(step_components_beta$inverse_info %*%  adjusted_grad_beta)

                if (failed_inversion_dispersion <- step_components_dispersion$failed_inversion) {
                    warning("failed to invert the information matrix")
                    break
                }
                if (failed_adjustment_dispersion<- step_components_dispersion$failed_adjustment) {
                    warning("failed to calculate score adjustment")
                    break
                }
                adjusted_grad_dispersion <- with(step_components_dispersion, {
                    grad + adjustment
                })
                step_dispersion <- as.vector(adjusted_grad_dispersion * step_components_dispersion$inverse_info)

                if (step_factor == 0 & iter == 1) {
                    testhalf <- TRUE
                } else {
                    s2 <- c(abs(step_beta), abs(step_dispersion))
                    s1 <- c(abs(step_beta_previous), abs(step_dispersion_previous))
                    testhalf <- sum(s2, na.rm = TRUE) > sum(s1,  na.rm = TRUE)
                }
                step_factor <- step_factor + 1
                if (control$trace) {
                    trace_iteration()
                }
            }
            failed <- failed_adjustment_beta | failed_inversion_beta |
                failed_adjustment_dispersion| failed_inversion_dispersion
            if (failed | sum(abs(c(step_beta, step_dispersion)),na.rm = TRUE) < control$epsilon) {
                break
            }
        }
    }
    adjusted_grad_all[betas_names] <- adjusted_grad_beta
    adjusted_grad_all["dispersion"] <- adjusted_grad_dispersion
    betas_all[betas_names] <- betas

    if(iter >= control$maxit) {
        convergence <- FALSE
        warning("optimization failed to converge")
    } else
    {
        convergence <- TRUE
    }
    if (!isTRUE(is_full_rank)) {
        x <- X_all
        betas <- betas_all
        betas[is.na(betas)] <- 0
        nvars <- nvars_all
    }

    par <- c(betas, dispersion)
    if (is_correction)
    {
        bias <- c(- step_components_beta$inverse_info%*%AS_mean_adjustment(par, level = 0, fit = quantities),
                  - step_components_dispersion$inverse_info%*%AS_mean_adjustment(par, level = 1, fit = quantities))
        par <- par - bias
        betas <- par[-(nvars + 1)]
        dispersion <-par[(nvars + 1)]
    }
    quantities <- key_quantities(par)
    step_components_beta <- compute_step_components(par, level = 0, fit = quantities)
    step_components_dispersion<- compute_step_components(par, level = 1, fit = quantities)

    qr.Wx <- quantities$qr_decomposition
    mus <- quantities$mus
    etas <- quantities$etas
    residuals <- with(quantities, (y - mus) / d1)
    working_weights <- quantities$working_weights
    wt <- rep.int(0, nobs)
    wt[keep] <- working_weights[keep]
    names(wt) <- names(residuals) <- names(mus) <- names(etas) <- names(weights) <- names(y) <- ynames

    fam <- do.call("negative.binomial",list(theta = 1 / linkinv_disp(dispersion),
                                            link = link))

    if (attr(Terms, "intercept") & missing_offset) {
        nullFit <- glm.fit(x[, "(Intercept)", drop = FALSE], y, weights,
                           offset = rep(0, nobs), family = fam, control = list(maxit = control$maxit,
                                                                               epsilon = control$epsilon, trace = control$trace > 1), intercept = TRUE)
        nullmus <- nullFit$fitted
    }
    if (!attr(Terms, "intercept")) {
        nullmus <- linkinv(offset)
    }
    if (attr(Terms, "intercept") & !missing_offset) {
        nullFit <- glm.fit(x[, "(Intercept)", drop = FALSE], y, weights,
                           offset = offset, family = fam, control = list(maxit = control$maxit,
                                                                         epsilon = control$epsilon, trace = control$trace > 1), intercept = TRUE)
        nullmus <- nullFit$fitted
    }

    th <- 1 / linkinv_disp(dispersion)
    Lm <- loglik(nobs, th, mus, y, weights)
    dev.resids <- fam$dev.resids
    aic <- fam$aic
    nulldev <- sum(dev.resids(y, nullmus, weights))
    nulldf <- nkeep - as.integer(attr(Terms, "intercept"))
    deviance <- sum(dev.resids(y, mus, weights))
    aic.model <- aic(y, nobs, mus, weights, deviance) + 2 * (rank + 1)

    vcov.mean <-step_components_beta$inverse_info
    rownames(vcov.mean) <- colnames(vcov.mean)  <-  betas_names

    vcov.dispersion <- step_components_dispersion$inverse_info
    names(vcov.dispersion) <-  paste0(control$transformation, "(dispersion)")
    out <- list(coefficients = betas,
                vcov.mean = vcov.mean,
                vcov.dispersion =  vcov.dispersion,
                residuals = residuals,
                fitted.values = mus,
                R = if (!EMPTY) qr.R(qr.Wx),
                rank = rank,
                qr = if (!EMPTY) structure(qr.Wx[c("qr", "rank", "qraux", "pivot", "tol")], class = "qr"),
                family = fam,
                twologlik = as.vector(2 * Lm),
                linear.predictors = etas,
                deviance = deviance,
                aic = aic.model,
                link=link,
                null.deviance = nulldev,
                iter = iter,
                weights = wt,
                prior.weights = weights,
                df.residual = df_residual,
                df.null = nulldf,
                y = y,
                converged = convergence,
                dispersion = dispersion,
                grad = adjusted_grad_all,
                transformation = control$transformation,
                theta = as.vector(th), ## dispersion as in glm.nb
                type = control$type)
    out$terms <- Terms
    out$formula <- as.vector(attr(Terms,"formula"))
    Call$link <- link
    out$call <- Call
    if (model) out$model <- mf
    out$na.action <- attr(mf,"na.action")
    if (return_x) out$x <- x
    if (!return_y) out$y <- NULL
    out$contrasts <- attr(x,"contrasts")
    out$xlevels <- .getXlevels(Terms,mf)
    out$control <- control
    out$offset <- offset
    class(out)<-c("brnb","negbin","glm")
    out
}


#' Extract model coefficients from [`"brnb"`][brnb] objects
#'
#' @inheritParams stats::coef
#' @param model one of `"mean"` (default), `"full"`, `"dispersion"`,
#'     to return the estimates of the parameters in the linear
#'     prediction only, or both, the estimate of the dispersion
#'     parameter only, respectively.
#'
#' @details
#'
#' See [coef()] for more details.
#'
#' @method coef brnb
#' @export
coef.brnb <- function(object, model = c("mean", "full", "dispersion"), ...) {
    object$transformed_dispersion <- object$dispersion
    coef.brglmFit(object, model, ...)
}

#' Extract model variance-covariance matrix from [`"brnb"`][brnb] objects
#'
#'
#' @inheritParams stats::vcov.glm
#' @param object an object of class [`"brnb"`][brnb], typically, a result of a call to [brnb()].
#' @param model character specifying for which component of the model variance-covariance matrix should be extracted.
#'
#' @details
#'
#' The options for `model` are `"mean"` for mean regression only
#'  (default), `"dispersion"` for the dispersion parameter (in a
#'  chosen transformation; see [brglmControl()], and `"full"` for both
#'  the mean regression and the (transformed) dispersion parameters.
#'  See [vcov()] for more details.
#'
#' @seealso [vcov()]
#'
#' @method vcov brnb
#' @export
vcov.brnb <- function(object, model = c("mean", "full", "dispersion"), complete = TRUE, ...) {
    object$info_transformed_dispersion <- 1/object$vcov.dispersion
    object$dispersion <- 1
    vcov.brglmFit(object, model , complete , ...)
}

#' [summary()] method for [`"brnb"`][brnb] objects
#'
#' @inheritParams stats::summary.glm
#' @param object an object of class [`"brnb"`][brnb], typically, a
#'     result of a call to [brnb()].
#' @param x an object of class [`"summary.brnb"`][summary.brnb],
#'     usually, a result of a call to [summary.brnb].
#' @details The interface of the summary method for [`"brnb"`][brnb]
#'     objects is similar to that of [`"brglmFit"`][brglmFit] objects
#'     with additional information.
#'
#'     p-values of the individual Wald statistics are based on the
#'     standard normal distribution.
#'
#' @seealso [summary.brglmFit()] and [glm()]
#'
#' @examples
#' # For examples see examples(brnb)
#'
#' @method summary brnb
#' @export
summary.brnb <- function(object, ...) {
    ## coefficient table
    cf <- object$coefficients
    cf.dispersion <- object$dispersion
    se.mean <- sqrt(diag(object$vcov.mean))
    se.dispersion <- sqrt(object$vcov.dispersion)
    cf <- cbind(cf, se.mean, cf/se.mean, 2 * pnorm(-abs(cf/se.mean)))
    colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(cf) <- names(object$coefficients)
    object$coefficients <- cf
    object$se.betas <- se.mean
    object$se.dispersion <- se.dispersion
    coef.dispersion <- cbind(cf.dispersion, se.dispersion)
    colnames(coef.dispersion) <- c("Estimate", "Std. Error")
    rownames(coef.dispersion) <- c("dispersion")
    object$coef.dispersion <- coef.dispersion
    ## number of iterations
    object$iterations <- object$iter
    ## delete some slots
    object$terms <- object$model <- object$y <-NULL
    object$x <- object$levels <- object$contrasts <-  NULL
    ## return
    class(object) <- "summary.brnb"
    object
}

#' @rdname summary.brnb
#' @method print summary.brnb
#' @export
print.summary.brnb <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
    if (!x$converged) {
        cat("model did not converge\n")
    }
    if (NROW(x$coefficients)) {
        cat(paste("\nCoefficients (mean model with ", x$link, " link):\n", sep = ""))
        printCoefmat(x$coefficients, digits = digits, signif.legend = FALSE)
    } else {
        cat("\nNo coefficients (in mean model)\n")
    }
    if (getOption("show.signif.stars") & any( x$coefficients[, 4L] < 0.1)) {
        cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
    }
    if (NROW(x$coef.dispersion)) {
        cat(paste("\n(Dispersion parameter  with ", x$transformation," transformation function: ",
                  round(x$coef.dispersion[,1], 7),")\n", sep = ""))
    } else {
        cat("\nNo coefficients (in precision model)\n")
    }
    cat("\n", apply(cbind(paste(format(c("Null",
                                         "Residual"), justify = "right"), "deviance:"),
                          format(unlist(x[c("null.deviance",
                                            "deviance")]), digits = max(5, digits + 1)), " on",
                          format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"),
                    1, paste, collapse = " "), sep = "")
    cat(paste("AIC:", round(x$aic,digits)))
    cat("\n\nType of estimator:", x$type, get_type_description(x$type))
    cat("\n", "Number of quasi-Fisher scoring iterations:", x$iter, "\n", sep = "")
    invisible(x)
}


#' @method print brnb
#' @export
print.brnb <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
    if (!x$converged) {
        cat("model did not converge\n")
    } else {
        if (length(x$coefficients)) {
            cat(paste("Coefficients (mean model with ", x$link, " link):\n", sep = ""))
            print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
            cat("\n")
        } else {
            cat("No coefficients (in mean model)\n\n")
        }
        cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
            x$df.residual, "Residual\n")
        if (nchar(mess <- naprint(x$na.action)))
            cat("  (", mess, ")\n", sep = "")
        cat("Null Deviance:\t", format(round(x$null.deviance, digits)), "\n")
        cat("Residual Deviance:\t", format(round(x$deviance, digits)),
            "\tAIC:", format(round(x$aic, digits)), "\n")
    }
}

#' Method for computing Wald confidence intervals for one or more
#' regression parameters in a [`"brnb"`][brnb] object
#'
#' @inheritParams stats::confint
#'
#' @method confint brnb
#' @export
confint.brnb <- function(object, parm, level = 0.95, ...) {
    confint.default(object, parm, level, ...)
}

#' Simulate Responses
#'
#' Simulate one or more responses from the distribution corresponding
#' to a fitted model [`"brnb"`][brnb] object.
#' @param object an object representing a fitted model.
#' @param nsim number of response vectors to simulate. Defaults to 1.
#' @param seed an object specifying if and how the random number
#'     generator should be initialized; see [set.seed()] for details.
#' @param ... extra arguments to be passed to methods. Not currently
#'     used.

#' @examples
#' # Example in Saha, K., & Paul, S. (2005). Bias-corrected maximum
#' # likelihood estimator of the negative binomial dispersion
#' # parameter.  Biometrics, 61, 179--185.
#' #
#' # Frequency distribution of red mites on apple leaves.
#' nomites <- 0:8
#' noleaves <- c(70, 38, 17, 10, 9, 3, 2, 1, 0)
#' fit_glmnb <- MASS::glm.nb(nomites~1,link="identity",weights = noleaves)
#' fit_brnb <- brnb(nomites ~ 1, link = "identity", transformation = "inverse",
#'                  type = "ML",weights = noleaves)
#' ## Let us simulate 10 response vectors
#' sim_glmnb <- simulate(fit_glmnb, nsim = 10, seed = 123)
#' sim_brnb <-  simulate(fit_brnb, nsim = 10, seed = 123)
#' # The results from glm.nb and brnb with type = "ML" are
#' # exactly the same
#' all.equal(sim_glmnb, sim_brnb, check.attributes = FALSE)


#' @method simulate brnb
#' @export
simulate.brnb <- function(object, nsim = 1, seed = NULL, ...) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    theta <- object$theta
    mus <- fitted(object)
    n <- length(mus)
    val <- rnegbin(n*nsim, mu = mus, theta = theta)
    dim(val) <- c(n, nsim)
    val <- as.data.frame(val)
    names(val) <- paste0("sim_", seq_len(nsim))
    attr(val, "seed") <- RNGstate
    val
}


