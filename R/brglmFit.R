# Copyright (C) 2016 - Ioannis Kosmidis, Oliver Clark

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


#' Fitting function for [glm()] for reduced-bias estimation and
#' inference
#'
#' [brglmFit()] is a fitting method for [glm()] that fits generalized
#' linear models using implicit and explicit bias reduction methods
#' (Kosmidis, 2014), and other penalized maximum likelihood
#' methods. Currently supported methods include the mean bias-reducing
#' adjusted scores approach in Firth (1993) and Kosmidis & Firth
#' (2009), the median bias-reduction adjusted scores approach in Kenne
#' Pagui et al. (2017), the correction of the asymptotic bias in
#' Cordeiro & McCullagh (1991), the mixed bias-reduction adjusted
#' scores approach in Kosmidis et al (2020), maximum penalized
#' likelihood with powers of the Jeffreys prior as penalty, and
#' maximum likelihood.
#'
#' @inheritParams stats::glm.fit
#' @aliases brglm_fit
#' @param x a design matrix of dimension `n * p`.
#' @param y a vector of observations of length `n`.
#' @param control a list of parameters controlling the fitting
#'     process. See [brglmControl()] for details.
#' @param start starting values for the parameters in the linear
#'     predictor. If `NULL` (default) then the maximum likelihood
#'     estimates are calculated and used as starting values.
#' @param mustart applied only when start is not `NULL`. Starting
#'     values for the vector of means to be passed to
#'     [glm.fit()] when computing starting values using maximum
#'     likelihood.
#' @param etastart applied only when start is not `NULL`. Starting
#'     values for the linear predictor to be passed to
#'     [glm.fit()] when computing starting values using maximum
#'     likelihood.
#' @param fixed_totals effective only when `family` is
#'     [poisson()]. Either `NULL` (no effect) or a vector that
#'     indicates which counts must be treated as a group. See Details
#'     for more information and [brmultinom()].
#' @param singular.ok logical. If `FALSE`, a singular model is an
#'     error.
#'
#' @details
#'
#' Trust-region fitting method for bias-reduced GLMs
#' Fits bias-reduced GLMs by trust-region optimization of the adjusted score equations.
#' Objective: f(beta) = ||adjusted_score_beta(beta)||^2 / 2.
#'
#' In the special case of generalized linear models for binomial,
#' Poisson and multinomial responses, the adjusted score equation
#' approaches for `type = "AS_mixed"`, `type = "AS_mean"`, and `type =
#' "AS_median"` (see below for what methods each `type` corresponds)
#' return estimates with improved frequentist properties, that are
#' also always finite, even in cases where the maximum likelihood
#' estimates are infinite (e.g. complete and quasi-complete separation
#' in multinomial regression). See, Kosmidis and Firth (2021) for a
#' proof for binomial-response GLMs with Jeffreys-prior penalties to
#' the log-likelihood, which is equivalent to mean bias reduction for
#' logistic regression. See, also,
#' [detectseparation::detect_separation()] and
#' [detectseparation::check_infinite_estimates()] for pre-fit and
#' post-fit methods for the detection of infinite estimates in
#' binomial response generalized linear models.
#'
#' The type of score adjustment to be used is specified through the
#' `type` argument (see [brglmControl()] for details). The available
#' options are
#'
#' * `type = "AS_mixed"`: the mixed bias-reducing score adjustments in
#' Kosmidis et al (2020) that result in mean bias reduction for the
#' regression parameters and median bias reduction for the dispersion
#' parameter, if any; default.
#'
#' * `type = "AS_mean"`: the mean bias-reducing score adjustments in
#' Firth, 1993 and Kosmidis & Firth, 2009. `type = "AS_mixed"` and
#' `type = "AS_mean"` will return the same results when `family` is
#' [binomial()] or [poisson()], i.e. when the dispersion is fixed
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
#' The null deviance is evaluated based on the fitted values using the
#' method specified by the `type` argument (see [brglmControl()]).
#'
#' The `family` argument of the current version of [brglmFit()] can
#' accept any combination of [`"family"`][family] objects and link functions,
#' including families with user-specified link functions, [mis()]
#' links, and [power()] links, but excluding [quasi()],
#' [quasipoisson()] and [quasibinomial()] families.
#'
#' The description of `method` argument and the `Fitting functions`
#' section in [glm()] gives information on supplying fitting
#' methods to [glm()].
#'
#' `fixed_totals` specifies groups of observations for which the sum
#' of the means of a Poisson model will be held fixed to the observed
#' count for each group. This argument is used internally in
#' [brmultinom()] and [bracl()] for baseline-category logit models and
#' adjacent category logit models, respectively.
#'
#' [brglm_fit()] is an alias to [brglmFit()].
#' 
#' @return
#'
#' An object inheriting from [`"brglmFit"`][brglmFit()] object, which
#' is a list having the same elements to the list that
#' [stats::glm.fit()] returns, with a few extra arguments.
#'
#' @author Ioannis Kosmidis `[aut, cre]` \email{ioannis.kosmidis@warwick.ac.uk}, Euloge Clovis Kenne Pagui `[ctb]` \email{kenne@stat.unipd.it}
#'
#' @seealso [brglmControl()], [glm.fit()], [glm()]
#'
#' @references
#'
#' Kosmidis I, Firth D (2021). Jeffreys-prior penalty, finiteness
#' and shrinkage in binomial-response generalized linear
#' models. *Biometrika*, **108**, 71-82. \doi{10.1093/biomet/asaa052}.
#'
#' Kosmidis I, Kenne Pagui E C, Sartori N (2020). Mean and median bias
#' reduction in generalized linear models. *Statistics and Computing*,
#' **30**, 43-59. \doi{10.1007/s11222-019-09860-6}.
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
#' Kosmidis I, Firth D (2009). Bias reduction in exponential family
#' nonlinear models. *Biometrika*, **96**, 793-804. \doi{10.1093/biomet/asp055}.
#'
#' Kosmidis I, Firth D (2010). A generic algorithm for reducing
#' bias in parametric estimation. *Electronic Journal of Statistics*,
#' **4**, 1097-1112. \doi{10.1214/10-EJS579}.
#'
#' Kosmidis I (2014). Bias in parametric estimation: reduction and
#' useful side-effects. *WIRE Computational Statistics*, **6**,
#' 185-196. \doi{10.1002/wics.1296}.
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
#' \donttest{
#' ## Another example from
#' ## King, Gary, James E. Alt, Nancy Elizabeth Burns and Michael Laver
#' ## (1990).  "A Unified Model of Cabinet Dissolution in Parliamentary
#' ## Democracies", _American Journal of Political Science_, **34**, 846-870
#'
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
#' \donttest{
#' ## An example with offsets from Venables & Ripley (2002, p.189)
#' data("anorexia", package = "MASS")
#'
#' anorexML <- glm(Postwt ~ Prewt + Treat + offset(Prewt),
#'                 family = gaussian, data = anorexia)
#' anorexBC <- update(anorexML, method = "brglmFit", type = "correction")
#' anorexBR_mean <- update(anorexML, method = "brglmFit")
#' anorexBR_median <- update(anorexML, method = "brglmFit", type = "AS_median")
#'
#' # All methods return the same estimates for the regression
#' # parameters because the maximum likelihood estimator is normally
#' # distributed around the `true` value under the model (hence, both
#' # mean and component-wise median unbiased). The Wald tests for
#' # anorexBC and anorexBR_mean differ from anorexML because the
#' # bias-reduced estimator of the dispersion is the unbiased, by
#' # degree of freedom adjustment (divide by n - p), estimator of the
#' # residual variance. The Wald tests from anorexBR_median are based
#' # on the median bias-reduced estimator of the dispersion that
#' # results from a different adjustment of the degrees of freedom
#' # (divide by n - p - 2/3)
#' summary(anorexML)
#' summary(anorexBC)
#' summary(anorexBR_mean)
#' summary(anorexBR_median)
#' }
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
#' TODO: Vignette write,  

#' @export
brglmFit <- function(x, y, weights = rep(1, nobs),
                     start = NULL, etastart = NULL,
                     mustart = NULL, offset = rep(0, nobs),
                     family = gaussian(),
                     control = list(), intercept = TRUE,
                     fixed_totals = NULL, singular.ok = TRUE) {

    ## compute_step_components does everything on the scale of the /transformed/ dispersion
    compute_step_components <- function(pars, fit, level = 0) {
        if (level == 0) {
            # Beta components
            grad <- fit$grad_beta
            inverse_info <- fit$inverse_info_beta
            adjustment <- adjustment_function(pars, fit = fit, level = 0,
                                                     x, nobs, nvars, weights)
            failed_adjustment <- any(is.na(adjustment))
            failed_inversion <- FALSE # Already computed successfully
        } else {
            # Dispersion components
            if (fit$no_dispersion || df_residual < 1) {
                grad <- adjustment <- inverse_info <- NA_real_
                failed_adjustment <- failed_inversion <- FALSE
            } else {
                d1zeta <- eval(d1_transformed_dispersion)
                d2zeta <- eval(d2_transformed_dispersion)
                grad <- fit$grad_zeta / d1zeta
                inverse_info <- fit$inverse_info_zeta * d1zeta^2
                adjustment <- adjustment_function(pars, fit = fit, level = 1, x, nobs, 
                                                  nvars, weights) / d1zeta - 0.5 * d2zeta / d1zeta^2
                failed_inversion  <- !is.finite(inverse_info)
                failed_adjustment <- is.na(adjustment)
            }
        }
        list(grad = grad, inverse_info = inverse_info, adjustment = adjustment,
             failed_adjustment = failed_adjustment, failed_inversion = failed_inversion)
    }




    customTransformation <- is.list(control$transformation) & length(control$transformation) == 2
    if (customTransformation) transformation0 <- control$transformation

    control <- do.call("brglmControl", control)

    adjustment_function <- switch(control$type,
        "correction"   = AS_mean_adjustment,
        "AS_mean"      = AS_mean_adjustment,
        "AS_median"    = AS_median_adjustment,
        "AS_mixed"     = AS_mixed_adjustment,
        "MPL_Jeffreys" = AS_Jeffreys_adjustment,
        "ML"           = function(pars, ...) 0)

    ## Some useful quantities
    is_ML         <- control$type == "ML"
    is_AS_median  <- control$type == "AS_median"
    is_AS_mixed   <- control$type == "AS_mixed"
    is_correction <- control$type == "correction"
    no_dispersion <- family$family %in% c("poisson", "binomial")

    if (is_ML | is_AS_median | is_AS_mixed) {
        transformation1 <- control$transformation
        Trans1 <- control$Trans
        inverseTrans1 <- control$inverseTrans
        ## Set the transformation to identity
        control$transformation <- "identity"
        control$Trans <- expression(dispersion)
        control$inverseTrans <- expression(transformed_dispersion)
    }

    ## If fixed_totals is specified the compute row_totals
    row_totals <- NULL
    if (is.null(fixed_totals)) {
        has_fixed_totals <- FALSE
    } else {
        if (family$family == "poisson") {
            row_totals  <- as.vector(tapply(y, fixed_totals, sum))[fixed_totals]
            has_fixed_totals <- TRUE
        } else {
            has_fixed_totals <- FALSE
        }
    }

    ## Ensure x is a matrix, extract variable names, observation
    ## names, nobs, nvars, and initialize weights and offsets if
    ## needed

    x <- as.matrix(x)
    betas_names <- dimnames(x)[[2L]]
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(betas_names) & !EMPTY)
        betas_names <- colnames(x) <- paste0("x", seq.int(nvars))
    ynames <- if (is.matrix(y)) rownames(y) else names(y)
    converged <- FALSE
    nobs   <- NROW(y)
    if (is.null(weights))              weights <- rep.int(1, nobs)
    if (missing_offset <- is.null(offset)) offset  <- rep.int(0, nobs)

    ok_links <- c("logit", "probit", "cauchit", "cloglog",
                  "identity", "log", "sqrt", "inverse")

    if (isTRUE(family$family %in% c("quasi", "quasibinomial", "quasipoisson")))
        stop("`brglmFit` does not currently support quasi families.")

    ## Enrich family
    family <- enrichwith::enrich(family, with = c("d1afun", "d2afun", "d3afun", "d1variance"))
    if ((family$link %in% ok_links) | grepl("mu\\^", family$link)) {
        ## Enrich the link object with d2mu.deta and update family object
        linkglm <- enrichwith::enrich(make.link(family$link), with = "d2mu.deta")
        ## Put everything into the family object
        family[names(linkglm)] <- linkglm
    }

    ## Annoying thing is that link-glm components other than the
    ## standard ones disappear when extra arguments are passed to a
    ## family functions... Anyway, we only require d2mu.deta here.

    ## Extract functions from the enriched family object

    variance <- family$variance
    d1variance <- family$d1variance
    linkinv <- family$linkinv
    linkfun <- family$linkfun
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object", call. = FALSE)
    mu.eta  <- family$mu.eta
    dev.resids <- family$dev.resids
    aic  <- family$aic

    ## If the family is custom then d2mu.deta cannot survive when
    ## passing throguh current family functions. But mu.eta does; so
    ## we compute d2mu.deta numerically; this allows also generality,
    ## as the users can then keep their custom link implementations
    ## unaltered. Issue is scalability, due to the need of evaluating
    ## n numerical derivatives
    if (is.null(family$d2mu.deta))
        family$d2mu.deta <- function(eta) numDeriv::grad(mu.eta, eta)

    d1_transformed_dispersion <- DD(control$Trans, "dispersion", order = 1)
    d2_transformed_dispersion <- DD(control$Trans, "dispersion", order = 2)

    ## Check for invalid etas and mus
    valid_eta <- unless_null(family$valideta, function(eta) TRUE)
    valid_mu <- unless_null(family$validmu,  function(mu)  TRUE)

    mustart <- NULL; etastart <- NULL

    ## Initialize as prescribed in family
    eval(family$initialize)

    ## If there are no covariates in the model then evaluate only the offset
    if (EMPTY) {
        etas <- rep.int(0, nobs) + offset
        if (!valid_eta(etas)) stop("invalid linear predictor values in empty model", call. = FALSE)
        mus  <- linkinv(etas)
        if (!valid_mu(mus))   stop("invalid fitted means in empty model", call. = FALSE)
        working_weights <- ((weights * mu.eta(etas)^2) / variance(mus))^0.5
        residuals <- (y - mus) / mu.eta(etas)
        boundary <- converged <- TRUE
        betas_all <- numeric()
        rank <- 0
        iter <- 0L
        keep  <- weights > 0
        nkeep <- sum(keep)
        df_residual <- nkeep
    } else {

        boundary <- converged <- FALSE
        ## Detect aliasing
        if (!isTRUE(control$check_aliasing)) {
            is_full_rank <- TRUE ## Assumption
            rank <- nvars_all <- nvars
            betas_names_all <- betas_names
        } else {
            qrx <- qr(x)
            rank <- qrx$rank
            is_full_rank <- rank == nvars
            if (!isTRUE(singular.ok) && !isTRUE(is_full_rank))
                stop("singular fit encountered")
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

        ## Handle starting values
        ## If start is NULL then start at the ML estimator else use start
        if (is.null(start)) {
            ## Adjust counts if binomial or Poisson in order to avoid infinite estimates
            adj <- control$response_adjustment
            if (is.null(adj)) adj <- nvars / nobs
            if (family$family == "binomial") {
                weights.adj <- weights + (!(is_correction)) * adj
                y.adj       <- (weights * y + (!(is_correction)) * 0.5 * adj) / weights.adj
            } else {
                weights.adj <- weights
                y.adj       <- y + if (family$family == "poisson") (!(is_correction)) * 0.5 * adj else 0
            }
            ## ML fit to get starting values
            ## Get startng values and kill warnings whilst doing that
            suppressWarnings(
                tempFit <- glm.fit(x = x, y = y.adj, weights = weights.adj,
                                   etastart = etastart, mustart = mustart,
                                   offset = offset, family = family,
                                   control = list(epsilon = control$epsilon,
                                                  maxit = 10000, trace = FALSE),
                                   intercept = intercept)
            )
            betas <- coef(tempFit); names(betas) <- betas_names
            dispList <- estimate_dispersion(betas = betas, y = y, x = x,
                                            weights = weights, offset = offset,
                                            family = family,
                                            fixed_totals = fixed_totals,
                                            row_totals = row_totals,
                                            no_dispersion = no_dispersion,
                                            nobs = nobs, nvars = nvars, keep = keep,
                                            df_residual = df_residual, control = control)
            dispersion <- dispList$dispersion
            if (is.na(dispersion))
                dispersion <- var(y) / variance(sum(weights * y) / sum(weights))
            # Enforce positivity constraint
            if (!is.na(dispersion) && dispersion <= 0) {
                warning(sprintf("Dispersion became non-positive (%.6g); resetting to small positive value.", dispersion))
                dispersion <- .Machine$double.eps
            }
            dispersion_ML <- dispList$dispersion_ML
            transformed_dispersion <- eval(control$Trans)
        } else {
            if ((length(start) == nvars_all) & is.numeric(start)) {
                betas_all <- start; names(betas_all) <- betas_names_all
                if (!isTRUE(is_full_rank)) {
                    betas_all[aliased] <- NA_real_
                    betas <- betas_all[-aliased]
                } else {
                    betas <- betas_all
                }
                ## Estimate dispersion based on current value for betas
                dispList <- estimate_dispersion(betas = betas, y = y, x = x,
                                                weights = weights, offset = offset,
                                                family = family,
                                                fixed_totals = fixed_totals,
                                                row_totals = row_totals,
                                                no_dispersion = no_dispersion,
                                                nobs = nobs, nvars = nvars, keep = keep,
                                                df_residual = df_residual, control = control)
                dispersion <- dispList$dispersion
                if (is.na(dispersion))
                    dispersion <- var(y) / variance(sum(weights * y) / sum(weights))
                # Enforce positivity constraint
                if (!is.na(dispersion) && dispersion <= 0) {
                    warning(sprintf("Dispersion became non-positive (%.6g); resetting to small positive value.", dispersion))
                    dispersion <- .Machine$double.eps
                }
                dispersion_ML          <- dispList$dispersion_ML
                transformed_dispersion <- eval(control$Trans)
            } else if ((length(start) == nvars_all + 1) & is.numeric(start)) {
                betas_all <- start[seq.int(nvars_all)]; names(betas_all) <- betas_names_all
                if (!isTRUE(is_full_rank)) {
                    betas_all[aliased] <- NA_real_
                    betas <- betas_all[-aliased]
                } else {
                    betas <- betas_all
                }
                transformed_dispersion <- start[nvars_all + 1]
                dispersion_ML <- NA_real_
                dispersion <- eval(control$inverseTrans)
                # Enforce positivity constraint
                if (!is.na(dispersion) && dispersion <= 0) {
                    warning(sprintf("Dispersion became non-positive (%.6g); resetting to small positive value.", dispersion))
                    dispersion <- .Machine$double.eps
                }
            } else {
                stop(gettextf(
                    "length of 'start' should equal %d (betas) or %d (betas + dispersion)",
                    nvars_all, nvars_all + 1), domain = NA)
            }
        }

        adjusted_grad_all <- rep(NA_real_, nvars_all + 1)
        names(adjusted_grad_all) <- c(betas_names_all, "Transformed dispersion")

        if (is_correction) {
            ## Needs original fisher scoring implementation to work for correction as this is specific to a single fisher scoring step
            if (control$maxit > 0) control$maxit <- 1
            control$slowit <- 1
            control$max_step_factor <- 1
        }

        theta                  <- c(betas, dispersion)
        transformed_dispersion <- eval(control$Trans)

        # Determine if we need to compute the inverse of the information matrix for the dispersion parameter
        needs_inverse <- !no_dispersion || control$type %in% c("AS_median", "AS_mixed")

        # Avoid unconstrained trust-region steps for dispersion in null/intercept-only models
        # Only update dispersion if model is not empty/null
        if (!EMPTY) {
        fit <- compute_fit(pars = theta, y = y, x = x, weights = weights,
                           offset = offset, family = family,
                           fixed_totals = fixed_totals, row_totals = row_totals,
                           no_dispersion = no_dispersion, nobs = nobs, nvars = nvars,
                           keep = keep, need_qr = TRUE, need_hatvalues = TRUE,
                           need_inverse = needs_inverse)
        } else {
            # For null/intercept-only models, use safe default for dispersion
            dispersion <- 1
            theta <- c(betas, dispersion)
            transformed_dispersion <- eval(control$Trans)
            fit <- NULL
        }

        step_components_beta <- compute_step_components(theta, level = 0, fit = fit)
        step_components_zeta <- compute_step_components(theta, level = 1, fit = fit)
        if (step_components_beta$failed_inversion)  warning("failed to invert the information matrix")
        if (step_components_beta$failed_adjustment) warning("failed to calculate score adjustment")

        adjusted_grad_beta <- with(step_components_beta, grad + adjustment)
        adjusted_grad_zeta <- if (no_dispersion) NA_real_ else
                              with(step_components_zeta, grad + adjustment)

        
        # trust region parameters (Nocedal & Wright 2006)
        Delta         <- 1.0
        Delta_max     <- 1e10
        radius_shrink <- 0.25
        radius_expand <- 0.75
        shrink_factor <- 0.25
        expand_factor <- 2.0
        cg_maxiter    <- min(nvars, 50)
        eta           <- 0.1 ## Add to control?

        # Attempt infrequent hatvalue updates
        hatvalues_cached <- fit$hatvalues
        hat_accept_count <- 0L
        hat_update_freq  <- control$hat_update_freq
        #hat_update_freq <- max(3L, min(10L, as.integer(floor(nobs / (5 * nvars)))))

        failed <- FALSE
        if (control$maxit == 0) {
            iter <- 0L
        } else {
            for (iter in seq.int(control$maxit)) {
                recomp <- FALSE # flag purely for trace output

                # Preconditioned CG-Steihaug trust-region step 
                # Minimise  -r_adj' p + 0.5 p' F p  s.t. ||p|| <= Delta
                # Uses stored fit$info_beta (p x p, O(p^2) matvec) with Jacobi
                # diagonal preconditioner M = diag(F) to reduce inner iterations.
                r_adj_sq  <- sum(adjusted_grad_beta^2)

                # Adaptive tolerance: loose early, tight near convergence
                #cg_tol_adapt <- min(0.5, sqrt(sqrt(r_adj_sq))) 

                # Eisenstat-Walker Choice 2 (Theorem 2.3)
                # https://softlib.rice.edu/pub/CRPC-TRs/reports/CRPC-TR94463.pdf
                gamma <- 0.9; alpha <- 2.0
                cg_tol_ew <- if (iter == 1) 0.5 else
                    min(0.5, gamma * (r_adj_sq / r_adj_sq_prev)^alpha)
                r_adj_sq_prev <- r_adj_sq   # store for next iteration

                #cg_tol = 0.1

                step <- cg_steihaug_pcg(
                    r_adj   = adjusted_grad_beta,
                    F_info  = fit$info_beta,
                    Delta   = Delta,
                    tol     = cg_tol_ew,
                    maxiter = cg_maxiter)

                # Predicted reduction (||r_adj||^2 objective) 
                # pred = 0.5(||r_adj||^2 - ||r_adj + (-F) p||^2), F stored as p x p.
                r_adj_model    <- adjusted_grad_beta - drop(fit$info_beta %*% step$p) # Matrix free?
                pred_reduction <- 0.5 * (r_adj_sq - sum(r_adj_model^2))

                betas_candidate <- betas + step$p
                theta_candidate <- c(betas_candidate, dispersion)

                fit_candidate <- try(
                    compute_fit(pars = theta_candidate, y = y, x = x, weights = weights,
                                offset = offset, family = family,
                                fixed_totals = fixed_totals, row_totals = row_totals,
                                no_dispersion = no_dispersion, nobs = nobs, nvars = nvars,
                                keep = keep, need_qr = TRUE, need_hatvalues = FALSE,
                                hatvalues_precomputed = hatvalues_cached,
                                need_inverse = needs_inverse),
                    silent = TRUE)

                if (inherits(fit_candidate, "try-error")) {
                    Delta <- shrink_factor * Delta
                    if (Delta < 1e-16) { warning("Trust region radius became too small"); break }
                    next
                }

                # Actual reduction (||r_adj||^2 objective) 
                # Compute adjusted score at the candidate point.
                adj_candidate   <- adjustment_function(theta_candidate, fit = fit_candidate,
                                                       level = 0, x, nobs, nvars, weights)
                r_adj_candidate <- fit_candidate$grad_beta + adj_candidate
                actual_reduction <- 0.5 * (r_adj_sq - sum(r_adj_candidate^2))

                # Reduction ratio (Nocedal & Wright 4.4) 
                # Avoids numerical issues when pred_reduction is small
                rho <- if (abs(pred_reduction) < 1e-16) {
                    if (actual_reduction > 0) 1.0 else 0.0
                } else {
                    actual_reduction / pred_reduction
                }

                # adapt radius (Nocedal & Wright Algo 4.1)
                if (rho < radius_shrink) {
                    Delta <- shrink_factor * Delta
                } else if (rho > radius_expand) {
                    Delta <- min(expand_factor * Delta, Delta_max)
                }

                # Accept/reject candidate, 0.1 threshold is fairly loose could be tightened to 0.25
                if (rho > eta) {
                    # Update with candidate 
                    betas <- betas_candidate
                    theta <- theta_candidate
                    fit <- fit_candidate
                    hat_accept_count <- hat_accept_count + 1L
 
                    # Periodic hat-value refresh: call qr.Q() on the already-stored
                    # qr_decomposition - no extra QR factorisation required.
                    if (hat_accept_count %% hat_update_freq == 0L) {
                        Q_tmp            <- qr.Q(fit$qr_decomposition)
                        hatvalues_cached <- .rowSums(Q_tmp * Q_tmp, nobs, nvars, TRUE)
                        Q_tmp            <- NULL   # free n×p matrix immediately
                        fit$hatvalues <- hatvalues_cached
                        recomp <- TRUE
                    }

                    step_components_beta <- compute_step_components(theta, level = 0, fit = fit)
                    adjusted_grad_beta <- with(step_components_beta, grad + adjustment)
                    accept_status <- "ACCEPT"
                } else {
                    accept_status <- "REJECT"
                }

                if (Delta < 1e-16) { warning("Trust region radius became too small"); break }

                step_components_zeta <- compute_step_components(theta, level = 1, fit = fit)
                adjusted_grad_zeta   <- if (no_dispersion) NA_real_ else
                                        with(step_components_zeta, grad + adjustment)

                # dispersion Newton step with step-halving
                if (!no_dispersion & df_residual > 0 &
                    !step_components_zeta$failed_inversion &
                    !step_components_zeta$failed_adjustment) {

                    step_zeta <- as.vector(adjusted_grad_zeta * step_components_zeta$inverse_info)
                    td_prev   <- transformed_dispersion
                    sf        <- 0L
                    
                    # Main step-halving loop: 
                    # keep halving the step until we get a positive dispersion 
                    # value or exceed max_step_factor
                    while (sf <= control$max_step_factor) {
                        td_new <- td_prev + 2^(-sf) * step_zeta
                        transformed_dispersion <- td_new
                        d_new  <- eval(control$inverseTrans)
                        if (is.finite(d_new) && d_new > 0) break
                        sf <- sf + 1L
                    }

                    # Recupte the fit at the new dispersion value but with the same betas
                    fit_new <- compute_fit(pars = c(betas, d_new), y = y, x = x,
                                           weights = weights, offset = offset,
                                           family = family,
                                           fixed_totals = fixed_totals,
                                           row_totals = row_totals,
                                           no_dispersion = no_dispersion,
                                           nobs = nobs, nvars = nvars, keep = keep,
                                           need_qr = FALSE, need_hatvalues = FALSE,
                                           need_inverse = FALSE)

                    # Rescale dispersion-dependent fields that aren't recomputed in fit_new
                    # Check these steps for completeness
                    scale                     <- d_new / dispersion
                    fit_new$info_beta         <- fit$info_beta         / scale
                    fit_new$inverse_info_beta <- fit$inverse_info_beta * scale
                    fit_new$hatvalues         <- fit$hatvalues
                    fit_new$qr_decomposition  <- fit$qr_decomposition
                    fit_new$R_matrix          <- fit$R_matrix
                    fit                       <- fit_new
                    dispersion                <- d_new
                    theta                     <- c(betas, dispersion)
                    step_components_beta      <- compute_step_components(theta, level = 0, fit = fit)
                    adjusted_grad_beta        <- with(step_components_beta, grad + adjustment)
                }

                step_zeta_conv <- if (no_dispersion || df_residual < 1 ||
                                      step_components_zeta$failed_inversion ||
                                      step_components_zeta$failed_adjustment) NA_real_ else
                                  as.vector(adjusted_grad_zeta * step_components_zeta$inverse_info)

                if (control$trace) {
                    cat(sprintf(
                        "Iter %3d: f=%.4e  ||r||=%.4e  ||p||=%.4e  dDisp=%.3e  Delta=%.3e  rho=%6.3f  %s%s%s\n",
                        iter,
                        0.5 * r_adj_sq,
                        sqrt(r_adj_sq),
                        sqrt(sum(step$p^2)),
                        if (is.na(step_zeta_conv)) 0 else abs(step_zeta_conv),
                        Delta, rho, accept_status,
                        if (step$on_boundary) " [BOUND]" else "",
                        if(recomp) " [HAT UPDATE]" else ""))
                }

                # Convergence: step-size criterion matching brglmFit_original.
                step_zeta_conv_abs <- if (is.na(step_zeta_conv)) NA_real_ else abs(step_zeta_conv)
                failed <- step_components_beta$failed_inversion ||
                          step_components_beta$failed_adjustment
                beta_converged <- max(abs(step$p)) < control$epsilon
                zeta_converged <- no_dispersion || df_residual < 1 ||
                                  is.na(step_zeta_conv_abs) ||
                                  step_zeta_conv_abs < control$epsilon
                if (failed || (beta_converged && zeta_converged)) break
            }
        }

        adjusted_grad_all[betas_names]              <- adjusted_grad_beta
        adjusted_grad_all["Transformed dispersion"] <- adjusted_grad_zeta
        betas_all[betas_names]                      <- betas

        ## Convergence analysis
        if ((failed | iter >= control$maxit) & !(is_correction)) {
            warning(paste("brglmFit: algorithm did not converge.",
                          "Try changing maxit, epsilon, slowit, or response_adjustment;",
                          "see ?brglm_control for defaults."), call. = FALSE)
            converged <- FALSE
        } else {
            converged <- TRUE
        }
        if (boundary) warning("brglmFit: algorithm stopped at boundary value", call. = FALSE)

        ## QR decomposition and fitted values are at the final value
        ## for the coefficients
        ## QR decomposition for cov.unscaled        
        if (!isTRUE(is_full_rank)) {
            x <- X_all
            betas_all[betas_names] <- betas
            betas <- betas_all
            betas[is.na(betas)] <- 0
            nvars <- nvars_all
        }

        ## If has_fixed_totals = TRUE, then scale fitted values before
        ## calculating QR decompositions, fitted values, etas,
        ## residuals and working_weights

        fit <- compute_fit(pars = c(betas, dispersion), y = y, x = x,
                           weights = weights, offset = offset, family = family,
                           fixed_totals = fixed_totals, row_totals = row_totals,
                           no_dispersion = no_dispersion, nobs = nobs, nvars = nvars,
                           keep = keep, need_qr = TRUE, need_hatvalues = FALSE)

        qr.Wx <- fit$qr_decomposition
        mus <- fit$mus
        etas <- fit$etas
        residuals <- with(fit, (y - mus) / d1mus)
        working_weights <- fit$working_weights

        ## info_transformed_dispersion will be NA if is_ML | is_AS_median | is_AS_mixed
        info_transformed_dispersion <- 1 / step_components_zeta$inverse_info
        if (is_ML | is_AS_median | is_AS_mixed) {
            transformed_dispersion <- eval(Trans1)
            d1zeta_1 <- eval(DD(Trans1, "dispersion", order = 1))
            adjusted_grad_all["Transformed dispersion"] <-
                adjusted_grad_all["Transformed dispersion"] / d1zeta_1
            info_transformed_dispersion <- info_transformed_dispersion / d1zeta_1^2
            control$transformation <- transformation1
            control$Trans          <- Trans1
            control$inverseTrans   <- inverseTrans1
        }

        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial" && (any(mus > 1 - eps) || any(mus < eps))) {
            warning("brglmFit: fitted probabilities numerically 0 or 1 occurred", call. = FALSE)
            boundary <- TRUE
        }
        if (family$family == "poisson" && any(mus < eps)) {
            warning("brglmFit: fitted rates numerically 0 occurred", call. = FALSE)
            boundary <- TRUE
        }
        if (df_residual == 0 & !no_dispersion) dispersion <- NA_real_
    }

    ## Working weights
    wt <- rep.int(0, nobs)
    wt[keep] <- working_weights[keep]
    names(wt) <- names(residuals) <- names(mus) <- names(etas) <-
        names(weights) <- names(y) <- ynames

    ## For the null deviance:
    ##
    ## If there is an intercept but not an offset then the ML fitted
    ## value is the weighted average and is calculated easily below if
    ## ML is used
    ##
    control0 <- control
    control0$maxit <- 1000
    if (customTransformation) control0$transformation <- transformation0

    if (intercept & missing_offset) {
        nullFit <- brglmFit(
            x = x[, "(Intercept)", drop = FALSE], y = y, weights = weights,
            offset = rep(0, nobs), family = family, intercept = TRUE,
            control = control0[c("epsilon", "maxit", "type", "transformation", "slowit")],
            start = if (no_dispersion) linkfun(mean(y)) else c(linkfun(mean(y)), 1))
        nullmus <- nullFit$fitted.values
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

    nulldev   <- sum(dev.resids(y, nullmus, weights))
    nulldf    <- nkeep - as.integer(intercept)
    deviance  <- sum(dev.resids(y, mus, weights))
    aic.model <- aic(y, n, mus, weights, deviance) + 2 * rank

    list(coefficients                = betas_all,
         residuals                   = residuals,
         fitted.values               = mus,
         R                           = if (!EMPTY) qr.R(qr.Wx),
         rank                        = rank,
         qr                          = if (!EMPTY) structure(
                                           qr.Wx[c("qr","rank","qraux","pivot","tol")],
                                           class = "qr"),
         family                      = family,
         linear.predictors           = etas,
         deviance                    = deviance,
         aic                         = aic.model,
         null.deviance               = nulldev,
         iter                        = iter,
         weights                     = wt,
         prior.weights               = weights,
         df.residual                 = df_residual,
         df.null                     = nulldf,
         y                           = y,
         converged                   = converged,
         boundary                    = boundary,
         dispersion                  = dispersion,
         dispersion_ML               = dispersion_ML,
         transformed_dispersion      = transformed_dispersion,
         info_transformed_dispersion = if (no_dispersion) NA_real_ else info_transformed_dispersion,
         grad                        = adjusted_grad_all,
         transformation              = control$transformation,
         type                        = control$type,
         control                     = control,
         class                       = "brglmFit")
}


#' Preconditioned CG-Steihaug trust-region subproblem solver
#'
#' Solves:  min_p { -r_adj' p + 0.5 p' F p }  s.t.  ||p|| <= Delta
#' where F = fit$info_beta (p x p Fisher information, already formed).
#'
#' Uses a Jacobi (diagonal) preconditioner  M = diag(F).  This costs
#' one O(p) divide per inner step and replaces the CG convergence rate
#' governed by kappa(F) with kappa(M^{-1} F), which is substantially
#' smaller when predictors have very different scales or block-correlated
#' structure.  In practice this halves or quarters the number of inner iterations k 
#' relative to unpreconditioned CG, directly reducing the k * O(p^2) inner-loop cost.
#'
#'
#' @param r_adj   Adjusted score vector (length p)
#' @param F_info  Fisher information matrix (p x p), i.e. fit$info_beta
#' @param Delta   Trust-region radius
#' @param tol     Relative residual tolerance (use Eisenstat-Walker adaptive
#'                value: min(0.5, sqrt(||r_adj||)) for best performance)
#' @param maxiter Maximum CG iterations
#'
#' @references Steihaug (1983) SIAM J. Numer. Anal. 20(3), 626-637.
#'             Nocedal & Wright (2006) Numerical Optimization, Alg. 7.2.
#'             Eisenstat & Walker (1996) SIAM J. Optim. 6(4), 1190-1206.
cg_steihaug_pcg <- function(r_adj, F_info, Delta, tol = 0.1, maxiter = 50) {
    # Jacobi preconditioner: M = diag(F), M^{-1} v = v / diag(F)
    # Guard against near-zero diagonal entries
    d_F <- diag(F_info)
    d_F <- pmax(d_F, .Machine$double.eps * max(d_F))
    Minv <- function(v) v / d_F              # O(p), just element-wise division

    z <- numeric(length(r_adj))
    r <- r_adj
    y <- Minv(r)                       # preconditioned residual
    d <- y
    ry <- sum(r * y)
    ry0 <- ry                            # for convergence check

    for (j in seq_len(maxiter)) {
        Fd  <- drop(F_info %*% d)            # O(p²), maybe matrix free matvec in future?
        dFd <- sum(d * Fd)

        if (dFd <= 0) {
            tau <- find_boundary_step(z, d, Delta)
            return(list(p = z + tau * d, on_boundary = TRUE, cg_iter = j))
        }

        alpha <- ry / dFd
        z_new <- z + alpha * d

        if (sqrt(sum(z_new^2)) >= Delta) {
            tau <- find_boundary_step(z, d, Delta)
            return(list(p = z + tau * d, on_boundary = TRUE, cg_iter = j))
        }

        z <- z_new
        r <- r - alpha * Fd
        y <- Minv(r)
        ry_new <- sum(r * y)

        # convergence in the M-norm of the residual, relative to initial
        if (sqrt(abs(ry_new)) < tol * sqrt(abs(ry0)))
            return(list(p = z, on_boundary = FALSE, cg_iter = j))

        d  <- y + (ry_new / ry) * d
        ry <- ry_new
    }
    list(p = z, on_boundary = FALSE, cg_iter = maxiter)
}

#' Boundary step: finds tau >= 0 s.t. ||z + tau*d|| = Delta
find_boundary_step <- function(z, d, Delta) {
    a    <- sum(d^2)
    b    <- 2 * sum(z * d)
    cc   <- sum(z^2) - Delta^2
    disc <- b^2 - 4 * a * cc
    if (disc < 0) return(0)
    t1  <- (-b + sqrt(disc)) / (2 * a)
    t2  <- (-b - sqrt(disc)) / (2 * a)
    pos <- c(t1, t2)[c(t1, t2) > 0]
    if (length(pos) == 0) 0 else min(pos)
}


#' Extract model coefficients from [`"brglmFit"`][brglmFit] objects
#'
#' @inheritParams stats::coef
#' @param model one of `"mean"` (default), `"dispersion"`, `"full"`,
#'     to return the estimates of the parameters in the linear
#'     prediction only, the estimate of the dispersion parameter only,
#'     or both, respectively.
#'
#' @details
#'
#' See [coef()] for more details.
#'
#' @seealso [coef()]
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
        transDisp  # always on the scale of the transformed dispersion
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

## Almost all code here is from stats:::print.summary.glm with minor modifications
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