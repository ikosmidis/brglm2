## Creates various objects that are necessary for the operations of
## brglmFit and puts them into an environment
set_brglmFit_env <- function(y, x, weights, offset,
                             family, fixed_totals,
                             control,
                             singular.ok,
                             start,
                             env) {
   
    ## Useful indicators
    is_ML <- control$type == "ML"
    is_AS_median <- control$type == "AS_median"
    is_AS_mixed <- control$type == "AS_mixed"
    is_correction <- control$type == "correction"
    has_sparse_x <- is(x, "sparseMatrix")
    has_no_dispersion <- family$family %in% c("poisson", "binomial")
    boundary <- converged <- FALSE

    ## Dispersion transformations
    if (is_ML | is_AS_median | is_AS_mixed) {
        transformation1 <- control$transformation
        Trans1 <- control$Trans
        inverseTrans1 <- control$inverseTrans
        control$transformation <- "identity"
        control$Trans <- expression(dispersion)
        control$inverseTrans <- expression(transformed_dispersion)
    }
    d1_transformed_dispersion <- DD(control$Trans, "dispersion", order = 1)
    d2_transformed_dispersion <- DD(control$Trans, "dispersion", order = 2)
    
    has_fixed_totals <- FALSE
    row_totals <- NULL 
    if (isTRUE(family$family == "poisson" && !is.null(fixed_totals))) {
        row_totals <-  as.vector(tapply(y, fixed_totals, sum))[fixed_totals]
        has_fixed_totals <- TRUE
    }

    ## Model matrix, response, names etc
    if (!(has_sparse_x | is.matrix(x))) {
        ## Ensure that x is a matrix        
        x <- as.matrix(x)
    }
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

    ## family and link functions
    ok_links <- c("logit", "probit", "cauchit",
                  "cloglog", "identity", "log",
                  "sqrt", "inverse")


    if (isTRUE(family$family %in% c("quasi", "quasibinomial", "quasipoisson"))) {
        stop("`brglmFit` does not currently support the `quasi`, `quasipoisson` and `quasibinomial` families.")
    }

    family <- enrichwith::enrich(family, with = c("d1afun", "d2afun", "d3afun", "d1variance"))
    if ((family$link %in% ok_links) | (grepl("mu\\^", family$link))) {
        ## Enrich the link object with d2mu.deta and update family object
        linkglm <- make.link(family$link)
        linkglm <- enrichwith::enrich(linkglm, with = "d2mu.deta")
        ## Put everything into the family object
        family[names(linkglm)] <- linkglm
    }

    variance <- family$variance
    d1variance <- family$d1variance
    linkinv <- family$linkinv
    linkfun <- family$linkfun
    if (!is.function(variance) || !is.function(linkinv)) {
        stop("'family' argument seems not to be a valid family object", call. = FALSE)
    }
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    ## If the family is custom then d2mu.deta cannot survive when
    ## passing throguh current family functions. But mu.eta does; so
    ## we compute d2mu.deta numerically; this allows also generality,
    ## as the users can then keep their custom link implementations
    ## unaltered. Issue is scalability, due to the need of evaluating
    ## n numerical derivatives
    if (is.null(family$d2mu.deta)) {
        d2mu.deta <- function(eta) {
            numDeriv::grad(mu.eta, eta)
        }
    }
    else {
        d2mu.deta <- family$d2mu.deta
    }
    d1afun <- family$d1afun
    d2afun <- family$d2afun
    d3afun <- family$d3afun
    simulate <- family$simulate

    valid_eta <- unless_null(family$valideta, function(eta) TRUE)
    valid_mu <- unless_null(family$validmu, function(mu) TRUE)

    mustart <- NULL
    etastart <- NULL

    ## Initialize as prescribed in family
    eval(family$initialize)
    
    ## Prepare the brglmFit environment
    Env <- named_list(control,
                      is_ML,
                      is_AS_median,
                      is_AS_mixed,
                      is_correction,
                      has_no_dispersion,
                      has_sparse_x,
                      row_totals,
                      has_fixed_totals,
                      singular.ok,
                      d1_transformed_dispersion,
                      d2_transformed_dispersion,
                      ## model matrix, response, offset, etc
                      y,
                      x,
                      betas_names,
                      ynames,
                      converged,
                      nobs,
                      nvars,
                      EMPTY,
                      weights,
                      missing_offset,
                      offset,
                      ## family and link functions
                      variance,
                      d1variance,
                      linkinv,
                      linkfun,
                      dev.resids,
                      aic,
                      mu.eta,
                      d2mu.deta,
                      d1afun,
                      d2afun,
                      d3afun,
                      simulate,
                      valid_mu,
                      valid_eta,
                      start)
    if (env) list2env(Env) else Env
}



## ## ## Benchmarks
## ## ### g2 environment creation is faster than g1, but g1 is cleaner when
## ## ### objects depends on objects that have been generated before
## g1 <- function(d) {
##     a1 <- 1
##     a2 <- 2
##     a3 <- 3
##     a4 <- 4
##     a5 <- 5
##     a6 <- 6
##     a7 <- 7
##     a8 <- 8
##     a9 <- 9
##     a10 <- 10
##     e <- named_list(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
##     list2env(e)
## }


## g2 <- function(d) {
##     e <- new.env()
##     e$a1 <- 1
##     e$a2 <- 2
##     e$a3 <- 3
##     e$a4 <- 4
##     e$a5 <- 5
##     e$a6 <- 6
##     e$a7 <- 7
##     e$a8 <- 8
##     e$a9 <- 9
##     e$a10 <- 10
##     e
## }

## g3 <- function(d) {
##     a1 <- 1
##     a2 <- 2
##     a3 <- 3
##     a4 <- 4
##     a5 <- 5
##     a6 <- 6
##     a7 <- 7
##     a8 <- 8
##     a9 <- 9
##     a10 <- 10
##     named_list(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
## }
