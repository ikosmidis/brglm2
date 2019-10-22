## Typically applied to the enviornment from create_brglmFit_env
prepare_empty_model <- function(Env) {
    nobs <- Env$nobs
    etas <- rep.int(0, nobs) + Env$offset
    if (!Env$valid_eta(etas)) {
        stop("invalid linear predictor values in empty model", call. = FALSE)
    }
    mus <- Env$linkinv(etas)
    if (!Env$valid_mu(mus)) {
        stop("invalid fitted means in empty model", call. = FALSE)
    }
    metas <- Env$mu.eta(etas)
    Env$working_weights <- ((Env$weights * metas^2) / Env$variance(mus))^0.5
    Env$residuals <- (Env$y - mus)/metas
    Env$keep <- rep(TRUE, nobs)
    Env$boundary <- Env$converged <- TRUE
    Env$betas_all <- numeric()
    Env$rank <- 0
    Env$iter <- 0L   
    Env$etas <- etas
    Env$mus <- mus
    Env
}
