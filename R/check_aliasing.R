## Typically applied to the enviornment from create_brglmFit_env
check_aliasing <- function(Env) {    
    qrx <- qr(Env$x)
    rank <- qrx$rank
    nvars <- Env$nvars
    betas_names <- Env$betas_names
    
    is_full_rank <- rank == nvars
    
    if (!Env$singular.ok && !is_full_rank) {
        stop("singular fit encountered")
    }
    if (!isTRUE(is_full_rank)) {
        aliased <- qrx$pivot[seq.int(qrx$rank + 1, nvars)]
        nvars_all <- nvars
        betas_names_all <- betas_names
        nvars <- nvars - sum(aliased)
        betas_names <- betas_names[-aliased]
    }
    else {
        aliased <- logical(nvars)
        nvars_all <- nvars
        betas_names_all <- betas_names
    }
    betas_all <- structure(rep(NA_real_, nvars_all), .Names = betas_names_all)
    keep <- Env$weights > 0

    ## Return enviornment
    Env$keep <- keep
    nkeep <- sum(keep)
    Env$df_residual <- nkeep - rank
    Env$nkeep <- nkeep
    Env$nvars <- nvars
    Env$nvars_all <- nvars_all
    Env$betas_names <- betas_names
    Env$betas_names_all <- betas_names_all
    Env$betas_all <- betas_all
    Env$aliased <- aliased
    Env$rank <- rank
    Env$is_full_rank <- is_full_rank
    Env
}
