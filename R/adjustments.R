AS_mean_adjustment <- function(pars, level = 0, fit = NULL) {
    if (is.null(fit)) {
        fit <- key_quantities(pars, y = y, level = level, qr = TRUE)
    }
    with(fit, {
        if (level == 0) {
            hatvalues <- hat_values(pars, fit = fit)
            ## Use only observations with keep = TRUE to ensure that no division with zero takes place
            return(.colSums(0.5 * hatvalues * d2mus/d1mus * x, nobs, nvars, TRUE))
        }
        if (level == 1) {
            s1 <- sum(weights^3 * d3afuns, na.rm = TRUE)
            s2 <- sum(weights^2 * d2afuns, na.rm = TRUE)
            return((nvars - 2)/(2 * dispersion) + s1/(2 * dispersion^2 * s2))
        }
    })
}

AS_Jeffreys_adjustment <- function(pars, level = 0, fit = NULL) {
    if (is.null(fit)) {
        fit <- key_quantities(pars, y = y, level = level, qr = TRUE)
    }
    with(fit, {
        if (level == 0) {
            hatvalues <- hat_values(pars, fit = fit)
            ## Use only observations with keep = TRUE to ensure that no division with zero takes place
            return(2 * control$a * .colSums(0.5 * hatvalues * (2 * d2mus/d1mus - d1varmus * d1mus / varmus) * x, nobs, nvars, TRUE))
        }
        if (level == 1) {
            s1 <- sum(weights^3 * d3afuns, na.rm = TRUE)
            s2 <- sum(weights^2 * d2afuns, na.rm = TRUE)
            return(2 * control$a * (-(nvars + 4)/(2 * dispersion) + s1/(2 * dispersion^2 * s2)))
        }
    })
}


## Implementation by Euloge Clovis Kenne Pagui, 20 April 2017 (kept here for testing)
## AS_median_adjustment <- function(pars, level = 0, fit = NULL) {
##     if (is.null(fit)) {
##         fit <- key_quantities(pars, y = y, level = level, qr = TRUE)
##     }
##     with(fit, {
##         if (level == 0) {
##             R_matrix <- qr.R(qr_decomposition)
##             info <- precision * crossprod(R_matrix)
##             inverse_info <- try(dispersion * tcrossprod(solve(R_matrix)))
##             nu_r_s_t <- nu_r_st <- array(0,c(nvars,nvars,nvars))
##             k1 <- k2 <- k3 <- b.beta <- rep(NA,nvars)
##             for (r in 1:nvars)
##             {
##               nu_r_s_t[r,,] <- t(x)%*%((working_weights*d1mus*d1varmus*x[,r]/varmus)*x)
##               nu_r_st[r,,] <- -t(x)%*%((working_weights*d1mus*(d1varmus/varmus-d2mus/d1mus^2)*x[,r])*x)
##             }
##             k2 <- 1/diag(inverse_info)
##             for (r in 1:nvars)
##             {
##               sum_s1 <- rep(0,nvars)
##               sum_s3 <- rep(0,nvars)
##               nu.tu <- inverse_info-outer(inverse_info[,r]*k2[r],inverse_info[,r])
##               for (s in 1:nvars){
##                 sum_s1[s] <- sum(diag(nu.tu%*%(nu_r_s_t[s,,]+nu_r_st[s,,])))
##                 sum_s3[s] <- sum((inverse_info[r,]%*%nu_r_s_t[s,,])*inverse_info[r,])
##               }
##               barb1r <- sum(sum_s1*inverse_info[r,])
##               barb2r <- k2[r]*sum(sum_s3*inverse_info[r,])
##               b.beta[r]<- barb1r/2+barb2r/6
##             }
##             return(info%*%b.beta)
##         }
##         if (level == 1) {
##               s1 <- sum(weights^3 * d3afuns, na.rm = TRUE)
##               s2 <- sum(weights^2 * d2afuns, na.rm = TRUE)
##               return(nvars/(2*dispersion) + s1/(6*dispersion^2*s2))
##         }
##     })
## }

## Implementation by Ioannis Kosmidis, 02 May 2017
AS_median_adjustment <- function(pars, level = 0, fit = NULL) {
    if (is.null(fit)) {
        fit <- key_quantities(pars, y = y, level = level, qr = TRUE)
    }
    with(fit, {
        if (level == 0) {
            hatvalues <- hat_values(pars, fit = fit)
            R_matrix <- qr.R(qr_decomposition)
            info_unscaled <- crossprod(R_matrix)
            inverse_info_unscaled <- chol2inv(R_matrix)
            ## FIXME: There is 1) definitely a better way to do this, 2) no time...
            b_vector <- numeric(nvars)
            for (j in seq.int(nvars)) {
                inverse_info_unscaled_j <- inverse_info_unscaled[j, ]
                vcov_j <- tcrossprod(inverse_info_unscaled_j) / inverse_info_unscaled_j[j]
                hats_j <- .rowSums((x %*% vcov_j) * x, nobs, nvars, TRUE) * working_weights
                b_vector[j] <- inverse_info_unscaled_j %*% .colSums(x * (hats_j * (d1mus * d1varmus / (6 * varmus) - 0.5 * d2mus/d1mus)), nobs, nvars, TRUE)
            }
            return(.colSums(0.5 * hatvalues * d2mus / d1mus * x, nobs, nvars, TRUE) +
                   info_unscaled %*% b_vector)
        }
        if (level == 1) {
            s1 <- sum(weights^3 * d3afuns, na.rm = TRUE)
            s2 <- sum(weights^2 * d2afuns, na.rm = TRUE)
            return(nvars/(2 * dispersion) + s1/(6 * dispersion^2 * s2))
        }
    })
}

AS_mixed_adjustment <- function(pars, level = 0, fit = NULL) {
    if (is.null(fit)) {
        fit <- key_quantities(pars, y = y, level = level, qr = TRUE)
    }
    with(fit, {
        if (level == 0) {
            hatvalues <- hat_values(pars, fit = fit)
            ## Use only observations with keep = TRUE to ensure that no division with zero takes place
            return(.colSums(0.5 * hatvalues * d2mus/d1mus * x, nobs, nvars, TRUE))
        }
        if (level == 1) {
            s1 <- sum(weights^3 * d3afuns, na.rm = TRUE)
            s2 <- sum(weights^2 * d2afuns, na.rm = TRUE)
            return(nvars/(2 * dispersion) + s1/(6 * dispersion^2 * s2))
        }
    })
}
