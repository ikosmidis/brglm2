#############################################################
### Median bias reduction in poisson regression  ############
#############################################################
mod1 <- function(X,mu,InfoInv) {
    X<-as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    nu_r_s_t  <- array(0,c(p,p,p))
    k1 <- k2 <- k3 <- rep(NA,p)
    for(r in 1:p) {
        nu_r_s_t[r,,] = t(X)%*%(X[,r]*mu*X)
    }
    k2 <- 1/diag(InfoInv)
    barA <- InfoInv*k2
    for (r in 1:p) {
        sum_s1 <- rep(0,p)
        sum_s3 <- rep(0,p)
        nu.tu <- InfoInv-outer(InfoInv[,r]*k2[r],InfoInv[,r])
        for (s in 1:p){
            sum_s1[s] <- sum(diag(nu.tu%*%(nu_r_s_t[s,,])))
            sum_s3[s] <- sum((barA[r,]%*%nu_r_s_t[s,,])*barA[r,])
        }
        k1[r] <- -0.5*sum(sum_s1*barA[r,])
        k3[r] <- sum(sum_s3*barA[r,])
    }
    return(-k1/k2+k3/(6*k2^2))
}


mbrpr <- function(par,y,X,eps=1e-06,maxit=500) {
    step <- .Machine$integer.max
    nIter <- 0
    test <- TRUE
    while ( test & (nIter < maxit)) {
        nIter <- nIter + 1
        eta = drop(X%*%par)
        mu <- exp(eta)
        info <- t(X)%*%(W<-diag(mu))%*%X
        score <- t(X)%*%(y-mu)
        InfoInv <- try(chol2inv(chol(info)),TRUE)
        if(failedInv <- inherits(InfoInv, "try-error"))  {
            warning("failed to invert the information matrix: iteration stopped prematurely")
            break
        }
        mod <- mod1(X, mu, InfoInv)
        modscore <- score + info%*%mod
        y.adj <- X%*%(par+mod)+ solve(diag(mu))%*%(y-mu)
        par <- InfoInv%*%t(X)%*%(W%*%y.adj)
        test <- sqrt(crossprod(drop(modscore))) > eps
    }
    converged <- nIter < maxit
    list(coefficients=drop(par),converged=converged,nIter=nIter,InfoInv=InfoInv)
}



