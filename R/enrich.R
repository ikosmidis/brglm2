
## quasi, quasibinomial and quasipoisson are not currently supported
enrichLink <- function(linkObject) {
    mu.eta <- linkObject$mu.eta
    linkinv <- linkObject$linkinv
    link <- linkObject$name

    if (grepl("mu\\^", link) & (link != "1/mu^2")) {
        linkObject$d2mu.deta <- function(eta) {
            (1/lambda) * (1/lambda - 1) * eta^(1/lambda - 2)
        }
        linkObject$d3mu.deta <- function(eta) {
            (1/lambda) * (1/lambda - 1) * (1/lambda - 2) * eta^(1/lambda - 3)
        }
        ## set the environment so lambda can be found
        environment(linkObject$d2mu.deta) <- environment(linkObject$linkfun)
        environment(linkObject$d3mu.deta) <- environment(linkObject$linkfun)
        return(linkObject)
    }
    linkObject$d2mu.deta <- switch(link,
                              "logit" = function(eta) {
                                  mu.eta(eta) * (1 - 2 * linkinv(eta))
                              },
                              "probit" = function(eta) {
                                  -eta * pmax(dnorm(eta),.Machine$double.eps)
                              },
                              "cauchit" = function(eta) {
                                  -2 * pi * eta * pmax(dcauchy(eta)^2, .Machine$double.eps)
                              },
                              "cloglog" = function(eta) {
                                  (1 - exp(eta)) * pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
                              },
                              "identity" = function(eta) {
                                  rep.int(0, length(eta))
                              },
                              "log" = function(eta) {
                                  pmax(exp(eta), .Machine$double.eps)
                              },
                              "sqrt" = function(eta) {
                                  rep.int(2, length(eta))
                              },
                              "1/mu^2" = function(eta) {
                                  3/(4 * eta^2.5)
                              },
                              "inverse" = function(eta) {
                               2/(eta^3)
                              },
                              ## else :
                              stop(sQuote(link), " link not recognised")
                              )# end switch(.)


    linkObject$d3mu.deta <- switch(link,
                                   "logit" = function(eta) {
                                       mu.eta(eta) * (1 - 6 * mu.eta(eta))
                                   },
                                   "probit" = function(eta) {
                                       (eta^2 - 1) * pmax(dnorm(eta),.Machine$double.eps)
                                   },
                                   "cauchit" = function(eta) {
                                       -2 * pi * pmax(dcauchy(eta)^2, .Machine$double.eps) + 8 * eta^2 * pi^2 * pmax(dcauchy(eta)^3, .Machine$double.eps)
                                   },
                                   "cloglog" = function(eta) {
                                       ((1 - exp(eta))^2 - exp(eta))  * pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
                                   },
                                   "identity" = function(eta) {
                                       rep.int(0, length(eta))
                                   },
                                   "log" = function(eta) {
                                       pmax(exp(eta), .Machine$double.eps)
                                   },
                                   "sqrt" = function(eta) {
                                       rep.int(0, length(eta))
                                   },
                                   "1/mu^2" = function(eta) {
                                       -15/(8 * eta^3.5)
                                   },
                                   "inverse" = function(eta) {
                                       -6/(eta^4)
                                   },
                                   ## else :
                                   stop(sQuote(link), " link not recognised")
                                   )# end switch(.)

    linkObject
}

enrichFamily <- function(familyObject) {
    family <- familyObject$family
    link <- familyObject$link
    stats <- enrichLink(make.link(link))
    familyObject$d2mu.deta <- stats$d2mu.deta
    familyObject$d3mu.deta <- stats$d3mu.deta
    ## Derivatives of the variance function
    familyObject$d1variance <- switch(family,
                                     "poisson" = function(mu) {
                                         1
                                     },
                                     "gaussian" = function(mu) {
                                         0
                                     },
                                     "binomial" = function(mu) {
                                         1 - 2 * mu
                                     },
                                     "Gamma" = function(mu) {
                                         2 * mu
                                     },
                                     "inverse.gaussian" = function(mu) {
                                         3 * mu^2
                                     },
                                     stop(sQuote(family), " family not supported"))
    familyObject$d2variance <- switch(family,
                                      "poisson" = function(mu) {
                                          0
                                      },
                                      "gaussian" = function(mu) {
                                          0
                                      },
                                      "binomial" = function(mu) {
                                          -2
                                      },
                                      "Gamma" = function(mu) {
                                          2
                                      },
                                      "inverse.gaussian" = function(mu) {
                                          6 * mu
                                      },
                                      stop(sQuote(family), " family not supported"))
    ## Derivatives of the afunction
    familyObject$d1afun <- switch(family,
                                  "gaussian" = function(zeta) {
                                      -1/zeta
                                  },
                                  "Gamma" = function(zeta) {
                                      -2*psigamma(-zeta, 0) + 2*log(-zeta)
                                  },
                                  "inverse.gaussian" = function(zeta) {
                                      -1/zeta
                                  })

    familyObject$d2afun <- switch(family,
                                  "gaussian" = function(zeta) {
                                      1/zeta^2
                                  },
                                  "Gamma" = function(zeta) {
                                      2*psigamma(-zeta, 1) + 2/zeta
                                  },
                                  "inverse.gaussian" = function(zeta) {
                                      1/zeta^2
                                  })

    familyObject$d3afun <- switch(family,
                                  "gaussian" = function(zeta) {
                                      -2/zeta^3
                                  },
                                  "Gamma" = function(zeta) {
                                      -2*psigamma(-zeta, 2) - 2/zeta^2
                                  },
                                  "inverse.gaussian" = function(zeta) {
                                      -2/zeta^3
                                  })
    familyObject
}




## power <- function(lambda = 1)
## {
##     if(!is.numeric(lambda) || is.na(lambda))
##         stop("invalid argument 'lambda'")
##     if(lambda <= 0) return(make.link("log"))
##     if(lambda == 1) return(make.link("identity"))
##     linkfun <- function(mu) mu^lambda
##     linkinv <- function(eta)
##         pmax(eta^(1/lambda), .Machine$double.eps)
##     mu.eta <- function(eta)
##         pmax((1/lambda) * eta^(1/lambda - 1), .Machine$double.eps)
##     dmu.deta <- function(eta)
##         (1/lambda) * (1/lambda - 1) * eta^(1/lambda - 2)
##     ddmu.ddeta <- function(eta)
##         (1/lambda) * (1/lambda - 1) * (1/lambda - 2) * eta^(1/lambda - 3)
##     valideta <- function(eta) all(eta>0)
##     link <- paste("mu^", round(lambda, 3), sep="")
##     structure(list(linkfun = linkfun, linkinv = linkinv,
##                    mu.eta = mu.eta,
##                    dmu.deta = dmu.deta,
##                    ddmu.ddeta = ddmu.ddeta,
##                    valideta = valideta, name = link),
##               class="link-glm")
## }

## ## Written by Simon Davies Dec 1995
## ## Modified by Thomas Lumley 26 Apr 97
## ## added valideta(eta) function..
## ## Modified by Ioannis Kosmidis
## ## added dmu.deta(eta) functions
## make.link <- function (link)
## {
##     switch(link,
##            "logit" = {
##                linkfun <- function(mu)
##                    .Call(stats:::C_logit_link, mu, PACKAGE="stats")
##                linkinv <- function(eta)
##                    .Call(stats:::C_logit_linkinv, eta, PACKAGE="stats")
##                mu.eta <- function(eta)
##                    .Call(stats:::C_logit_mu_eta, eta, PACKAGE="stats")
##                dmu.deta <- function(eta)
##                  .Call(stats:::C_logit_mu_eta, eta, PACKAGE="stats")*
##                    (1 - 2 * .Call(stats:::C_logit_linkinv, eta, PACKAGE="stats"))
##                ## temporarily
##                ddmu.ddeta <- function(eta)
##                    exp(eta)*(1 - 4*exp(eta) + exp(2*eta))/(1 + exp(eta))^4
##                valideta <- function(eta) TRUE
##            },
##            "probit" = {
##                linkfun <- function(mu) qnorm(mu)
##                linkinv <- function(eta) {
##                    thresh <- - qnorm(.Machine$double.eps)
##                    eta <- pmin(pmax(eta, -thresh), thresh)
##                    pnorm(eta)
##                }
##                mu.eta <- function(eta)
##                    pmax(dnorm(eta),.Machine$double.eps)
##                dmu.deta <- function(eta) {
##                  -eta * pmax(dnorm(eta),.Machine$double.eps)
##                }
##                ddmu.ddeta <- function(eta) {
##                  (eta^2 - 1) * pmax(dnorm(eta),.Machine$double.eps)
##                }
##                valideta <- function(eta) TRUE
##            },
##            "cauchit" = {
##                linkfun <- function(mu) qcauchy(mu)
##                linkinv <- function(eta) {
##                    thresh <- -qcauchy(.Machine$double.eps)
##                    eta <- pmin(pmax(eta, -thresh), thresh)
##                    pcauchy(eta)
##                }
##                mu.eta <- function(eta)
##                    pmax(dcauchy(eta), .Machine$double.eps)
##                dmu.deta <- function(eta)
##                    -2 * pi * eta * pmax(dcauchy(eta)^2, .Machine$double.eps)
##                ddmu.ddeta <- function(eta)
##                    -2 * pi * pmax(dcauchy(eta)^2, .Machine$double.eps) + 8 * eta^2 * pi^2 * pmax(dcauchy(eta)^3, .Machine$double.eps)
##                valideta <- function(eta) TRUE
##            },
##            "cloglog" = {
##                linkfun <- function(mu) log(-log(1 - mu))
##                linkinv <- function(eta)
##                    pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps),
##                         .Machine$double.eps)
##                mu.eta <- function(eta) {
##                    eta <- pmin(eta, 700)
##                    pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
##                }
##                dmu.deta <- function(eta) {
##                  (1 - exp(eta)) * pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
##                }
##                ddmu.ddeta <- function(eta) {
##                  ((1 - exp(eta))^2 - exp(eta))  * pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
##                }
##                valideta <- function(eta) TRUE
##            },
##            "identity" = {
##                linkfun <- function(mu) mu
##                linkinv <- function(eta) eta
##                mu.eta <- function(eta) rep.int(1, length(eta))
##                dmu.deta <- function(eta) rep.int(0, length(eta))
##                ddmu.ddeta <- function(eta) rep.int(0, length(eta))
##                valideta <- function(eta) TRUE
##            },
##            "log" = {
##                linkfun <- function(mu) log(mu)
##                linkinv <- function(eta)
##                    pmax(exp(eta), .Machine$double.eps)
##                mu.eta <- function(eta)
##                    pmax(exp(eta), .Machine$double.eps)
##                dmu.deta <- function(eta)
##                    pmax(exp(eta), .Machine$double.eps)
##                ddmu.ddeta <- function(eta)
##                    pmax(exp(eta), .Machine$double.eps)
##                valideta <- function(eta) TRUE
##            },
##            "sqrt" = {
##                linkfun <- function(mu) sqrt(mu)
##                linkinv <- function(eta) eta^2
##                mu.eta <- function(eta) 2 * eta
##                dmu.deta <- function(eta) rep.int(2, length(eta))
##                ddmu.ddeta <- function(eta) rep.int(0, length(eta))
##                valideta <- function(eta) all(eta>0)
##            },
##            "1/mu^2" = {
##                linkfun <- function(mu) 1/mu^2
##                linkinv <- function(eta) 1/sqrt(eta)
##                mu.eta <- function(eta) -1/(2 * eta^1.5)
##                dmu.deta <- function(eta) 3/(4 * eta^2.5)
##                ddmu.ddeta <- function(eta) -15/(8 * eta^3.5)
##                valideta <- function(eta) all(eta>0)
##            },
##            "inverse" = {
##                linkfun <- function(mu) 1/mu
##                linkinv <- function(eta) 1/eta
##                mu.eta <- function(eta) -1/(eta^2)
##                dmu.deta <- function(eta) 2/(eta^3)
##                ddmu.ddeta <- function(eta) -6/(eta^4)
##                valideta <- function(eta) all(eta!=0)
##            },
##            ## else :
##            stop(sQuote(link), " link not recognised")
##            )# end switch(.)
##     structure(list(linkfun = linkfun, linkinv = linkinv,
##                    mu.eta = mu.eta, dmu.deta = dmu.deta,
##                    ddmu.ddeta = ddmu.ddeta,
##                    valideta = valideta, name = link),
##               class="link-glm")
## }


## quasipoisson <- function (link = "log")
## {
##     linktemp <- substitute(link)
##     if (!is.character(linktemp)) {
## 	linktemp <- deparse(linktemp)
##         ## the idea here seems to be that we can have a character variable
##         ## 'link' naming the link.
## 	if (linktemp == "link") {
##             warning("use of quasipoisson(link=link) is deprecated\n",
##                     domain = NA)
## 	    linktemp <- eval(link)
##             if(!is.character(linktemp) || length(linktemp) != 1L)
##                 stop("'link' is invalid", domain=NA)
##         }
##     }

##     okLinks <- c("log", "identity", "sqrt")
##     if (linktemp %in% okLinks)
##         stats <- make.link(linktemp)
##     else if (is.character(link)) {
##         stats <- make.link(link)
##         linktemp <- link
##     } else {
##         ## what else shall we allow?  At least objects of class link-glm.
##         if(inherits(link, "link-glm")) {
##             stats <- link
##             if(!is.null(stats$name)) linktemp <- stats$name
##         } else {
## 	    stop(gettextf('link "%s" not available for quasipoisson family; available links are %s',
## 			  linktemp, paste(sQuote(okLinks), collapse =", ")),
## 		 domain = NA)
##         }
##     }
##     variance <- function(mu) mu
##     dvariance <- function(mu) rep.int(1, length(mu))
##     ddvariance <- function(mu) rep.int(0, length(mu))
##     validmu <- function(mu) all(mu>0)
##     dev.resids <- function(y, mu, wt)
## 	2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
##     aic <- function(y, n, mu, wt, dev) NA
##     initialize <- expression({
## 	if (any(y < 0))
## 	    stop("negative values not allowed for the quasiPoisson family")
## 	n <- rep.int(1, nobs)
## 	mustart <- y + 0.1
##     })
##     structure(list(family = "quasipoisson",
## 		   link = linktemp,
## 		   linkfun = stats$linkfun,
## 		   linkinv = stats$linkinv,
## 		   variance = variance,
##                    dvariance = dvariance,
##                    ddvariance = ddvariance,
## 		   dev.resids = dev.resids,
## 		   aic = aic,
## 		   mu.eta = stats$mu.eta,
##                    dmu.deta = stats$dmu.deta,
## 		   initialize = initialize,
## 		   validmu = validmu,
## 		   valideta = stats$valideta),
## 	      class = "family")
## }

## quasibinomial <- function (link = "logit")
## {
##     linktemp <- substitute(link)
##     if (!is.character(linktemp)) {
## 	linktemp <- deparse(linktemp)
##         ## the idea here seems to be that we can have a character variable
##         ## 'link' naming the link.
## 	if (linktemp == "link") {
##             warning("use of quasibinomial(link=link) is deprecated\n", domain = NA)
## 	    linktemp <- eval(link)
##             if(!is.character(linktemp) || length(linktemp) != 1L)
##                 stop("'link' is invalid", domain=NA)
##         }
##     }

##     okLinks <- c("logit", "probit", "cloglog", "cauchit", "log")
##     if (linktemp %in% okLinks)
##         stats <- make.link(linktemp)
##     else if (is.character(link)) {
##         stats <- make.link(link)
##         linktemp <- link
##     } else {
##         ## what else shall we allow?  At least objects of class link-glm.
##         if(inherits(link, "link-glm")) {
##             stats <- link
##             if(!is.null(stats$name)) linktemp <- stats$name
##         } else {
## 	    stop(gettextf('link "%s" not available for quasibinomial family; available links are %s',
## 			  linktemp, paste(sQuote(okLinks), collapse =", ")),
## 	     domain = NA)
##         }
##     }
##     variance <- function(mu) mu * (1 - mu)
##     dvariance <- function(mu) 1 - 2 * mu
##     ddvariance <- function(mu) -2
##     validmu <- function(mu) all(mu>0) && all(mu<1)
##     dev.resids <- function(y, mu, wt)
## 	2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) +
## 		  (1 - y) * log(ifelse(y == 1, 1, (1 - y)/(1 - mu))))
##     aic <- function(y, n, mu, wt, dev) NA
##     initialize <- expression({
## 	if (NCOL(y) == 1) {
## 	    if (is.factor(y)) y <- y != levels(y)[1L]
## 	    n <- rep.int(1, nobs)
## 	    if (any(y < 0 | y > 1))
## 		stop("y values must be 0 <= y <= 1")
##             mustart <- (weights * y + 0.5)/(weights + 1)
## 	}
## 	else if (NCOL(y) == 2) {
## 	    n <- y[, 1] + y[, 2]
## 	    y <- ifelse(n == 0, 0, y[, 1]/n)
## 	    weights <- weights * n
##             mustart <- (n * y + 0.5)/(n + 1)
## 	}
## 	else stop("for the quasibinomial family, y must be a vector of 0 and 1\'s\n",
##                   "or a 2 column matrix where col 1 is no. successes and col 2 is no. failures")
##     })
##     structure(list(family = "quasibinomial",
## 		   link = linktemp,
## 		   linkfun = stats$linkfun,
## 		   linkinv = stats$linkinv,
## 		   variance = variance,
##                    dvariance = dvariance,
##                    ddvariance = ddvariance,
## 		   dev.resids = dev.resids,
## 		   aic = aic,
## 		   mu.eta = stats$mu.eta,
##                    dmu.deta = stats$dmu.deta,
## 		   initialize = initialize,
## 		   validmu = validmu,
## 		   valideta = stats$valideta),
## 	      class = "family")
## }

## quasi <- function (link = "identity", variance = "constant")
## {
##     linktemp <- substitute(link)
##     if (!is.character(linktemp)) linktemp <- deparse(linktemp)
##     if (linktemp %in% c("logit", "probit", "cloglog", "identity",
##                         "inverse", "log", "1/mu^2", "sqrt"))
##         stats <- make.link(linktemp)
##     else if (is.character(link)) {
##         stats <- make.link(link)
##         linktemp <- link
##     } else {
##         stats <- link
##         linktemp <- if(!is.null(stats$name)) stats$name else deparse(linktemp)
##     }
##     vtemp <- substitute(variance)
##     if (!is.character(vtemp)) vtemp <- deparse(vtemp)
##     variance_nm <- vtemp
##     switch(vtemp,
##            "constant" = {
##                varfun <- function(mu) rep.int(1, length(mu))
##                dev.resids <- function(y, mu, wt) wt * ((y - mu)^2)
##                    validmu <- function(mu) TRUE
##                initialize <- expression({n <- rep.int(1, nobs); mustart <- y})
##            },
##            "mu(1-mu)" = {
##                varfun <- function(mu) mu * (1 - mu)
##                validmu <- function(mu) all(mu>0) && all(mu<1)
##                dev.resids <- function(y, mu, wt)
##                    2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) +
##                              (1 - y) * log(ifelse(y == 1, 1, (1 - y)/(1 - mu))))
##                initialize <- expression({n <- rep.int(1, nobs)
##                                          mustart <- pmax(0.001, pmin(0.999, y))})
##            },
##            "mu" = {
##                varfun <- function(mu) mu
##                validmu <- function(mu) all(mu>0)
##                dev.resids <- function(y, mu, wt)
##                    2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
##                ## 0.1 fudge here matches poisson: S has 1/6.
##                initialize <- expression({n <- rep.int(1, nobs)
##                                          mustart <- y + 0.1 * (y == 0)})
##            },
##            "mu^2" = {
##                varfun <- function(mu) mu^2
##                validmu <- function(mu) all(mu>0)
##                dev.resids <- function(y, mu, wt)
## 		   pmax(-2 * wt * (log(ifelse(y == 0, 1, y)/mu) - (y - mu)/mu), 0)
##                initialize <- expression({n <- rep.int(1, nobs)
##                                          mustart <- y + 0.1 * (y == 0)})
##            },
##            "mu^3" = {
##                varfun <- function(mu) mu^3
##                validmu <- function(mu) all(mu>0)
##                dev.resids <- function(y, mu, wt)
##                    wt * ((y - mu)^2)/(y * mu^2)
##                initialize <- expression({n <- rep.int(1, nobs)
##                                          mustart <- y + 0.1 * (y == 0)})
##            },
##            variance_nm <- NA
##            )# end switch(.)

##     if(is.na(variance_nm)) {
##         if(is.character(variance))
##             stop(gettextf('\'variance\' "%s" is invalid: possible values are "mu(1-mu)", "mu", "mu^2", "mu^3" and "constant"', variance_nm), domain = NA)
##         ## so we really meant the object.
##         varfun <- variance$varfun
##         validmu <- variance$validmu
##         dev.resids <- variance$dev.resids
##         initialize <- variance$initialize
##         variance_nm <- variance$name
##     }
##     aic <- function(y, n, mu, wt, dev) NA
##     structure(list(family = "quasi",
## 		   link = linktemp,
## 		   linkfun = stats$linkfun,
## 		   linkinv = stats$linkinv,
## 		   variance = varfun,
## 		   dev.resids = dev.resids,
## 		   aic = aic,
## 		   mu.eta = stats$mu.eta,
##                    dmu.deta = stats$dmu.deta,
## 		   initialize = initialize,
## 		   validmu = validmu,
## 		   valideta = stats$valideta,
##                    ## character form of the var fun is needed for gee
##                    varfun = variance_nm),
## 	      class = "family")
## }
