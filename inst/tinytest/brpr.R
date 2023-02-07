## Bias-reduced Poisson and multinomial log-linear models
## using maximum penalized likelihood as in Firth (1993) Biometrika.
##
## Author: David Firth, d.firth@warwick.ac.uk, 2 Sep 2010
##
## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE!  Provided "as is",
## licensed under GPL2.  NO WARRANTY OF FITNESS FOR ANY PURPOSE!
##
## This code still has lots of rough edges, some of which are
## indicated by the embedded comments.
##
## The intention is to include a more polished implementation as part
## of the CRAN package "brglm", in due course.

## Arguments are as for a Poisson log-linear glm,
## with one new argument:
## -- "fixed.totals", a factor indexing groups of observations
##    whose response totals are fixed (for example, this might be
##    the rows of a table).
## If fixed.totals is NULL (the default), pure Poisson sampling (with
## no totals fixed) is assumed.

brpr <- function(formula,
                  data,
                  subset,
                  na.action,
                  offset = NULL,
                  control = glm.control(...),
                  contrasts = NULL,
                  fixed.totals = NULL,
                  weights = NULL,
                  ...) {
    ## The "weights" and "offset" arguments here can only be NULL -- anything
    ## else is an error (for now, at least)
    if (!is.null(weights)) stop("The weights argument here can only be NULL")
    if (!is.null(offset)) stop("The offset argument here can only be NULL")
    ##  work needed here, to allow the use of offset as an argument?

    ## Initial setup as in glm (except that `prior weights' are not used here,
    ## and an offset (if any) must be specified through the formula rather
    ## than (as is also allowed by glm) as a separate argument:
    call <- match.call()
    theFormula <- formula
    if (missing(data))
        data <- environment(formula)
    if (is.null(call$epsilon)) control$epsilon <- 1e-8
    ##  Is 1e-8 too stringent?  Work needed here!
    fixed.totals <- eval(substitute(fixed.totals), envir = data)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action",
                 "offset", "fixed.totals"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "numeric")

    ## Correct the fixed.totals variable for subsetting, NA treatment, etc.
    if (!is.null(fixed.totals)) {
        fixed.totals <-  mf$"(fixed.totals)"
        if (!is.factor(fixed.totals)) stop("fixed.totals must be a factor")
        rowTotals <-  as.vector(tapply(Y, fixed.totals, sum))[fixed.totals]
    }

    ## Do an initial fit of the model, with constant-adjusted counts
    adj <- 0.5
    X <- model.matrix(formula, data = mf)
    offset <- model.offset(mf)
    nFixed <- if (is.null(fixed.totals)) 0 else nlevels(fixed.totals)
    mf$y.adj <- Y + adj * (ncol(X) - nFixed)/nrow(X)
    formula <- update(formula, y.adj ~ .)
    op <- options(warn = -1)
    fit <- glm.fit(X, mf$y.adj, family = poisson(), offset = offset,
                   control = glm.control(maxit = 1))
    options(op)
    ## Update the model iteratively, refining the adjustment at each iteration
    criterion <- 1
    iter <- 1
    coefs <- coef(fit)
    Xwork <- X[ , !is.na(coefs), drop = FALSE]
    muAdj <- fitted(fit)
    coefs <- na.omit(coefs)
    while (criterion > control$epsilon && iter < control$maxit) {
        iter <- iter + 1
        if (!is.null(fixed.totals)){
            muTotals <-  as.vector(tapply(muAdj, fixed.totals,
                                          sum))[fixed.totals]
            mu <- muAdj * rowTotals / muTotals
        } else { ## case fixed.totals is NULL
            mu <- muAdj
        }
        W.X <- sqrt(mu) * Xwork
        XWXinv <- chol2inv(chol(crossprod(W.X)))
        coef.se <- sqrt(diag(XWXinv))
        h <- diag(Xwork %*% XWXinv %*% t(mu * Xwork))
        mf$y.adj <- Y + h * adj
        z <- log(muAdj) + (mf$y.adj - muAdj)/muAdj
        lsfit <- lm.wfit(Xwork, z, muAdj, offset = offset)
        criterion <- max((abs(lsfit$coefficients - coefs))/coef.se)
        if (control$trace) cat("Iteration ", iter,
                               ": largest (scaled) coefficient change is ",
                               criterion, "\n", sep = "")
        coefs <- lsfit$coefficients
        muAdj <- exp(lsfit$fitted.values)
    }

    ## The object `lsfit' uses *adjusted* counts, fitted values, etc. --
    ## so we need to correct various things before returning (so that, for example,
    ## deviance and standard errors are correct).
    fit$coefficients[!is.na(fit$coefficients)] <- coefs
    fit$y <- Y
    fit$fitted.values <- fit$weights <- mu
    fit$linear.predictors <- log(mu)
    dev.resids <- poisson()$dev.resids
    wt <- fit$prior.weights
    fit$deviance <- sum(dev.resids(Y, mu, wt))
    fit$null.deviance <- sum(dev.resids(Y, mean(Y), wt))
    fit$aic <- NA # sum((poisson()$aic)(Y, 1, mu, wt, 0)) ??  work needed!
    fit$residuals <- Y/mu - 1

    ## Next bit is straight from glm -- not sure it applies here (work needed!)
    if (any(offset) && attr(mt, "intercept") > 0) {
        fit$null.deviance <- glm.fit(x = X[, "(Intercept)", drop = FALSE],
                                     y = Y, weights = rep(1, nrow(mf)),
                                     offset = offset, family = poisson(),
                                     intercept = TRUE)$deviance
    }

    fit$iter <- iter
    fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    fit$x <- X
    fit$qr <- qr(sqrt(mu) * X)
    fit$R <- qr.R(fit$qr, complete = TRUE)
    ## The "effects" component of the fit object is almost certainly not
    ## correct as is -- but where does it actually get used?  (work needed!)
    ##  rownames(fit$R) <- colnames(fit$R)  ##  FIX THIS!
    fit <- c(fit, list(call = call, formula = theFormula, terms = mt,
                       data = data, offset = offset, method = "brpr",
                       contrasts = attr(X, "contrasts"),
                       xlevels = .getXlevels(mt, mf)))
    class(fit) <- c("glm", "lm")
    return(fit)
}
