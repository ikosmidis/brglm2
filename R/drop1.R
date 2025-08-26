## Essentially stats:::drop1.glm() with `brglmFit()` calls in place of `glm.fit()` calls
#' @export
drop1.brglmFit <- function(object, scope, scale = 0, test = c("none", "Rao", "LRT",  "Chisq", "F"), k = 2, ...) {
    test <- match.arg(test)
    if (test == "Chisq")
        test <- "LRT"
    x <- model.matrix(object)
    n <- nrow(x)
    asgn <- attr(x, "assign")
    tl <- attr(object$terms, "term.labels")
    if (missing(scope))
        scope <- drop.scope(object)
    else {
        if (!is.character(scope))
            scope <- attr(terms(update.formula(object, scope)),
                          "term.labels")
        if (!all(match(scope, tl, 0L) > 0L))
            stop("scope is not a subset of term labels")
    }
    ndrop <- match(scope, tl)
    ns <- length(scope)
    rdf <- object$df.residual
    chisq <- object$deviance
    dfs <- numeric(ns)
    dev <- numeric(ns)
    score <- numeric(ns)
    y <- object$y
    if (is.null(y)) {
        y <- model.response(model.frame(object))
        if (!is.factor(y))
            storage.mode(y) <- "double"
    }
    wt <- object$prior.weights %||% rep.int(1, n)
    for (i in seq_len(ns)) {
        ii <- seq_along(asgn)[asgn == ndrop[i]]
        jj <- setdiff(seq(ncol(x)), ii)
        z <- brglmFit(x[, jj, drop = FALSE], y, wt, offset = object$offset,
                      family = object$family, control = object$control)
        dfs[i] <- z$rank
        dev[i] <- z$deviance
        if (test == "Rao") {
            r <- z$residuals
            w <- z$weights
            zz <- glm.fit(x, r, w)
            score[i] <- zz$null.deviance - zz$deviance
        }
    }
    scope <- c("<none>", scope)
    dfs <- c(object$rank, dfs)
    dev <- c(chisq, dev)
    if (test == "Rao") {
        score <- c(NA, score)
    }
    dispersion <- if (is.null(scale) || scale == 0)
                      summary(object, dispersion = NULL)$dispersion
    else scale
    fam <- object$family$family
    loglik <- if (fam == "gaussian") {
                  if (scale > 0)
                      dev/scale - n
                  else n * log(dev/n)
              }
              else dev/dispersion
    aic <- loglik + k * dfs
    dfs <- dfs[1L] - dfs
    dfs[1L] <- NA
    aic <- aic + (extractAIC(object, k = k)[2L] - aic[1L])
    aod <- data.frame(Df = dfs, Deviance = dev, AIC = aic, row.names = scope,
                      check.names = FALSE)
    if (all(is.na(aic)))
        aod <- aod[, -3]
    if (test == "LRT") {
        dev <- pmax(0, loglik - loglik[1L])
        dev[1L] <- NA
        nas <- !is.na(dev)
        LRT <- if (dispersion == 1)
                   "LRT"
               else "scaled dev."
        aod[, LRT] <- dev
        dev[nas] <- safe_pchisq(dev[nas], aod$Df[nas], lower.tail = FALSE)
        aod[, "Pr(>Chi)"] <- dev
    }
    else if (test == "Rao") {
        dev <- pmax(0, score)
        nas <- !is.na(dev)
        SC <- if (dispersion == 1)
                  "Rao score"
              else "scaled Rao sc."
        dev <- dev/dispersion
        aod[, SC] <- dev
        dev[nas] <- safe_pchisq(dev[nas], aod$Df[nas], lower.tail = FALSE)
        aod[, "Pr(>Chi)"] <- dev
    }
    else if (test == "F") {
        if (fam == "binomial" || fam == "poisson")
            warning(gettextf("F test assumes 'quasi%s' family",
                             fam), domain = NA)
        dev <- aod$Deviance
        rms <- dev[1L]/rdf
        dev <- pmax(0, dev - dev[1L])
        dfs <- aod$Df
        rdf <- object$df.residual
        Fs <- (dev/dfs)/rms
        Fs[dfs < 1e-04] <- NA
        P <- Fs
        nas <- !is.na(Fs)
        P[nas] <- safe_pf(Fs[nas], dfs[nas], rdf, lower.tail = FALSE)
        aod[, c("F value", "Pr(>F)")] <- list(Fs, P)
    }
    head <- c("Single term deletions", "\nModel:", deparse(formula(object)),
              if (!is.null(scale) && scale > 0) paste("\nscale: ",
                                                      format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

#' @export
drop1.mdyplFit <- function(object, scope, scale = 0, test = c("none", "Rao", "LRT",  "Chisq", "F"), k = 2, ...) {
    test <- match.arg(test)
    if (test == "Chisq")
        test <- "LRT"
    x <- model.matrix(object)
    n <- nrow(x)
    asgn <- attr(x, "assign")
    tl <- attr(object$terms, "term.labels")
    if (missing(scope))
        scope <- drop.scope(object)
    else {
        if (!is.character(scope))
            scope <- attr(terms(update.formula(object, scope)),
                          "term.labels")
        if (!all(match(scope, tl, 0L) > 0L))
            stop("scope is not a subset of term labels")
    }
    ndrop <- match(scope, tl)
    ns <- length(scope)
    rdf <- object$df.residual
    chisq <- object$deviance
    dfs <- numeric(ns)
    dev <- numeric(ns)
    score <- numeric(ns)
    y <- object$y
    if (is.null(y)) {
        y <- model.response(model.frame(object))
        if (!is.factor(y))
            storage.mode(y) <- "double"
    }
    wt <- object$prior.weights %||% rep.int(1, n)
    for (i in seq_len(ns)) {
        ii <- seq_along(asgn)[asgn == ndrop[i]]
        jj <- setdiff(seq(ncol(x)), ii)
        z <- mdyplFit(x[, jj, drop = FALSE], y, wt, offset = object$offset,
                      family = object$family, control = object$control)
        dfs[i] <- z$rank
        dev[i] <- z$deviance
        if (test == "Rao") {
            r <- z$residuals
            w <- z$weights
            zz <- glm.fit(x, r, w)
            score[i] <- zz$null.deviance - zz$deviance
        }
    }
    scope <- c("<none>", scope)
    dfs <- c(object$rank, dfs)
    dev <- c(chisq, dev)
    if (test == "Rao") {
        score <- c(NA, score)
    }
    dispersion <- if (is.null(scale) || scale == 0)
                      summary(object, dispersion = NULL)$dispersion
    else scale
    fam <- object$family$family
    loglik <- if (fam == "gaussian") {
                  if (scale > 0)
                      dev/scale - n
                  else n * log(dev/n)
              }
              else dev/dispersion
    aic <- loglik + k * dfs
    dfs <- dfs[1L] - dfs
    dfs[1L] <- NA
    aic <- aic + (extractAIC(object, k = k)[2L] - aic[1L])
    aod <- data.frame(Df = dfs, Deviance = dev, AIC = aic, row.names = scope,
                      check.names = FALSE)
    if (all(is.na(aic)))
        aod <- aod[, -3]
    if (test == "LRT") {
        dev <- pmax(0, loglik - loglik[1L])
        dev[1L] <- NA
        nas <- !is.na(dev)
        LRT <- if (dispersion == 1)
                   "LRT"
               else "scaled dev."
        aod[, LRT] <- dev
        dev[nas] <- safe_pchisq(dev[nas], aod$Df[nas], lower.tail = FALSE)
        aod[, "Pr(>Chi)"] <- dev
    }
    else if (test == "Rao") {
        dev <- pmax(0, score)
        nas <- !is.na(dev)
        SC <- if (dispersion == 1)
                  "Rao score"
              else "scaled Rao sc."
        dev <- dev/dispersion
        aod[, SC] <- dev
        dev[nas] <- safe_pchisq(dev[nas], aod$Df[nas], lower.tail = FALSE)
        aod[, "Pr(>Chi)"] <- dev
    }
    else if (test == "F") {
        if (fam == "binomial" || fam == "poisson")
            warning(gettextf("F test assumes 'quasi%s' family",
                             fam), domain = NA)
        dev <- aod$Deviance
        rms <- dev[1L]/rdf
        dev <- pmax(0, dev - dev[1L])
        dfs <- aod$Df
        rdf <- object$df.residual
        Fs <- (dev/dfs)/rms
        Fs[dfs < 1e-04] <- NA
        P <- Fs
        nas <- !is.na(Fs)
        P[nas] <- safe_pf(Fs[nas], dfs[nas], rdf, lower.tail = FALSE)
        aod[, c("F value", "Pr(>F)")] <- list(Fs, P)
    }
    head <- c("Single term deletions", "\nModel:", deparse(formula(object)),
              if (!is.null(scale) && scale > 0) paste("\nscale: ",
                                                      format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}
