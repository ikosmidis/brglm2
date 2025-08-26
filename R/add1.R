## Essentially stats:::add1.glm() with `brglmFit()` calls in place of `glm.fit()` calls
#' @export
add1.brglmFit <- function(object, scope, scale = 0,
                          test = c("none", "Rao", "LRT", "Chisq", "F"), x = NULL, k = 2, ...)  {

    Fstat <- function(table, rdf) {
        dev <- table$Deviance
        df <- table$Df
        diff <- pmax(0, (dev[1L] - dev)/df)
        Fs <- diff/(dev/(rdf - df))
        Fs[df < .Machine$double.eps] <- NA
        P <- Fs
        nnas <- !is.na(Fs)
        P[nnas] <- safe_pf(Fs[nnas], df[nnas], rdf - df[nnas],
                           lower.tail = FALSE)
        list(Fs = Fs, P = P)
    }
    test <- match.arg(test)
    if (test == "Chisq")
        test <- "LRT"
    if (!is.character(scope))
        scope <- add.scope(object, update.formula(object, scope))
    if (!length(scope))
        stop("no terms in scope for adding to object")
    oTerms <- attr(object$terms, "term.labels")
    int <- attr(object$terms, "intercept")
    ns <- length(scope)
    dfs <- numeric(ns + 1L)
    names(dfs) <- c("<none>", scope)
    dev <- score <- dfs
    add.rhs <- eval(str2lang(paste("~ . +", paste(scope, collapse = "+"))))
    new.form <- update.formula(object, add.rhs)
    Terms <- terms(new.form)
    y <- object$y
    if (is.null(x)) {
        fc <- object$call
        fc$formula <- Terms
        fob <- list(call = fc, terms = Terms)
        class(fob) <- oldClass(object)
        m <- model.frame(fob, xlev = object$xlevels)
        offset <- model.offset(m)
        wt <- model.weights(m)
        x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
        oldn <- length(y)
        y <- model.response(m)
        if (!is.factor(y))
            storage.mode(y) <- "double"
        if (NCOL(y) == 2) {
            n <- y[, 1] + y[, 2]
            y <- ifelse(n == 0, 0, y[, 1]/n)
            wt <- (wt %||% rep.int(1, length(y))) * n
        }
        newn <- length(y)
        if (newn < oldn)
            warning(sprintf(ngettext(newn, "using the %d/%d row from a combined fit",
                                     "using the %d/%d rows from a combined fit"),
                            newn, oldn), domain = NA)
    }
    else {
        wt <- object$prior.weights
        offset <- object$offset
    }
    n <- nrow(x)
    if (is.null(wt))
        wt <- rep.int(1, n)
    Terms <- attr(Terms, "term.labels")
    asgn <- attr(x, "assign")
    ousex <- match(asgn, match(oTerms, Terms), 0L) > 0L
    if (int)
        ousex[1L] <- TRUE
    X <- x[, ousex, drop = FALSE]
    z <- brglmFit(X, y, wt, offset = offset, family = object$family,
                  control = object$control)
    dfs[1L] <- z$rank
    dev[1L] <- z$deviance
    r <- z$residuals
    w <- z$weights
    sTerms <- sapply(strsplit(Terms, ":", fixed = TRUE), function(x) paste(sort(x),
                                                                           collapse = ":"))
    for (tt in scope) {
        stt <- paste(sort(strsplit(tt, ":")[[1L]]), collapse = ":")
        usex <- match(asgn, match(stt, sTerms), 0L) > 0L
        X <- x[, usex | ousex, drop = FALSE]
        z <- brglmFit(X, y, wt, offset = offset, family = object$family,
                      control = object$control)
        dfs[tt] <- z$rank
        dev[tt] <- z$deviance
        if (test == "Rao") {
            zz <- glm.fit(X, r, w, offset = offset)
            score[tt] <- zz$null.deviance - zz$deviance
        }
    }
    dispersion <- if (scale == 0)
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
    aic <- aic + (extractAIC(object, k = k)[2L] - aic[1L])
    dfs <- dfs - dfs[1L]
    dfs[1L] <- NA
    aod <- data.frame(Df = dfs, Deviance = dev, AIC = aic, row.names = names(dfs),
                      check.names = FALSE)
    if (all(is.na(aic)))
        aod <- aod[, -3]
    test <- match.arg(test)
    if (test == "LRT") {
        dev <- pmax(0, loglik[1L] - loglik)
        dev[1L] <- NA
        LRT <- if (dispersion == 1)
                   "LRT"
               else "scaled dev."
        aod[, LRT] <- dev
        nas <- !is.na(dev)
        dev[nas] <- safe_pchisq(dev[nas], aod$Df[nas], lower.tail = FALSE)
        aod[, "Pr(>Chi)"] <- dev
    }
    else if (test == "Rao") {
        dev <- pmax(0, score)
        dev[1L] <- NA
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
            warning(gettextf("F test assumes quasi%s family",
                             fam), domain = NA)
        rdf <- object$df.residual
        aod[, c("F value", "Pr(>F)")] <- Fstat(aod, rdf)
    }
    head <- c("Single term additions", "\nModel:", deparse(formula(object)),
              if (scale > 0) paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

#' @export
add1.mdyplFit <- function(object, scope, scale = 0,
                          test = c("none", "Rao", "LRT", "Chisq", "F"), x = NULL, k = 2, ...)  {

    Fstat <- function(table, rdf) {
        dev <- table$Deviance
        df <- table$Df
        diff <- pmax(0, (dev[1L] - dev)/df)
        Fs <- diff/(dev/(rdf - df))
        Fs[df < .Machine$double.eps] <- NA
        P <- Fs
        nnas <- !is.na(Fs)
        P[nnas] <- safe_pf(Fs[nnas], df[nnas], rdf - df[nnas],
                           lower.tail = FALSE)
        list(Fs = Fs, P = P)
    }
    test <- match.arg(test)
    if (test == "Chisq")
        test <- "LRT"
    if (!is.character(scope))
        scope <- add.scope(object, update.formula(object, scope))
    if (!length(scope))
        stop("no terms in scope for adding to object")
    oTerms <- attr(object$terms, "term.labels")
    int <- attr(object$terms, "intercept")
    ns <- length(scope)
    dfs <- numeric(ns + 1L)
    names(dfs) <- c("<none>", scope)
    dev <- score <- dfs
    add.rhs <- eval(str2lang(paste("~ . +", paste(scope, collapse = "+"))))
    new.form <- update.formula(object, add.rhs)
    Terms <- terms(new.form)
    y <- object$y
    if (is.null(x)) {
        fc <- object$call
        fc$formula <- Terms
        fob <- list(call = fc, terms = Terms)
        class(fob) <- oldClass(object)
        m <- model.frame(fob, xlev = object$xlevels)
        offset <- model.offset(m)
        wt <- model.weights(m)
        x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
        oldn <- length(y)
        y <- model.response(m)
        if (!is.factor(y))
            storage.mode(y) <- "double"
        if (NCOL(y) == 2) {
            n <- y[, 1] + y[, 2]
            y <- ifelse(n == 0, 0, y[, 1]/n)
            wt <- (wt %||% rep.int(1, length(y))) * n
        }
        newn <- length(y)
        if (newn < oldn)
            warning(sprintf(ngettext(newn, "using the %d/%d row from a combined fit",
                                     "using the %d/%d rows from a combined fit"),
                            newn, oldn), domain = NA)
    }
    else {
        wt <- object$prior.weights
        offset <- object$offset
    }
    n <- nrow(x)
    if (is.null(wt))
        wt <- rep.int(1, n)
    Terms <- attr(Terms, "term.labels")
    asgn <- attr(x, "assign")
    ousex <- match(asgn, match(oTerms, Terms), 0L) > 0L
    if (int)
        ousex[1L] <- TRUE
    X <- x[, ousex, drop = FALSE]
    z <- mdyplFit(X, y, wt, offset = offset, family = object$family,
                  control = object$control)
    dfs[1L] <- z$rank
    dev[1L] <- z$deviance
    r <- z$residuals
    w <- z$weights
    sTerms <- sapply(strsplit(Terms, ":", fixed = TRUE), function(x) paste(sort(x),
                                                                           collapse = ":"))
    for (tt in scope) {
        stt <- paste(sort(strsplit(tt, ":")[[1L]]), collapse = ":")
        usex <- match(asgn, match(stt, sTerms), 0L) > 0L
        X <- x[, usex | ousex, drop = FALSE]
        z <- mdyplFit(X, y, wt, offset = offset, family = object$family,
                      control = object$control)
        dfs[tt] <- z$rank
        dev[tt] <- z$deviance
        if (test == "Rao") {
            zz <- glm.fit(X, r, w, offset = offset)
            score[tt] <- zz$null.deviance - zz$deviance
        }
    }
    dispersion <- if (scale == 0)
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
    aic <- aic + (extractAIC(object, k = k)[2L] - aic[1L])
    dfs <- dfs - dfs[1L]
    dfs[1L] <- NA
    aod <- data.frame(Df = dfs, Deviance = dev, AIC = aic, row.names = names(dfs),
                      check.names = FALSE)
    if (all(is.na(aic)))
        aod <- aod[, -3]
    test <- match.arg(test)
    if (test == "LRT") {
        dev <- pmax(0, loglik[1L] - loglik)
        dev[1L] <- NA
        LRT <- if (dispersion == 1)
                   "LRT"
               else "scaled dev."
        aod[, LRT] <- dev
        nas <- !is.na(dev)
        dev[nas] <- safe_pchisq(dev[nas], aod$Df[nas], lower.tail = FALSE)
        aod[, "Pr(>Chi)"] <- dev
    }
    else if (test == "Rao") {
        dev <- pmax(0, score)
        dev[1L] <- NA
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
            warning(gettextf("F test assumes quasi%s family",
                             fam), domain = NA)
        rdf <- object$df.residual
        aod[, c("F value", "Pr(>F)")] <- Fstat(aod, rdf)
    }
    head <- c("Single term additions", "\nModel:", deparse(formula(object)),
              if (scale > 0) paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}


safe_pf <- utils::getFromNamespace("safe_pf", "stats")
safe_pchisq <- utils::getFromNamespace("safe_pchisq", "stats")
