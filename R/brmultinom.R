#' Bias reduction for multinomial response models using the
#' "Poisson trick". See Kosmidis & Firth (2011) for details
#'
#' @details
#'
#' This function is a wrapper for brglmFit.R that can be used to get
#' bias reduced estimates for the parameters of multinomial regression
#' models.
#'
brmultinom <- function(formula, data, weights, subset, na.action, contrasts = NULL, control = list(...), ...) {
    call <- match.call()
    if (missing(data)) {
        data <- environment(formula)
    }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")

    ## Exceptions
    stopifnot(is.factor(Y))


##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@15@"]]));##:ess-bp-end:##

    ## From multinom
    if (is.factor(Y)) {
        counts <- table(Y)
        if (any(counts == 0L)) {
            empty <- lev[counts == 0L]
            warning(sprintf(ngettext(length(empty), "group %s is empty",
                                     "groups %s are empty"), paste(sQuote(empty),
                                                                   collapse = " ")), domain = NA)
            Y <- factor(Y, levels = lev[counts > 0L])
            lev <- lev[counts > 0L]
        }
        if (length(lev) < 2L)
            stop("need two or more classes to fit a multinom model")
        if (length(lev) == 2L)
            Y <- as.integer(Y) - 1
        else Y <- class.ind(Y)
    }

    if (is.factor(Y))

    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0L)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")

}



