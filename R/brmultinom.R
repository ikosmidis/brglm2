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
    mf <- eval.parent(mf)
    Terms <- attr(mf, "terms")
    X <- model.matrix(Terms, mf, contrasts)
    Xcontrasts <- attr(X, "contrasts")
    Y <- model.response(mf, "any")

    ## The chunk of code between +BEGIN and +END has been adopted from
    ## nnet::multinom
    ##+BEGIN
    if (!is.matrix(Y))
        Y <- as.factor(Y)
    w <- model.weights(mf)
    if (length(w) == 0L)
        if (is.matrix(Y))
            w <- rep(1, dim(Y)[1L])
        else w <- rep(1, length(Y))
    lev <- levels(Y)
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
        ## Consider removing the below
        if (length(lev) < 2L)
            stop("need two or more classes to fit a multinom model")
        if (length(lev) == 2L)
            Y <- as.integer(Y) - 1
        else Y <- nnet::class.ind(Y)
    }
    ##+END

    ncat <- length(lev)
    nvar <- ncol(X)

    ## Set up the model matrix for the poisson fit
    Xnuisance <- Matrix::Diagonal(nrow(X))
    Xextended <- cbind(kronecker(rep(1, ncat), Xnuisance),
                       kronecker(Matrix::Diagonal(ncat)[, -1], X))
    int <- seq.int(nrow(X))
    nd <- paste0("%0", nchar(max(int)), "d")
    colnames(Xextended) <- c(paste0(".nuisance", sprintf(nd, int)),
                             ## CHECK: lev[-1] contrasts?
                             paste(rep(lev[-1], each = nvar),
                                   rep(colnames(X), ncat - 1), sep = ":"))

    ## Set up the extended response
    Yextended <- c(Y * w)

    ## TODO:
    ## + starting values
    ## + subset
    ## + na.action
    ## + control
    fit <- brglmFit(x = Xextended, y = Yextended, start = NULL,
            family = poisson("log"))


    ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@3@"]]));##:ess-bp-end:##

    1




}



