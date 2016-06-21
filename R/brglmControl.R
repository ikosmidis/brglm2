brglmControl <- function(epsilon = 1e-07, maxit = 100,
                         trace = FALSE,
                         qr = TRUE,
                         type = c("adjusted_scores", "correction", "maximum_likelihood"),
                         dispTrans = "identity",
                         slowit = 1,
                         maxStepFactor = 1) {
    type <- match.arg(type)
    Trans <- switch(dispTrans,
                    identity = expression(disp),
                    sqrt = expression(disp^0.5),
                    inverse = expression(1/disp),
                    log = expression(log(disp)),
                    inverseSqrt = expression(1/sqrt(disp)),
                    custom = customTrans[[1]])
    inverseTrans <- switch(dispTrans,
                           identity = expression(transdisp),
                           sqrt = expression(transdisp^2),
                           inverse = expression(1/transdisp),
                           log = expression(exp(transdisp)),
                           inverseSqrt = expression(1/transdisp^2),
                           custom = customTrans[[2]])
    if (!(dispTrans %in% c("identity",
                           "sqrt",
                           "inverse",
                           "log",
                           "custom",
                           "inverseSqrt"))) {
    stop(dispTrans, " is not a supported transofrmation of the dispersion")
    }
    ## if (dispTrans == "custom") {
    ##     warning("Please note that bo check has been made whether Trans and inverseTrans agree")
    ## }
    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("value of 'epsilon' must be > 0")
    if (!is.numeric(maxit) || maxit <= 0)
        stop("maximum number of iterations must be > 0")
    list(epsilon = epsilon, maxit = maxit, trace = trace,
       type = type,
       Trans = Trans,
       inverseTrans = inverseTrans,
       dispTrans = dispTrans,
       slowit = slowit,
       maxStepFactor = maxStepFactor,
       qr = qr)
}


