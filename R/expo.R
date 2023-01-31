#' @export
expo.brglmFit <- function(object, type = c("correction*", "correction+", "Lylesetal2012", "AS_median", "ML"), level = 0.95) {
    type <- match.arg(type)
    object_type <- object$type
    to_correct <- type %in% c("correction*", "correction+", "Lylesetal2012")
    if (to_correct) {
        if (object_type != "AS_mixed") {
            object <- update(object, type = "AS_mixed", data = object$data)
        }
    } else {
        ## AS_mean and ML are invariant to transformation
        if (object_type != type) {
            object <- update(object, type = type, data = object$data)
        }
    }
    info_fun <- get_information_function(object)
    dispersion <- object$dispersion
    coefs <- coef(object, model = "mean")
    ses <- sqrt(diag(solve(info_fun(coefs, dispersion))))
    trans_coefs <- exp(coefs) * switch(type,
                                       "ML" = 1,
                                       "AS_median" = 1,
                                       "Lylesetal2012" = exp(- ses^2/2),
                                       "correction*" = 1 / (1 + ses^2/2),
                                       "correction+" = (1 - ses^2/2))
    if (to_correct) {
        ## Recompute coefs and ses on the original scale from trans_coefs if method is not invariant
        coefs <- log(trans_coefs)
        ses <- sqrt(diag(solve(info_fun(coefs, dispersion))))
    }
    trans_ses <- trans_coefs * ses
    a <- (1 - level)/2
    ci <- coefs + qnorm(a) * cbind(ses, -ses)
    colnames(ci) <- paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE, digits = 3), "%")
    out <- list(coef = trans_coefs, se = trans_ses, ci = exp(ci), type = type)
    out$call <- match.call()
    class(out) <- c("brglmFit_expo", class(out))
    out
}

#' @method print brglmFit_expo
#' @export
print.brglmFit_expo <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cmat <- cbind(x$coef, x$se, x$ci)
    colnames(cmat) <- c("Estimate", "Std. Error", colnames(x$ci))
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    printCoefmat(cmat, digits = digits, ...)
    cat("\n\nType of estimator:", x$type, get_type_description(x$type), "\n")
}

#' @method coef brglmFit_expo
#' @export
coef.brglmFit_expo <- function(object, ...) {
    object$coef
}
