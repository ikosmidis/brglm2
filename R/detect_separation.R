# Copyright (C) 2017-2020 Ioannis Kosmidis

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


#' Method for \code{\link{glm}} that tests for data separation and
#' finds which parameters have infinite maximum likelihood estimates
#' in generalized linear models with binomial responses
#'
#' \code{\link{detect_separation}} is a method for \code{\link{glm}}
#' that tests for the occurrence of complete or quasi-complete
#' separation in datasets for binomial response generalized linear
#' models, and finds which of the parameters will have infinite
#' maximum likelihood estimates. \code{\link{detect_separation}}
#' relies on the linear programming methods developed in Konis (2007).
#'
#' @inheritParams stats::glm.fit
#'
#' @aliases detectSeparation print.detect_separation
#'
#' @param x \code{x} is a design matrix of dimension \code{n * p}.
#' @param y \code{y} is a vector of observations of length \code{n}.
#' @param control a list of parameters controlling separation
#'     detection. See \code{\link{detect_separation_control}} for
#'     details.
#' @param start currently not used.
#' @param mustart currently not used.
#' @param etastart currently not used.
#' @param singular.ok logical. If \code{FALSE}, a singular model is an
#'     error.
#'
#' @details
#'
#' For the definition of complete and quasi-complete separation, see
#' Albert and Anderson (1984).
#'
#' \code{\link{detect_separation}} is a wrapper to the \code{separator}
#' function from the **safeBinaryRegression** R package, that can be
#' passed directly as a method to the \code{\link{glm}} function. See,
#' examples.
#'
#' The interface to \code{separator} was designed by Ioannis Kosmidis
#' after correspondence with Kjell Konis, and a port of
#' \code{separator} has been included in **brglm2** under the
#' permission of Kjell Konis.
#'
#' \code{detectSeparation} is an alias for \code{detect_separation}.
#'
#' @note
#' 
#' \code{detect_separation} will be removed from \pkg{brglm2} at
#' version 0.8. A new version of \code{detect_separation} is now
#' maintained in the \pkg{detectseparation} R package at
#' \url{https://cran.r-project.org/package=detectseparation}. In order
#' to use the version in \code{detect_separation} load first
#' \pkg{brglm2} and then \pkg{detectseparation}, i.e. 
#' \code{library(brglm2); library(detectseparation)}.
#'
#' @author Ioannis Kosmidis [aut, cre] \email{ioannis.kosmidis@warwick.ac.uk}, Kjell Konis [ctb] \email{kjell.konis@me.com}
#'
#' @seealso \code{\link{brglm_fit}}, \code{\link{glm.fit}} and \code{\link{glm}}
#'
#' @references
#'
#' Konis K. (2007). *Linear Programming Algorithms for Detecting
#' Separated Data in Binary Logistic Regression
#' Models*. DPhil. University of Oxford.
#' \url{https://ora.ox.ac.uk/objects/uuid:8f9ee0d0-d78e-4101-9ab4-f9cbceed2a2a}
#'
#' Konis K. (2013). safeBinaryRegression: Safe Binary Regression. R
#' package version 0.1-3.
#' \url{https://CRAN.R-project.org/package=safeBinaryRegression}
#'
#' @examples
#'
#' ## endometrial data from Heinze \& Schemper (2002) (see ?endometrial)
#' data("endometrial", package = "brglm2")
#' endometrial_sep <- glm(HG ~ NV + PI + EH, data = endometrial,
#'                        family = binomial("logit"),
#'                        method = "detect_separation")
#' endometrial_sep
#' ## The maximum likelihood estimate for NV is infinite
#' summary(update(endometrial_sep, method = "glm.fit"))
#'
#' \dontrun{
#' ## Example inspired by unpublished microeconometrics lecture notes by
#' ## Achim Zeileis https://eeecon.uibk.ac.at/~zeileis/
#' ## The maximum likelihood estimate of sourhernyes is infinite
#' data("MurderRates", package = "AER")
#' murder_sep <- glm(I(executions > 0) ~ time + income +
#'                   noncauc + lfp + southern, data = MurderRates,
#'                   family = binomial(), method = "detect_separation")
#' murder_sep
#' ## which is also evident by the large estimated standard error for NV
#' murder_glm <- update(murder_sep, method = "glm.fit")
#' summary(murder_glm)
#' ## and is also reveal by the divergence of the NV column of the
#' ## result from the more computationally intensive check
#' check_infinite_estimates(murder_glm)
#' ## Mean bias reduction via adjusted scores results in finite estimates
#' update(murder_glm, method = "brglm_fit")
#' }
#' @export
detect_separation <- function(x, y, weights = rep(1, nobs),
                             start = NULL, etastart = NULL,  mustart = NULL,
                             offset = rep(0, nobs), family = gaussian(),
                             control = list(), intercept = TRUE, singular.ok = TRUE) {

    function_moves_to_new_package(gsub("\\(|\\)", "", deparse(match.call()[1])),
                                  "0.8",
                                  "brglm2",
                                  "detectseparation")
    
    if (family$family != "binomial") {
        warning("detect_separation has been developed for use with binomial-response models")
    }
    control <- do.call("detect_separation_control", control)
    ## ensure x is a matrix
    x <- as.matrix(x)
    betas_names <- dimnames(x)[[2L]]
    ##
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) {
        weights <- rep.int(1, nobs)
    }
    if (missingOffset <- is.null(offset)) {
        offset <- rep.int(0, nobs)
    }
    ## Initialize as prescribed in family
    eval(family$initialize)
    if (EMPTY) {
        out <- list(separation = FALSE)
    }
    else {
        ## as in brglmFit
        boundary <- converged <- FALSE
        ## Detect aliasing
        qrx <- qr(x)
        rank <- qrx$rank
        is_full_rank <- rank == nvars

        if (!singular.ok && !is_full_rank) {
            stop("singular fit encountered")
        }

        if (!isTRUE(is_full_rank)) {
            aliased <- qrx$pivot[seq.int(qrx$rank + 1, nvars)]
            X_all <- x
            x <- x[, -aliased]
            nvars_all <- nvars
            nvars <- ncol(x)
            betas_names_all <- betas_names
            betas_names <- betas_names[-aliased]
        }
        else {
            nvars_all <- nvars
            betas_names_all <- betas_names
        }
        betas_all <- structure(rep(NA_real_, nvars_all), .Names = betas_names_all)
        ## Observations with zero weight do not enter calculations so ignore
        keep <- weights > 0
        x <- x[keep, , drop = FALSE]
        y <- y[keep]
        ## Reshape data set: keep 0 and 1, and replace anything in (0,
        ## 1) with one zero and one 1
        ones <- y == 1
        zeros <- y == 0
        non_boundary <- !(ones | zeros)        
        x <- x[c(which(ones), which(zeros), rep(which(non_boundary), 2)), , drop = FALSE]
        y <- c(y[ones], y[zeros], rep(c(0., 1.), each = sum(non_boundary)))
        ## Run linear program
        out <- separator(x = x, y = y, linear_program = control$linear_program, purpose = control$purpose, beta_tolerance = control$beta_tolerance)
        if (is.null(out$beta)) {
            betas_all <- NULL
        }
        else {
            betas <- out$beta
            names(betas) <- betas_names
            inds <- abs(betas) < control$beta_tolerance
            betas <- Inf * betas
            betas[inds] <- 0
            betas_all[betas_names] <- betas
        }
        out <- list(x = x, y = y, betas = betas_all, separation = out$separation)
    }
    out$linear_program <- control$linear_program
    out$purpose <- control$purpose
    out$class <- "detect_separation"
    class(out) <- "detect_separation_core"
    return(out)
}

#' Auxiliary function for the \code{\link{glm}} interface when
#' \code{method} is \code{\link{detect_separation}}.
#'
#' Typically only used internally by \code{\link{detect_separation}}
#' but may be used to construct a \code{control} argument.
#'
#' @aliases detectSeparationControl
#' @param linear_program should \code{\link{detect_separation}} solve
#'     the \code{"primal"} or \code{"dual"} linear program for
#'     separation detection?
#' @param purpose should \code{\link{detect_separation}} simply
#'     \code{"test"} for separation or also \code{"find"} which
#'     parameters are infinite?
#' @param beta_tolerance maximum absolute variable value from the
#'     linear program, before separation is declared.
#'
#' @export
detect_separation_control <- function(linear_program = c("primal", "dual"),
                                      purpose = c("find", "test"),
                                      beta_tolerance = sqrt(.Machine$double.eps)) {
    linear_program <- match.arg(linear_program)
    purpose <- match.arg(purpose)
    list(linear_program = linear_program, purpose = purpose, beta_tolerance = beta_tolerance)
}


#' @method print detect_separation
#' @export
print.detect_separation <- function(x, digits = max(5L, getOption("digits") - 3L), ...) {
    cat("Separation:", x$separation, "\n")
    if (!is.null(x$betas)) {
        cat("Existence of maximum likelihood estimates\n")
        print(x$betas)
        cat("0: finite value, Inf: infinity, -Inf: -infinity\n")
    }
}

print.detect_separation_core <- function(x, digits = max(5L, getOption("digits") - 3L), ...) {
    cat("Separation:", x$separation, "\n")
    if (!is.null(x$betas)) {
        cat("Existence of maximum likelihood estimates\n")
        print(x$betas)
        cat("0: finite value, Inf: infinity, -Inf: -infinity\n")
    }
}

