#' Defunct Functions in package \pkg{brglm2}
#'
#'  The functions or variables listed here are no longer part of
#'  \pkg{brglm2}.
#'
#' \itemize{
#'
#' \item \code{\link{detect_separation}}: This function is defunct
#' from \pkg{brglm2} since version 0.8.0. A new version of
#' \code{detect_separation} is now maintained in the
#' \pkg{detectseparation} R package at
#' \url{https://cran.r-project.org/package=detectseparation}.
#'
#' \item \code{check_infinite_estimates} is defunct from
#' \pkg{brglm2} since version 0.8.0. An new version of
#' \code{check_infinite_estimates} is now maintained in the
#' \pkg{detectseparation} R package at
#' \url{https://cran.r-project.org/package=detectseparation}.
#'
#' }
#' @name brglm2-defunct
NULL

#' @rdname brglm2-defunct
#' @export
check_infinite_estimates <- function(object, ...) {
    function_moved_to_new_package(gsub("\\(|\\)", "", deparse(match.call()[1])),
                                  "0.8.0",
                                  "brglm2",
                                  "detectseparation")
}

#' @rdname brglm2-defunct
#' @export
detect_separation <- function(x, y, weights = rep(1, nobs),
                             start = NULL, etastart = NULL,  mustart = NULL,
                             offset = rep(0, nobs), family = gaussian(),
                             control = list(), intercept = TRUE, singular.ok = TRUE) {
    function_moved_to_new_package(gsub("\\(|\\)", "", deparse(match.call()[1])),
                                  "0.8.0",
                                  "brglm2",
                                  "detectseparation")
}
