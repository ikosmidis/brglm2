#' Defunct Functions in package \pkg{brglm2}
#'
#'  The functions or variables listed here are no longer part of
#'  \pkg{brglm2}.
#'
#'
#' @param ... arguments to be passed to functions and methods.
#'
#'
#' @details
#'
#' * [detect_separation()]: This function is defunct from \pkg{brglm2}
#' since version 0.8.0. A new version of [detect_separation()] is now
#' maintained in the
#' [\pkg{detectseparation}](https://cran.r-project.org/package=detectseparation)
#' R package.
#'
#' * [check_infinite_estimates()] is defunct from \pkg{brglm2} since
#' version 0.8.0. An new version of [check_infinite_estimates()] is
#' now maintained in the
#' [\pkg{detectseparation}](https://cran.r-project.org/package=detectseparation)
#' R package.
#'
#' @name brglm2-defunct
NULL

#' @rdname brglm2-defunct
#' @export
check_infinite_estimates <- function(...) {
    function_moved_to_new_package(gsub("\\(|\\)", "", deparse(match.call()[1])),
                                  "0.8.0",
                                  "brglm2",
                                  "detectseparation")
}

#' @rdname brglm2-defunct
#' @export
detect_separation <- function(...) {
    function_moved_to_new_package(gsub("\\(|\\)", "", deparse(match.call()[1])),
                                  "0.8.0",
                                  "brglm2",
                                  "detectseparation")
}
