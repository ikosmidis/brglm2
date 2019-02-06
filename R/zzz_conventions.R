## Some convenience names

#' @rdname brglmControl
#' @export
brglm_control <- brglmControl

#' @rdname brglmFit
#' @export
brglm_fit <- brglmFit

#' @rdname detect_separation
#' @export
detectSeparation <- detect_separation

#' @rdname detect_separation_control
#' @export
detectSeparationControl <- detect_separation_control

#' @export
checkInfiniteEstimates <- check_infinite_estimates

#### Method conventions

#' @rdname residuals.brmultinom
#' @method residuals bracl
#' @export
residuals.bracl <- residuals.brmultinom

