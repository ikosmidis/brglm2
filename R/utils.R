#' A simple utility that checks whether its argument x is null and
#' returns if.null if it is and x if it is not
unless_null <- function(x, if_null) {
    if (is.null(x)) {
        if_null
    }
    else {
        x
    }
}

