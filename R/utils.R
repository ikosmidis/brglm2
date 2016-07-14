unless_null <- function(x, if_null) {
    if (is.null(x)) {
        if_null
    }
    else {
        x
    }
}

