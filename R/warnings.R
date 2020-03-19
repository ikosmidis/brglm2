function_moves_to_new_package <- function(function_name,
                                          removal_version,
                                          current_pkg,
                                          new_pkg,
                                          extra_message = NULL) {
    function_name <- paste0("'", function_name, "'")
    current_pkg <- paste0("'", current_pkg, "'")
    new_pkg <- paste0("'", new_pkg, "'")
    msg <- paste(function_name,
                 "will be removed from",
                 current_pkg,
                 "at version",
                 paste0(removal_version, "."),
                 "A new version of",
                 function_name,
                 "is now maintained in the",
                 new_pkg,
                 "package.",
                 extra_message)
    .Deprecated(msg = msg, package = new_pkg)
}
