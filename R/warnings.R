# Copyright (C) 2020- Ioannis Kosmidis

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
    .Deprecated(msg = msg, package = current_pkg)
}

function_moved_to_new_package <- function(function_name,
                                          removal_version,
                                          current_pkg,
                                          new_pkg,
                                          extra_message = NULL) {
    function_name <- paste0("'", function_name, "'")
    current_pkg <- paste0("'", current_pkg, "'")
    new_pkg <- paste0("'", new_pkg, "'")
    msg <- paste(function_name,
                 "has been removed from",
                 current_pkg,
                 "at version",
                 paste0(removal_version, "."),
                 "A new version of",
                 function_name,
                 "is now maintained in the",
                 new_pkg,
                 "package.",
                 extra_message)
    .Defunct(msg = msg, package = current_pkg)
}
