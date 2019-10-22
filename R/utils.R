# Copyright (C) 2016-2019 Ioannis Kosmidis

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

unless_null <- function(x, if_null) {
    if (is.null(x)) {
        if_null
    }
    else {
        x
    }
}

## named_list is the function namedList from utils.R of the ordinal R
## package version 2019.4-25 by Rune Haubo Bojesen Christensen
named_list <- function(...) {
    setNames(list(...), nm=sapply(as.list(match.call()), deparse)[-1])
}
