# Copyright (C) 2016-2021 Ioannis Kosmidis

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



get_type_description <- function(type, parenthesized = TRUE) {
    pp <- function(txt) {
        ifelse(parenthesized, paste0("(", txt, ")"), txt)
    }
    switch(type,
           "ML" = pp("maximum likelihood"),
           "correction" = pp("bias correction"),
           "AS_mean" = pp("mean bias-reducing adjusted score equations"),
           "AS_median" = pp("median bias-reducing adjusted score equations"),
           "AS_mixed" = pp("mixed bias-reducing adjusted score equations"),
           "MPL_Jeffreys" = pp("maximum penalized likelihood with Jeffreys'-prior penalty")
           )
}
