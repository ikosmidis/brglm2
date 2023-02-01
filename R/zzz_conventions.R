# Copyright (C) 2017- Ioannis Kosmidis

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

## Some convenience names

#' @rdname brglmControl
#' @export
brglm_control <- brglmControl

#' @rdname brglmFit
#' @export
brglm_fit <- brglmFit

#### Method conventions

#' @rdname residuals.brmultinom
#' @method residuals bracl
#' @export
residuals.bracl <- residuals.brmultinom

