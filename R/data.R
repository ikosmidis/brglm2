#' Habitat Preferences of Lizards
#'
#' The lizards data frame has 23 rows and 6 columns. Variables
#' \code{grahami} and \code{opalinus} are counts of two lizard species
#' at two different perch heights, two different perch diameters, in
#' sun and in shade, at three times of day.
#'
#' \itemize{
#'   \item grahami. count of grahami lizards
#'   \item opalinus. count of opalinus lizards
#'   \item height. a factor with levels \code{<5ft}, \code{>=5ft}
#'   \item diameter. a factor with levels \code{<=2in}, \code{>2in}
#'   \item light. a factor with levels \code{sunny}, \code{shady}
#'   \item time. a factor with levels \code{early}, \code{midday}, \code{late}
#' }
#'
#'
#' @source
#'
#'   McCullagh, P. and Nelder, J. A. (1989) _Generalized Linear
#'   Models_ (2nd Edition).  London: Chapman and Hall.
#'
#' Originally from
#'
#'     Schoener, T. W. (1970) Nonsynchronous spatial overlap of lizards
#'     in patchy habitats.  _Ecology_ *51*, 408-418.
#'
"lizards"

#' Histology grade and risk factors for 79 cases of endometrial cancer
#'
#' @format A data frame with 79 rows and 4 variables:
#' \describe{
#'
#' \item{NV}{neovasculation with coding 0 for absent and 1 for present}
#'
#' \item{PI}{pulsality index of arteria uterina}
#'
#' \item{EH}{endometrium heigh}
#'
#' \item{HG}{histology grade with coding 0 for low grade and 1 for high grade}
#'
#' }
#'
#' @source The packaged data set was downloaded in \code{.dat} format
#'     from \url{http://www.stat.ufl.edu/~aa/glm/data}. The latter
#'     link provides the data sets used in Agresti (2015).
#'
#'     The endometrial data set was first analysed in Heinze and
#'     Schemper (2002), and was originally provided by Dr
#'     E. Asseryanis from the Medical University of Vienna.
#'
#' @references
#'
#' Agresti, A. (2015). *Foundations of Linear and Generalized Linear
#' Models*.  Wiley Series in Probability and Statistics. Wiley
#'
#' Heinze, G., & Schemper, M. (2002). A Solution to the Problem of
#' Separation in Logistic Regression. *Statistics in Medicine*, **21**, 2409–2419
#'
"endometrial"


#' Coalition data
#'
#' @note
#'
#' Data is as provided in the Zeilig R package (\url{https://cran.r-project.org/package=Zelig})
#'
#' @references
#'
#'  King, G., Alt, J. E., Burns, N. E. and Laver, M. (1990). A Unified
#'  Model of Cabinet Dissolution in Parliamentary
#'  Democracies. *American Journal of Political Science*, **34**,
#'  846-870.
#'
#'  King, G., Alt, J. E., Burns, N. E. and Laver, M. ICPSR
#'  Publication Related Archive, 1115.
#'
"coalition"
