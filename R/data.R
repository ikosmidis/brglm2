#' Habitat preferences of lizards
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
#' @seealso
#'
#' \code{\link{brglm_fit}}

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
#' \item{NV}{neovasculization with coding 0 for absent and 1 for present}
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
#' @seealso
#'
#' \code{\link{brglm_fit}}
#' 
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
#' @seealso
#'
#' \code{\link{brglm_fit}}
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


#' Alligator food choice data
#'
#' @format A data frame with 80 rows and 5 variables:
#'
#' \describe{
#'
#' \item{foodchoice}{primary food type, in volume, found in an alligator’s stomach, with levels \code{fish}, \code{invertebrate}, \code{reptile}, \code{bird}, \code{other}}
#'
#' \item{lake}{lake of capture with levels \code{Hancock}, \code{Oklawaha}, \code{Trafford}, \code{George}}
#'
#' \item{gender}{gender of the alligator with levels \code{Male} and \code{Female}}
#'
#' \item{size}{size of the alligator with levels \code{<=2.3} meters long and \code{>2.3} meters long}
#'
#' \item{freq}{number of alligators for each foodchoice, lake, gender and size combination}
#'
#' }
#'
#' @source
#'
#' The alligators data set is analysed in Agresti (2002, Subsection~7.1.2).
#'
#' @seealso
#'
#' \code{\link{brmultinom}}
#' 
#' @references
#'
#' Agresti, A. (2002). *Categorical Data Analysis*.  Wiley Series in
#' Probability and Statistics. Wiley
#'
"alligators"



#' Opinion on Stem Cell Research and Religious Fundamentalism
#'
#' A data set from the 2006 General Social Survey that shows the
#' relationship in the United States between opinion about funding
#' stem cell research and the fundamentalism/liberalism of one’s
#' religious beliefs, stratified by gender.
#'
#' @format A data frame with 24 rows and 4 variables:
#' \describe{
#'
#' \item{research}{opinion about funding stem cell research with levels \code{definitely}, \code{probably}, \code{probably not}, \code{definitely not}}
#'
#' \item{gender}{the gender of the respondent with levels \code{female} and \code{male}}
#'
#' \item{religion}{the fundamentalism/liberalism of one’s religious
#' beliefs with levels \code{fundamentalist}, \code{moderate},
#' \code{liberal}}
#'
#' \item{frequency}{the number of times a respondent fell in each of the combinations of levels for \code{research}, \code{religion} and \code{gender}}
#'
#' }
#'
#' @seealso
#'
#' \code{\link{bracl}}
#' 
#' @source
#'
#' The \code{stemcell} data set is analysed in Agresti (2010, Subsection~4.1.5).
#'
#' @references
#'
#' Agresti, A. (2010). *Analysis of Ordinal Categorical Data* (2nd edition).  Wiley Series in
#' Probability and Statistics. Wiley
#'
"stemcell"
