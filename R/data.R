#' Habitat preferences of lizards
#'
#' @format A data frame with 23 rows and 6 columns:
#'
#' * `grahami`. count of grahami lizards
#' * `opalinus`. count of opalinus lizards
#' * `height`. a factor with levels `<5ft`, `>=5ft`
#' * `diameter`. a factor with levels `<=2in`, `>2in`
#' * `light`. a factor with levels `sunny`, `shady`
#' * `time`. a factor with levels `early`, `midday`, `late`
#'
#' The variables `grahami` and `opalinus` are counts of two lizard
#' species at two different perch heights, two different perch
#' diameters, in sun and in shade, at three times of day.
#'
#' @seealso
#'
#' [brglm_fit()]

#'
#' @source
#'
#'   McCullagh P, Nelder J A (1989) _Generalized Linear
#'   Models_ (2nd Edition).  London: Chapman and Hall.
#'
#' Originally from
#'
#'     Schoener T W (1970) Nonsynchronous spatial overlap of lizards
#'     in patchy habitats.  _Ecology_ *51*, 408-418.
#'
"lizards"

#' Histology grade and risk factors for 79 cases of endometrial cancer
#'
#' @format A data frame with 79 rows and 4 variables:
#'
#' * `NV`: neovasculization with coding 0 for absent and 1 for present
#' * `PI`: pulsality index of arteria uterina
#' * `EH`: endometrium height
#' * `HG` histology grade with coding 0 for low grade and 1 for high grade
#'
#' @source The packaged data set was downloaded in `.dat` format from
#'     \url{http://www.stat.ufl.edu/~aa/glm/data}. The latter link
#'     provides the data sets used in Agresti (2015).
#'
#'     The endometrial data set was first analyzed in Heinze and
#'     Schemper (2002), and was originally provided by Dr
#'     E. Asseryanis from the Medical University of Vienna.
#'
#' @seealso
#'
#' [brglm_fit()]
#'
#'
#' @references
#'
#' Agresti A (2015). *Foundations of Linear and Generalized Linear
#' Models*.  Wiley Series in Probability and Statistics. Wiley.
#'
#' Heinze G, Schemper M (2002). A Solution to the Problem of
#' Separation in Logistic Regression. *Statistics in Medicine*,
#' **21**, 2409–2419. \doi{10.1002/sim.1047}.
#'
#' Kosmidis I, Firth D (2021). Jeffreys-prior penalty, finiteness
#' and shrinkage in binomial-response generalized linear
#' models. *Biometrika*, **108**, 71-82. \doi{10.1093/biomet/asaa052}.
#'
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
#' [brglm_fit()]
#'
#' @references
#'
#'  King G, Alt J E, Burns N E, Laver M. (1990). A Unified
#'  Model of Cabinet Dissolution in Parliamentary
#'  Democracies. *American Journal of Political Science*, **34**,
#'  846-870. \doi{10.2307/2111401}.
#'
#'  King G, Alt J E, Burns N E, Laver M. ICPSR
#'  Publication Related Archive, 1115.
#'
"coalition"


#' Alligator food choice data
#'
#' @format A data frame with 80 rows and 5 variables:
#'
#' * `foodchoice`: primary food type, in volume, found in an alligator’s stomach, with levels `fish`, `invertebrate`,`reptile`, `bird`, `other`
#' * `lake`: lake of capture with levels `Hancock`, `Oklawaha`, `Trafford`, `George`.
#' * `gender`: gender of the alligator with levels `Male` and `Female`
#' * `size`: size of the alligator with levels `<=2.3` meters long and `>2.3` meters long
#' * `freq`: number of alligators for each foodchoice, lake, gender and size combination
#'
#'
#' @source
#'
#' The alligators data set is analyzed in Agresti (2002, Subsection 7.1.2).
#'
#' @seealso
#'
#' [brmultinom()]
#'
#' @references
#'
#' Agresti A (2002). *Categorical Data Analysis*.  Wiley Series in
#' Probability and Statistics. Wiley.
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
#'
#' * `research`: opinion about funding stem cell research with levels `definitely`, `probably`, `probably not`, `definitely not`
#' * `gender`: the gender of the respondent with levels `female` and `male`
#' * `religion`: the fundamentalism/liberalism of one’s religious beliefs with levels `fundamentalist`, `moderate`,
#' `liberal`
#' `frequency`: the number of times a respondent fell in each of the combinations of levels for `research`, `religion` and `gender`
#'
#'
#' @seealso
#'
#' [bracl()]
#'
#' @source
#'
#' The `stemcell` data set is analyzed in Agresti (2010, Subsection 4.1.5).
#'
#' @references
#'
#' Agresti A (2010). *Analysis of Ordinal Categorical Data* (2nd edition).  Wiley Series in
#' Probability and Statistics. Wiley.
#'
"stemcell"
