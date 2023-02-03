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
#'     \url{https://users.stat.ufl.edu/~aa/glm/data/}. The latter link
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
#' This data set contains survival data on government coalitions in
#' parliamentary democracies (Belgium, Canada, Denmark, Finland,
#' France, Iceland, Ireland, Israel, Italy, Netherlands, Norway,
#' Portugal, Spain, Sweden, and the United Kingdom) for the period
#' 1945-1987.  For parsimony, country indicator variables are omitted
#' in the sample data.
#'
#' @format
#'
#' A data frame with 314 rows and the 7 variables "duration",
#' "ciep12", "invest", "fract", "polar", "numst2", and "crisis".  For
#' variable descriptions, please refer to King et al (1990).
#'
#' @note
#'
#' Data is as it was provided in the
#' [\pkg{Zeilig}](https://cran.r-project.org/package=Zelig) R package.
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

#' The effects of AZT in slowing the development of AIDS symptoms
#'
#' The data is from a 3-year study on the effects of AZT in slowing the
#' development of AIDS symptoms. 338 veterans whose immune systems
#' were beginning to falter after infection with the AIDS virus were
#' randomly assigned either to receive AZT immediately or to wait
#' until their T cells showed severe immune weakness.
#'
#' @format A data frame with 4 rows and 4 variables:
#'
#' * `symptomatic`: counts of veterans showing AIDS symptoms during the 3-year study
#'
#' * `asymptomatic`: counts of veterans not showing AIDS symptoms during the 3-year study
#'
#' * `race`: the race of the veterans with levels `"White"` and `"Black"`
#'
#' * `AZT_use`: whether the veterans received AZT immediately (`"Yes"`)
#' or waited until their T cells showed severe immune weakness (`"No"`)
#'
#' @source
#'
#' The data set is analyzed in Agresti (2002, Subsection 5.4.2). Its
#' original source is New York Times, Feb. 15, 1991.
#'
#' @seealso
#'
#' [brmultinom()]
#'
#' @references
#'
#' Agresti A (2002). *Categorical Data Analysis*.  Wiley Series in
#' Probability and Statistics. Wiley.
"aids"

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
