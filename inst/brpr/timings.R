##  Some examples of the use of "brpr" to fit a multinomial logit
##  model by penalized (reduced-bias) maximum likelihood.
##
##  Examples prepared by I. Kosmidis, University College London

source(system.file("inst", "brpr/brpr.R", package = "brglm2"))

library("pmlr")
library("MASS") # for data(housing)
library("gnm") # for pickCoef and expandCategorical used in coefinds and expandMult

## coefinds takes as input the glm object that brpr results and
## returns the location of the coefficients for the fitted multinomial
## model.
coefinds <- function(obj) {
  require(gnm)
  cats <- as.character(obj$call$formula[3])
  cats <- strsplit(gsub(" ", "", cats), "\\*\\(")[[1]][1]
  cats <- rev(strsplit(cats, "\\+")[[1]])[1]
  cc <- pickCoef(obj, cats)
  cc
}

## expandMult will epxand the multinomial data set so that all
## available categories are represented in each covariate class. If
## there are non-available categories in the original data set they
## are included in the resulted data frame and are assigned frequency
## zero.
expandMult <- function(data, catvar, countsvar) {
  require(gnm)
  temp <- expandCategorical(data = data, catvar = catvar,
                            sep = ".", countvar = "CcCounts",
                            idvar = "id", as.ordered = FALSE,
                            group = TRUE)
  ## Get the correct category counts
  temp[[countsvar]] <- temp[[countsvar]] * temp$CcCounts
  temp$CcCounts <-   temp$id <- NULL
  temp
}

timing <- function(..., R = 10) {
    system.time(replicate(R, ...))
}

#####################################################################
## Analysis of the Coppenhagen housing data set in the MASS library
#####################################################################
data("housing", package = "MASS")
contrasts(housing$Sat) <- contr.treatment(3, base = 1)

timing(housbrpr <- brpr(Freq ~ Infl * Type * Cont +
                            Sat * (Infl + Type + Cont), data = housing,
                        fixed.totals = Infl:Type:Cont, epsilon = 1e-06), R = 50)
timing(houspmlr <- pmlr(Sat ~ Infl + Type + Cont, weights = Freq,
                        data = housing, method = "wald"), R = 50)
timing(housbrmultinom <- brmultinom(Sat ~ Infl + Type + Cont, weights = Freq,
                                    data = housing, epsilon = 1e-06), R = 50)

#####################################################################
## Analysis of the hepatitis data set in Bull et. al. (2007)
#####################################################################
data("hepatitis", package = "pmlr")
## Construct a variable with the multinomial categories according to
## the HCV and nonABC columns
hepat <- hepatitis
hepat$type <- with(hepat, factor(1 - HCV*nonABC + HCV + 2 * nonABC))
hepat$type <- factor(hepat$type, labels = c("noDisease", "C", "nonABC"))
contrasts(hepat$type) <- contr.treatment(3, base = 1)
## compare the result with the one that pmlr gives
hepatnew <- expandMult(data = hepat, catvar = "type", countsvar = "counts")

timing(heppmlr <- pmlr(type ~ group + time + group:time,
                       data = hepatnew, weights = counts, method = "wald",
                       penalized = TRUE), R = 50)
timing(hepbrmultinom <- brmultinom(type ~ group + time + group:time,
                                   data = hepatnew, weights = counts,
                                   epsilon = 1e-06), R = 50)
timing(hepbrpr <- brpr(counts ~ group*time + type*(group*time),
                       fixed.totals = group:time, data = hepatnew,
                       epsilon = 1e-06), R = 50)

#####################################################################
## Analysis of the enzymes data set in ?pmlr
#####################################################################
data("enzymes", package = "pmlr")
## Exclude patients in Group 4 (post-necrotic cirrhosis)
enzymes <- enzymes[enzymes$Group != 4,]
## Center and scale covariates
AST <- scale(log(enzymes$AST))
ALT <- scale(log(enzymes$ALT))
GLDH <- scale(log(enzymes$GLDH))
OCT <- scale(log(enzymes$OCT))
enzymes <- data.frame(Patient = enzymes$Patient,
                      Group = enzymes$Group, AST, ALT, GLDH, OCT)
## Remove 10 observations to create separation
enzymes <- enzymes[-c(9, 18, 33, 58, 61, 77, 94, 97, 99, 100),]
## Multinomial: acute viral hepatitis and aggressive chronic hepatits
## vs. persistent chronic hepatitis
## Assign Group 2 (persistent chronic hepatitis) as baseline category
enzymes$Group <- factor(enzymes$Group, levels=c("2","1","3"))
## Re-express data set in an appropriate form for brpr
enzymes$counts <- rep(1, nrow(enzymes))
enzymes$ind <- factor(1:131)
enz <- expandMult(enzymes, "Group", "counts")
timing(enzpmlr <- pmlr(Group ~ AST + GLDH, weights = counts,
                       data = enz, method = "wald"), R = 10)
timing(enzbrpr <- brpr(counts ~ -1 + ind + Group * (AST + GLDH),
                       fixed.totals = ind, data = enz,
                       epsilon = 1e-14), R = 10)
timing(enzbrmultinom <- brmultinom(Group ~ AST + GLDH, weights = counts,
                                   data = enz,
                                   epsilon = 1e-14), R = 10)
## max(abs(matrix(enzbrpr$fitted.values, ncol = 3, byrow = TRUE) - enzbrmultinom$fitted.values))


timing <- function(..., R = 10) {
    system.time(replicate(R, ...))
}
