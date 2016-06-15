##  Some examples of the use of "brpr" to fit a multinomial logit
##  model by penalized (reduced-bias) maximum likelihood.
##
##  Examples prepared by I. Kosmidis, University College London

source(system.file("inst", "brpr/brpr.R", package = "brglm2"))

library(pmlr)
library(MASS) # for data(housing)
library(brglm) # for the penalizedDeviance profileModel objective
library(gnm) # for pickCoef and expandCategorical used in coefinds and expandMult

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


#####################################################################
## Analysis of the Coppenhagen housing data set in the MASS library
#####################################################################
data(housing)
contrasts(housing$Sat) <- contr.treatment(3, base = 1)

housbrpr <- brpr(Freq ~ Infl * Type * Cont +
                  Sat * (Infl + Type + Cont), data = housing,
                  fixed.totals = Infl:Type:Cont)
houspmlr <- pmlr(Sat ~ Infl + Type + Cont, weights = Freq,
                 data = housing, method = "wald")

## The estimates are the same at least up to 5 decimals
coef(housbrpr)[coefinds(housbrpr)]
houspmlr$coefficients


#####################################################################
## Analysis of the hepatitis data set in Bull et. al. (2007)
#####################################################################
data(hepatitis)
## Construct a variable with the multinomial categories according to
## the HCV and nonABC columns
hepat <- hepatitis
hepat$type <- with(hepat, factor(1 - HCV*nonABC + HCV + 2 * nonABC))
hepat$type <- factor(hepat$type, labels = c("noDisease", "C", "nonABC"))
contrasts(hepat$type) <- contr.treatment(3, base = 1)
## compare the result with the one that pmlr gives
hepatnew <- expandMult(data = hepat, catvar = "type", countsvar = "counts")
heppmlr <- pmlr(type ~ group + time + group:time,
                data = hepatnew, weights = counts, method = "wald",
                penalized = TRUE)
hepbrpr <- brpr(counts ~ group*time + type*(group*time),
                 fixed.totals = group:time, data = hepatnew,
                 epsilon = 1e-14) ## very strict epsilon

## The estimates are the same at least up to 5 decimals
coef(hepbrpr)[coefinds(hepbrpr)]
heppmlr$coefficients

## Speed comparison (it takes a while to run the next 2 lines!)
## Shows that brpr with epsilon 1e-14 is roughly 5 times faster
## than pmlr for this data set.
system.time(for (i in 1:100) {
  pmlr(type ~ group + time + group:time,
       data = hepatnew, weights = counts, method = "wald",
       penalized = TRUE) })
system.time(for (i in 1:100) {
  brpr(counts ~ group*time + type*(group*time),
        fixed.totals = group:time, data = hepatnew,
        epsilon = 1e-14) })

## Profile confidence intervals based on the penalized likelihood
hepbrprCIs <- confintModel(hepbrpr, objective = "penalizedDeviance",
                           quantile = qchisq(0.95, 1),
                           X = model.matrix(hepbrpr),
                           which = coefinds(hepbrpr),
                           method = "zoom",
                           endpoint.tolerance = 1e-04)
heppmlr <- pmlr(type ~ group + time + group:time,
                data = hepatnew, weights = counts,
                method = "likelihood", penalized = TRUE)
nn <- dimnames(heppmlr$coefficients)
coefnames <- c(paste(nn[[2]], nn[[3]][1], sep = ":"),
               paste(nn[[2]], nn[[3]][2], sep = ":"))
heppmlrCIs <- cbind(heppmlr$CI.lower, heppmlr$CI.upper)
rownames(heppmlrCIs) <- coefnames

# The confidence intervals from pmlr and profileModel seem to agree
# only with slight differences. But as far as those differences are
# concerned, I trust confintModel's result; those confidence intervals
# are obtained through binary search.
heppmlrCIs
hepbrprCIs


#####################################################################
## Analysis of the enzymes data set in ?pmlr
#####################################################################
data(enzymes)
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
enzpmlr <- pmlr(Group ~ AST + GLDH, weights = counts,
                data = enz, method = "wald")
enzbrpr <- brpr(counts ~ -1 + ind + Group * (AST + GLDH),
                 fixed.totals = ind, data = enz,
                 epsilon = 1e-14)
## brpr appears slower than pmlr in this case, most probably because
## of the large number of nuisances.

## The estimates are the same at least up to 5 decimals
coef(enzbrpr)[coefinds(enzbrpr)]
enzpmlr$coefficients
