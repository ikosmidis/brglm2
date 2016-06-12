context("test correct cacluation of null deviace")

## For the null deviance:
##
##
## If there is an intercept but not an offset then the fitted
## value is the weighted average and is calculated easily below
##
## If there is an offset but not an intercept then the fitted
## value is the inverse link evaluated at the offset (APPLIES)
##
## If there is neither an offset nor an intercept then the fitted
## values is the inverse link at zero (and hence covered by
## linkinv(zero) (APPLIES)
## If there is an intercept and an offset then, for calculating
## the null deviance glm will make a call to the fitter to fit the
## glm with intercept and the offset Make sure that glm.fit is
## being used here

data(coalition, package = "Zelig")
## The maximum likelihood fit with log link
coalitionMLoff <- glm(duration ~ fract + offset(numst2), family = Gamma("log"), data = coalition)
## The bias-reduced fit
coalitionBRoff <- update(coalitionMLoff, method = "brglmFit")
## The bias-corrected fit
coalitionBCoff <- update(coalitionMLoff, method = "brglmFit", type = "correction")

