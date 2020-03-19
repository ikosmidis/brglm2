# brglm 0.6.1
## Bug fixes
* Fixed bug in AIC reported by `print.summary` for `brmultinom` and `bracl`
* `detect_separation` now handles one-column model matrices correctly

## Other improvements, updates and additions
* Documentation improvements and typo fixes

# brglm 0.6
## New functionality
* `brglmFit` can now do maximum penalized likelihood with powers of the Jeffreys prior as penalty (`type = "MPL_Jeffreys`) for all supported generalized linear models. See `brglmControl` and `brglmFit` for details.

## Other improvements, updates and additions
* Documentation updates and improvements
* Updated vignettes to include maximum penalized likelihood with powers of the Jeffreys prior as penalty
* New examples in `?brglmFit`

# brglm 0.5.2
## Bug fixes
* `print.brmultinom` is now exported, so `bracl` and `brmultinom` fits print correctly

## New functionality
* Added `response_adjustment` argument in `brglmControl` to allow for more fine-tuning of the starting values when `brglmFit` is called with `start = NULL`

## Other improvements, updates and additions
* Documentation updates and improvements
* Added Kosmidis et al (2019) in the description file
* Added tests for `brglmControl`

# brglm 0.5.1

## Other improvements, updates and additions
* Fixed typos in vignettes and documentation
* Added ORCHID for Ioannis Kosmidis in DESCRIPTION

# brglm 0.5.0
## Bug fixes
* `brglmFit` now works as expected with custom link functions (mean and median bias reduction)
* `brglmFit` respects the specification of the transformation argument in `brglmControl`
* Fixed bug in the computation of the QR decomposition under aliasing in `brglmFit`
* Other minor bug fixes and performance improvements
* Protection against use of `quasi`, `quasibinomial` and `quasibinomial` families and documentation update.

## New functionality
* Added `bracl` for fitting adjacent category logit models for ordinal responses using maximum likelihood, mean bias reduction, and median bias reduction and associated methods (`logLik`, `summary` and so on)
* Added `predict` methods for `brmultinom` and `bracl`
@ Added `residuals` methods for `brmultinom` and `bracl` (residuals of the equivalent Poisson log-linear model)
* Added the `mis` link functions for accounting for misclassification in binomial response models (Neuhaus, 1999, Biometrika)

## Other improvements, updates and additions
* Improved `summary` method for `brmultinom` objects
* Better starting values for null fits
* Added references to [arxiv:1804.04085](https://arxiv.org/abs/1804.04085) in documentation
* Updated reference to [Kenne Pagui et al (2017)](https://doi.org/10.1093/biomet/asx046)

# brglm 0.1.8
## Other improvements, updates and additions
* Improved documentation examples
* Removed warning about observations with non-positive weights from brmultinom
* Updated email address for Ioannis Kosmidis in brglmFit

## Bug fixes
* brmultinom returns a fitted values matrix that respects the dimension of data
* Fixed bug on condition for NA dispersion for models with 0 df resid

# brglm 0.1.7

## Other improvements, updates and additions
* Eliminated errors from markdown chunks in multinomial vignette

# brglm2 0.1.6

## Bug fixes
* Compatibility with new version of enrichwith

## Other improvements, updates and additions
* New email for Ioannis Kosmidis

# brglm2 0.1.5

## Bug fixes

## New functionality
* Added `type = AS_mixed` as an option to use **mean-bias reducing score functions** for the regression parameters and **median-bias reducing score functions** for the dispersion in models with unknown dispersion
* `check_infinite_estimates` now accepts `brmultinom` objects
* Added `singular.ok` argument to `brglmFit` and `detect_separation` methods in line with the update of `glm.fit`

## Other improvements, updates and additions
* less strict tolerance in `brglm_control`
* Updates to help files
* Fixed typos in iteration vignette
* Added URL and bugreports in Description
* Added new tests

# brglm2 0.1.4

## Bug fixes
* `brglmControl` is now exported
* `slowit` did nothing; now included in iteration

## New functionality
* The `detect_separation` `method` for the `glm` function can be used to check for separation in binomial response settings without fitting the model. This relies on a port of Kjell Konis' `safeBinaryRegression:::separator` function (see ?detect_separation)
* brglm2 provides estimation via **median-bias reducing score functions** with `type = "AS_median"`
* brglm2 provides camel and underscored aliases for basic methods (`brglmFit`, `brglm_fit`, `detectSeparation`, `detect_separation`, `brglm_control`, `brglmControl`, `detectSeparationControl`, `detect_separation_control`, `checkInfiniteEstimates`, `check_infinite_estimates`)

## Other improvements, updates and additions
* Minor enhancements in the codebase
* The inverse expected information matrix is computed internally using cho2inv
* Internal changes to have more meaningful variable names
* Renamed detect_infinite* to check_infinite

# brglm2 0.1.3

## Bug fixes

## New functionality

## Other improvements, updates and additions
* Fixed typo in f_{Y_i}(y) in iteration vignette (thanks to Eugene
  Clovis Kenne Pagui for spotting)

# brglm2 0.1.2

* First release



