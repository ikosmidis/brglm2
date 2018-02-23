# brglm 0.1.7

## Other improvements, updates and addition
* Eliminated errors from markdown chunks in multinomial vignette

# brglm2 0.1.6

## Bug fixes
* Compatibility with new version of enrichwith

## New functionality

## Other improvements, updates and addition
* New email for Ioannis Kosmidis


# brglm2 0.1.5

## Bug fixes

## New functionality
* Added `type = AS_mixed` as an option to use **mean-bias reducing score functions** for the regression parameters and **median-bias reducing score functions** for the dispersion in models with uknown dispersion
* `check_infinite_estimates` now accepts `brmultinom` objects
* Added `singular.ok` argument to `brglmFit` and `detect_separation` methods in line with the update of `glm.fit`

## Other improvements, updates and addition
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



