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



