[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/brglm2)](https://cran.r-project.org/package=brglm2)
[![Travis-CI Build Status](https://travis-ci.org/ikosmidis/brglm2.svg?branch=master)](https://travis-ci.org/ikosmidis/brglm2)
[![Coverage Status](https://img.shields.io/codecov/c/github/ikosmidis/brglm2/master.svg)](https://codecov.io/github/ikosmidis/brglm2?branch=master)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

brglm2
======

[**brglm2**](https://github.com/ikosmidis/brglm2) provides tools for the estimation and inference from generalized linear models using various methods for bias reduction ([Kosmidis, 2014](https://doi.org/10.1002/wics.1296)). **brglm2** supports all generalized linear models supported in R, and provides methods for multinomial logistic regression (nominal responses) and adjacent category models (ordinal responses). 

Reduction of estimation bias is achieved by solving either the mean-bias reducing adjusted score equations in [Firth (1993)](https://doi.org/10.1093/biomet/80.1.27) and [Kosmidis & Firth (2009)](https://doi.org/10.1093/biomet/asp055) or the median-bias reducing adjusted score equations in [Kenne et al (2016)](https://arxiv.org/abs/1604.04768), or through the direct subtraction of an estimate of the bias of the maximum likelihood estimator from the maximum likelihood estimates as prescribed in [Cordeiro and McCullagh (1991)](http://www.jstor.org/stable/2345592). [Kosmidis et al (2019)](https://doi.org/10.1007/s11222-019-09860-6) provides a unifying framework and algorithms for mean and median bias reduction for the estimation of generalized linear models. 

In the special case of generalized linear models for binomial and multinomial responses (both ordinal and nomial), the adjusted score equations return estimates with improved frequentist properties, that are also always finite, even in cases where the maximum likelihood estimates are infinite (e.g. complete and quasi-complete separation). See, [Kosmidis & Firth (2019)](http://arxiv.org/abs/1812.01938) for the proof of the latter result in the case of mean bias reduction for logistic regression (and, for more general binomial-response models where the likelihood is penalized by a power of the Jeffreys invariant prior).

**brglm2** also provides *pre-fit* and *post-fit* methods for the detection of separation and of infinite maximum likelihood estimates in binomial response generalized linear models (see `?detect_separation` and `?check_infinite_estimates`).

### Installation

Install the development version from github:

``` r
# install.packages("devtools")
devtools::install_github("ikosmidis/brglm2")
```

### Solving adjusted score equations using quasi-Fisher scoring

The workhorse function in **brglm2** is
[`brglmFit`](https://github.com/ikosmidis/brglm2/blob/master/R/brglmFit.R),
which can be passed directly to the `method` argument of the `glm`
function. `brglmFit` implements a quasi [Fisher
scoring](https://en.wikipedia.org/wiki/Scoring_algorithm) procedure,
whose special cases result in a range of explicit and implicit bias
reduction methods for generalized linear models. Bias reduction for
multinomial logistic regression (nominal responses) can be performed
using the function `brmultinom`, and for adjacent category models
(ordinal responses) using the function `bracl`. Both `brmultinom` and
`bracl` rely on `brglmFit`.

The [iteration vignette](https://cran.r-project.org/package=brglm2/vignettes/iteration.html) and [Kosmidis et al
(2019)](https://doi.org/10.1007/s11222-019-09860-6) apresent the iteration and give mathematical details for the bias-reducing adjustments to the score functions for generalized linear models.

The classification of bias reduction methods into explicit and implicit is as given in [Kosmidis (2014)](https://doi.org/10.1002/wics.1296).

### References and resources

**brglm2** was presented by [Ioannis Kosmidis](http://www.ikosmidis.com) at the useR! 2016 international conference at University of Stanford on 16 June 2016. The presentation was titled "Reduced-bias inference in generalized linear models" and can be watched online at this [link](https://channel9.msdn.com/Events/useR-international-R-User-conference/useR2016/brglm-Reduced-bias-inference-in-generalized-linear-models).


Motivation, details and discussion on the methods that **brglm2** implements are provided in

Kosmidis, I, Kenne Pagui, E C, Sartori N. (2017). Mean and median bias reduction in generalized linear models. To appear in [*Statistics and Computing*](https://doi.org/10.1007/s11222-019-09860-6). *arXiv*, [arXiv:1710.11217](http://arxiv.org/abs/1804.04085). 


