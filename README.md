brglm2
======

[**brglm2**](https://github.com/ikosmidis/brglm2) provides tools for the estimation and inference from generalized linear models using various methods for bias reduction ([Kosmidis, 2014](https://doi.org/10.1002/wics.1296)). Reduction of estimation bias is achieved either through the adjusted score equations approach in [Firth (1993)](https://doi.org/10.1093/biomet/80.1.27) and [Kosmidis & Firth (2009)](https://doi.org/10.1093/biomet/asp055) or through the direct subtraction of an estimate of the bias of the maximum likelihood estimator from the maximum likelihood estimates as prescribed in [Cordeiro and McCullagh (1991)](http://www.jstor.org/stable/2345592)

In the special case of generalized linear models for binomial and multinomial responses, the adjusted score equations approach returns estimates with improved frequentist properties, that are also always finite, even in cases where the maximum likelihood estimates are infinite (e.g. complete and quasi-complete separation).

### Installation

Install the development version from github:

``` r
# install.packages("devtools")
devtools::install_github("ikosmidis/brglm2")
```

### Quasi Fisher scoring for solving adjusted score equations

The workhorse function in **brglm2** is [`brglmFit`](https://github.com/ikosmidis/brglm2/blob/master/R/brglmFit.R), which can be passed directly to the `method` argument of the `glm` function and . `brglmFit` implements a quasi
[Fisher scoring](https://en.wikipedia.org/wiki/Scoring_algorithm)
procedure, whose special cases result in various explicit and implicit
bias reduction methods for generalized linear models.

The [iteration vignette](https://github.com/ikosmidis/brglm2/blob/master/vignettes/iteration.pdf) describes the iteration and gives the mathematical details for the bias-reducing adjustments to the score functions for generalized linear models.

The classification of bias reduction methods into explicit and implicit is
as given in [Kosmidis (2014)](https://doi.org/10.1002/wics.1296).

### References and resources

**brglm2** was presented by [Ioannis Kosmidis](https://www.ucl.ac.uk/~ucakiko/) at the [useR! 2016 international R User conference](http://user2016.org) at University of Stanford on 16 June 2016. The presentation was titled "Reduced-bias inference in generalized linear models" and can be watched online at this [link](https://channel9.msdn.com/Events/useR-international-R-User-conference/useR2016/brglm-Reduced-bias-inference-in-generalized-linear-models).
