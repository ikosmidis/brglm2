brglm2
======

Estimation and inference from generalized linear models based on implicit methods for bias reduction (see Kosmidis, 2014, WIRE Computational Statistics). brglm2 can achieve reduction of estimation bias either through the adjusted score equations approach in Firth (1993, Biometrika) and Kosmidis and Firth (2009, Biometrika), or through the direct subtraction of an estimate of the bias of the maximum likelihood estimator from the maximum likelihood estimates. In the special case of generalized linear models for binomial and multinomial responses, the adjusted score equations approach returns estimates with improved frequentist properties, that are also always finite, even in cases where the maximum likelihood estimates are infinite (e.g. complete and quasi-complete separation). Estimation in all cases takes place via a modified Fisher scoring algorithm.

### Installation

Install the development version from github:

``` r
# install.packages("devtools")
devtools::install_github("ikosmidis/brglm2")
```
