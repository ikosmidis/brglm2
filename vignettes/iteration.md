The **brglm2** package
======================

[**brglm2**](https://github.com/ikosmidis/brglm2) provides tools for the
estimation and inference from [generalized linear
models](https://en.wikipedia.org/wiki/Generalized_linear_model) using
various methods for bias reduction (Kosmidis 2014). Reduction of
estimation bias is achieved either through the adjusted score equations
approach in Firth (1993) and Kosmidis and Firth (2009), or through the
direct subtraction of an estimate of the bias of the maximum likelihood
estimator from the maximum likelihood estimates as prescribed in
Cordeiro and McCullagh (1991).

In the special case of generalized linear models for binomial and
multinomial responses, the adjusted score equations approach returns
estimates with improved frequentist properties, that are also always
finite, even in cases where the maximum likelihood estimates are
infinite, like in complete and quasi-complete separation as defined in
Albert and Anderson (1984).

The workhorse function is
[`brglmFit`](https://github.com/ikosmidis/brglm2/blob/master/R/brglmFit.R),
which can be passed directly to the `method` argument of the `glm`
function. `brglmFit` implements a quasi [Fisher
scoring](https://en.wikipedia.org/wiki/Scoring_algorithm) procedure,
whose special cases result in various explicit and implicit bias
reduction methods for generalized linear models (the classification of
bias reduction methods into explicit and implicit is given in Kosmidis
2014).

This vignette
=============

This vignette

-   presents the bias-reducing adjustments to the score functions for
    generalized linear models
-   describes the fitting algorithm at the core of **brglm2**

Other resources
===============

The bias-reducing quasi Fisher scoring iteration is also described in
detail in the [bias
vignette](https://cran.r-project.org/package=enrichwith/vignettes/bias.html)
of the [**enrichwith**](https://cran.r-project.org/package=enrichwith) R
package. Kosmidis and Firth (2010) describe a parallel quasi
[Newton-Raphson](https://en.wikipedia.org/wiki/Newton%27s_method)
procedure.

Most of the material in this vignette comes from a presentation by
[Ioannis Kosmidis](https://www.ucl.ac.uk/~ucakiko/) at the [useR! 2016
international R User conference](http://user2016.org) at University of
Stanford on 16 June 2016. The presentation was titled "Reduced-bias
inference in generalized linear models" and can be watched online at
this
[link](https://channel9.msdn.com/Events/useR-international-R-User-conference/useR2016/brglm-Reduced-bias-inference-in-generalized-linear-models).

Generalized linear models
=========================

### Model

Suppose that *y*<sub>1</sub>, …, *y*<sub>*n*</sub> are observations on
independent random variables *Y*<sub>1</sub>, …, *Y*<sub>*n*</sub>, each
with probability density/mass function of the form
$$
f\_{Y\_i}(y) = \\exp\\left\\{\\frac{y \\theta\_i - b(\\theta\_i) - c\_1(y)}{\\phi/m\_i} - \\frac{1}{2}a\\left(-\\frac{m\_i}{\\phi}\\right) + c\_2(y) \\right\\}
$$
 for some sufficiently smooth functions *b*(.), *c*<sub>1</sub>(.),
*a*(.) and *c*<sub>2</sub>(.), and fixed observation weights
*m*<sub>1</sub>, …, *m*<sub>*n*</sub>. The expected value and the
variance of *Y*<sub>*i*</sub> are then
Hence, in this parameterization, *ϕ* is a dispersion parameter.

A generalized linear model links the mean *μ*<sub>*i*</sub> to a linear
predictor *η*<sub>*i*</sub> as
$$
g(\\mu\_i) = \\eta\_i = \\sum\_{t=1}^p \\beta\_t x\_{it}
$$
 where *g*(.) is a monotone, sufficiently smooth link function, taking
values on ℜ, *x*<sub>*i**t*</sub> is the (*i*, *t*)*t**h* component of a
model matrix *X*, and
*β* = (*β*<sub>1</sub>, …, *β*<sub>*p*</sub>)<sup>⊤</sup>.

### Score functions and information matrix

The derivatives of the log-likelihood about *β* and *ϕ* (score
functions) are
with *y* = (*y*<sub>1</sub>, …, *y*<sub>*n*</sub>)<sup>⊤</sup>,
*μ* = (*μ*<sub>1</sub>, …, *μ*<sub>*n*</sub>)<sup>⊤</sup>,
$W = {\\rm diag}\\left\\{w\_1, \\ldots, w\_n\\right\\}$ and
$D = {\\rm diag}\\left\\{d\_1, \\ldots, d\_n\\right\\}$, where
*w*<sub>*i*</sub> = *m*<sub>*i*</sub>*d*<sub>*i*</sub><sup>2</sup>/*V*(*μ*<sub>*i*</sub>)
is the *i*th working weight, and
*d*<sub>*i*</sub> = *d**μ*<sub>*i*</sub>/*d**η*<sub>*i*</sub>.
Furthermore,
*q*<sub>*i*</sub> = −2*m*<sub>*i*</sub>{*y*<sub>*i*</sub>*θ*<sub>*i*</sub> − *b*(*θ*<sub>*i*</sub>)−*c*<sub>1</sub>(*y*<sub>*i*</sub>)}
and *ρ*<sub>*i*</sub> = *m*<sub>*i*</sub>*a*′( − *m*<sub>*i*</sub>/*ϕ*)
are the *i*th deviance residual (e.g. as implemented in the `dev.resid`
component of a `family` object) and its expectation, respectively.

The expected information matrix about *β* and *ϕ* is
$$
i(\\beta, \\phi) =
\\left\[
\\begin{array}{cc}
i\_{\\beta\\beta}(\\beta, \\phi) & 0\_p \\\\
0\_p^\\top & i\_{\\phi\\phi}(\\beta, \\phi)
\\end{array}
\\right\]
=
\\left\[
\\begin{array}{cc}
\\frac{1}{\\phi} X^\\top W X & 0\_p \\\\
0\_p^\\top & \\frac{1}{2\\phi^4}\\sum\_{i = 1}^n m\_i^2 a''(-m\_i/\\phi)
\\end{array}
\\right\]\\,,
$$
 where 0<sub>*p*</sub> is a *p*-vector of zeros.

### Maximum likelihood estimation

The maximum likelihood estimator of *β* and *ϕ* satisfies
$s\_\\beta(\\hat\\beta,\\hat\\phi) = 0\_p$ and
$s\_\\phi(\\hat\\beta, \\hat\\phi) = 0$.

### Bias-reducing adjusted score functions

Let
*A*<sub>*β*</sub>(*β*, *ϕ*)= − *i*<sub>*β*</sub>(*β*, *ϕ*)*b*<sub>*β*</sub>(*β*, *ϕ*)
and
*A*<sub>*ϕ*</sub>(*β*, *ϕ*)= − *i*<sub>*ϕ*</sub>(*β*, *ϕ*)*b*<sub>*ϕ*</sub>(*β*, *ϕ*),
where *b*<sub>*β*</sub>(*β*, *ϕ*) and *b*<sub>*ϕ*</sub>(*β*, *ϕ*) are
the first terms in the expansion of the bias of the maximum likelihood
estimator of the regression parameters *β* and dispersion *ϕ*,
respectively. The results in Firth (1993) can be used to show that the
solution of the adjusted score equations
results in estimators $\\tilde\\beta$ and $\\tilde\\phi$ with bias of
smaller asymptotic order than the maximum likelihood estimator.

The results in either Kosmidis and Firth (2009) or Cordeiro and
McCullagh (1991) can then be used to re-express the adjustments in forms
that are convenient for implementation. In particular, and after some
algebra the bias-reducing adjustments for generalized linear models are
where *ξ* = (*ξ*<sub>1</sub>, …, *ξ*<sub>*n*</sub>)<sup>*T*</sup> with
*ξ*<sub>*i*</sub> = −*h*<sub>*i*</sub>*d*<sub>*i*</sub>′/(2*d*<sub>*i*</sub>*w*<sub>*i*</sub>),
*d*<sub>*i*</sub>′=*d*<sup>2</sup>*μ*<sub>*i*</sub>/*d**η*<sub>*i*</sub><sup>2</sup>
and *h*<sub>*i*</sub> is the "hat" value for the *i*th observation (see,
e.g. `?hatvalues`).

Fitting algorithm in `brglmFit`
===============================

`brglmFit` implements a quasi Fisher scoring procedure for solving the
adjusted score equations
*s*<sub>*β*</sub>(*β*, *ϕ*)+*A*<sub>*β*</sub>(*β*, *ϕ*)=0<sub>*p*</sub>
and *s*<sub>*ϕ*</sub>(*β*, *ϕ*)+*A*<sub>*ϕ*</sub>(*β*, *ϕ*)=0. The
iteration consists of an outer loop and an inner loop that implements
step-halving. The algorithm is as follows:

### Input

-   *s*<sub>*β*</sub>(*β*, *ϕ*), *i*<sub>*β**β*</sub>(*β*, *ϕ*),
    *A*<sub>*β*</sub>(*β*, *ϕ*)
-   *s*<sub>*ϕ*</sub>(*β*, *ϕ*), *i*<sub>*ϕ**ϕ*</sub>(*β*, *ϕ*),
    *A*<sub>*ϕ*</sub>(*β*, *ϕ*)
-   Starting values *β*<sup>(0)</sup> and *ϕ*<sup>(0)</sup>
-   *ϵ* &gt; 0: tolerance for the *L*1 norm of the direction before
    reporting convergence
-   *M*: maximum number of halving steps that can be taken

### Output

-   $\\tilde\\beta$, $\\tilde\\phi$

### Iteration

*Initialize outer loop*

1.  *k* ← 0

2.  *υ*<sub>*β*</sub><sup>(0)</sup> ← {*i*<sub>*β**β*</sub>(*β*<sup>(0)</sup>,*ϕ*<sup>(0)</sup>)}<sup>−1</sup>{*s*<sub>*β*</sub>(*β*<sup>(0)</sup>,*ϕ*<sup>(0)</sup>)+*A*<sub>*β*</sub>(*β*<sup>(0)</sup>,*ϕ*<sup>(0)</sup>)}

3.  *υ*<sub>*ϕ*</sub><sup>(0)</sup> ← {*i*<sub>*ϕ**ϕ*</sub>(*β*<sup>(0)</sup>,*ϕ*<sup>(0)</sup>)}<sup>−1</sup>{*s*<sub>*ϕ*</sub>(*β*<sup>(0)</sup>,*ϕ*<sup>(0)</sup>)+*A*<sub>*ϕ*</sub>(*β*<sup>(0)</sup>,*ϕ*<sup>(0)</sup>)}

*Initialize inner loop*

1.  *m* ← 0

2.  *b*<sup>(*m*)</sup> ← *β*<sup>(*k*)</sup>

3.  *f*<sup>(*m*)</sup> ← *ϕ*<sup>(*k*)</sup>

4.  *v*<sub>*β*</sub><sup>(*m*)</sup> ← *υ*<sub>*β*</sub><sup>(*k*)</sup>

5.  *v*<sub>*ϕ*</sub><sup>(*m*)</sup> ← *υ*<sub>*ϕ*</sub><sup>(*k*)</sup>

6.  *d* ← |*v*<sub>*β*</sub><sup>(*m*)</sup>|<sub>1</sub> + |*v*<sub>*ϕ*</sub><sup>(*m*)</sup>|

*Update parameters*

1.  *b*<sup>(*m* + 1)</sup> ← *b*<sup>(*m*)</sup> + 2<sup>−*m*</sup>*v*<sub>*β*</sub><sup>(*m*)</sup>

2.  *f*<sup>(*m* + 1)</sup> ← *f*<sup>(*m*)</sup> + 2<sup>−*m*</sup>*v*<sub>*ϕ*</sub><sup>(*m*)</sup>

*Update direction*

1.  *v*<sub>*β*</sub><sup>(*m* + 1)</sup> ← {*i*<sub>*β**β*</sub>(*b*<sup>(*m* + 1)</sup>,*f*<sup>(*m* + 1)</sup>)}<sup>−1</sup>{*s*<sub>*β*</sub>(*b*<sup>(*m* + 1)</sup>,*f*<sup>(*m* + 1)</sup>)+*A*<sub>*β*</sub>(*b*<sup>(*m* + 1)</sup>,*f*<sup>(*m* + 1)</sup>)}

2.  *v*<sub>*ϕ*</sub><sup>(*m* + 1)</sup> ← {*i*<sub>*ϕ**ϕ*</sub>(*b*<sup>(*m* + 1)</sup>,*f*<sup>(*m* + 1)</sup>)}<sup>−1</sup>{*s*<sub>*ϕ*</sub>(*b*<sup>(*m* + 1)</sup>,*f*<sup>(*m* + 1)</sup>)+*A*<sub>*ϕ*</sub>(*b*<sup>(*m* + 1)</sup>,*f*<sup>(*m* + 1)</sup>)}

*Continue or break halving within inner loop*

1.  if *m* + 1 &lt; *M* and
    |*v*<sub>*β*</sub><sup>(*m* + 1)</sup>|<sub>1</sub> + |*v*<sub>*ϕ*</sub><sup>(*m* + 1)</sup>| &gt; *d*

    14.1. *m* ← *m* + 1

    14.2. GO TO 10

2.  else

    15.1. *β*<sup>(*k* + 1)</sup> ← *b*<sup>(*m* + 1)</sup>

    15.2. *ϕ*<sup>(*k* + 1)</sup> ← *f*<sup>(*m* + 1)</sup>

    15.3.
    *υ*<sub>*β*</sub><sup>(*k* + 1)</sup> ← *v*<sub>*b*</sub><sup>(*m* + 1)</sup>

    15.4.
    *υ*<sub>*ϕ*</sub><sup>(*k* + 1)</sup> ← *v*<sub>*f*</sub><sup>(*m* + 1)</sup>

*Continue or break outer loop*

1.  if *k* + 1 &lt; *K* and
    |*υ*<sub>*β*</sub><sup>(*k* + 1)</sup>|<sub>1</sub> + |*υ*<sub>*ϕ*</sub><sup>(*k* + 1)</sup>| &gt; *ϵ*

    16.1 *k* ← *k* + 1

    16.2. GO TO 4

2.  else

    17.1. $\\tilde\\beta \\leftarrow \\beta^{(k + 1)}$

    17.2. $\\tilde\\phi \\leftarrow \\phi^{(k + 1)}$

    17.3. STOP

Notes
=====

-   For *K* = *M* = 1, $\\beta^{(0)} = \\hat\\beta$ and
    $\\phi^{(0)} = \\hat\\phi$, the above iteration computes the
    bias-corrected estimates proposed in Cordeiro and McCullagh (1991).

-   The steps where *ϕ* and the *ϕ* direction are updated are ignored
    for generalized linear models with known dispersion parameter, like
    in models with binomial and poisson responses. Also, in that case,
    *v*<sub>*ϕ*</sub><sup>(.)</sup> and *υ*<sub>*ϕ*</sub><sup>(.)</sup>
    in steps 9, 14 and 16 are set to zero.

-   The implementation of the adjusted score functions requires ready
    implementations of
    *d*<sup>2</sup>*μ*<sub>*i*</sub>/*d**η*<sub>*i*</sub><sup>2</sup>,
    *a*′(.), *a*″(.) and *a*‴(.). The
    [**enrichwith**](https://cran.r-project.org/package=enrichwith) R
    package is used internally to enrich the base `family` and
    `link-glm` objects with implementations of those functions (see
    `?enrich.family` and `?enrich.link-glm`).

-   The above iteration can be used to implement a variety of additive
    adjustments to the score function, by supplying the algorithm with
    appropriate adjustment functions *A*<sub>*β*</sub>(*β*, *ϕ*) and
    *A*<sub>*ϕ*</sub>(*β*, *ϕ*)

References
==========

Albert, A., and J. A. Anderson. 1984. “On the Existence of Maximum
Likelihood Estimates in Logistic Regression Models.” *Biometrika* 71
(1): 1–10.

Cordeiro, G. M., and P. McCullagh. 1991. “Bias Correction in Generalized
Linear Models.” *Journal of the Royal Statistical Society, Series B:
Methodological* 53 (3): 629–43.

Firth, D. 1993. “Bias Reduction of Maximum Likelihood Estimates.”
*Biometrika* 80 (1): 27–38.

Kosmidis, I. 2014. “Bias in Parametric Estimation: Reduction and Useful
Side-Effects.” *Wiley Interdisciplinary Reviews: Computational
Statistics* 6 (3). John Wiley & Sons, Inc.: 185–96.
doi:[10.1002/wics.1296](https://doi.org/10.1002/wics.1296).

Kosmidis, I., and D. Firth. 2009. “Bias Reduction in Exponential Family
Nonlinear Models.” *Biometrika* 96 (4): 793–804.
doi:[10.1093/biomet/asp055](https://doi.org/10.1093/biomet/asp055).

———. 2010. “A Generic Algorithm for Reducing Bias in Parametric
Estimation.” *Electronic Journal of Statistics* 4: 1097–1112.
doi:[10.1214/10-EJS579](https://doi.org/10.1214/10-EJS579).
