[![Build Status](https://travis-ci.org/mpadge/paretoconv.svg)](https://travis-ci.org/mpadge/paretoconv) [![Build status](https://ci.appveyor.com/api/projects/status/github/mpadge/paretoconv?svg=true)](https://ci.appveyor.com/project/mpadge/paretoconv) [![codecov](https://codecov.io/gh/mpadge/paretoconv/branch/master/graph/badge.svg)](https://codecov.io/gh/mpadge/paretoconv) [![Project Status: WIP](http://www.repostatus.org/badges/0.1.0/wip.svg)](http://www.repostatus.org/#wip) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/paretoconv)](http://cran.r-project.org/web/packages/paretoconv)

paretoconv
==========

An `R` package to calculate the *n*-fold convolution of Pareto distributions. *f(x)=a x<sup>1-a</sup>*, for *x&gt;0*, *a&gt;1* using the techniques devised by Colin Ramsay in

1.  'The Distribution of Sums of Certain I.I.D. Pareto Variates' (*Communications in Statistics - Theory and Methods* **35**:395-405, 2006); and

2.  'The Distribution of Sums of I.I.D. Pareto Random Variables with Arbitrary Shape Parameter' (*Communications in Statistics - Theory and Methods* **37**:2177-2184, 2008).

The package contains only one function:

    paretoconv (x, a, n, cdf=FALSE)

where `n` specifies the number of convolutions. Both this and `a` must be single-valued, while `x` can be a vector. `cdf` generates the cumulative distribution function, otherwise the probability density function is returned.

------------------------------------------------------------------------

### News

-   Initial working version

### Installation

``` r
devtools::install_github("mpadge/paretoconv")
```

### Usage

``` r
library(paretoconv)
packageVersion("paretoconv") 
#> [1] '0.0.0'
```

------------------------------------------------------------------------

Example
-------

Solid lines in the figure below are a reproduction of Ramsay's (2006) Figure 2 of probability density functions for the first 5 convolutions of the Pareto pdf with shape parameter of *a=5*. Dashed lines are analogous values for the non-integer value of *a=4.5*.

``` r
x <- 1:50 / 10
n <- 1:5
yint <- lapply (n, function (i) paretoconv (x=x, a=5, n=i))
ynon <- lapply (n, function (i) paretoconv (x=x, a=4.5, n=i))
cols <- rainbow (length (n))
plot (NULL, NULL, xlim=range (x), ylim=range (yint, na.rm=TRUE), xlab="x", ylab="p")
for (i in n) {
    lines (x, yint [[i]], col=cols [i])
    lines (x, ynon [[i]], col=cols [i], lty=2)
}
legend ("topright", lwd=1, col=cols, bty="n", 
        legend=sapply (seq (n), function (i) paste0 ("n=", i)))
```

![](./fig/README-example.png)

------------------------------------------------------------------------

### Test Results

``` r
date()
#> [1] "Wed Dec 21 11:29:55 2016"
testthat::test_dir("tests/")
#> testthat results ========================================================================================================
#> OK: 21 SKIPPED: 0 FAILED: 0
#> 
#> DONE ===================================================================================================================
```
