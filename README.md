[![Build Status](https://travis-ci.org/mpadge/paretoconv.svg)](https://travis-ci.org/mpadge/paretoconv) [![Project Status: Concept - Minimal or no implementation has been done yet.](http://www.repostatus.org/badges/0.1.0/concept.svg)](http://www.repostatus.org/#concept) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/paretoconv)](http://cran.r-project.org/web/packages/paretoconv) [![codecov](https://codecov.io/gh/mpadge/paratoconv/branch/master/graph/badge.svg)](https://codecov.io/gh/mpadge/paretoconv)

`paretoconv` contains only one function:

    paretoconv (x, a, n)

which calculates the *n*-fold convolution of two Pareto distributions,
\begin{equation}
f(x)=a x^{1-a},
\end{equation}
for (*x*, *a*)&gt;0 using the techniques devised by Colin Ramsay in

1.  'The Distribution of Sums of Certain I.I.D. Pareto Variates' (Communications in Statistics - Theory and Methods 35:395-405, 2006); and

2.  'The Distribution of Sums of I.I.D. Pareto Random Variables with Arbitrary Shape Parameter' (Communications in Statistics - Theory and Methods 37:2177-2184, 2008).

### News

-   Initial working version

### Installation

``` r
devtools::install_github("mpadge/paretoconv")
```

### Usage

``` r
library(paretoconv)
# current verison
packageVersion("paretoconv")
#> [1] '0.0.0'
paretoconv (x=0:5, a=1, n=5)
#> [1] 0.000000000 0.001966075 0.020217706 0.060057727 0.113677829 0.172729184
```

Note that `x` can be a vector, while both `a` and `n` must single-valued.

### Test Results

``` r
library(paretoconv)
library(testthat)

date()

test_dir("tests/")
```
