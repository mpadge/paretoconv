(Note that the LaTeX equations are not rendered here but are kept that way for ready copy-paste into other rendering environments.)

``` r
devtools::load_all (".")
```

    ## Loading paretoconv

------------------------------------------------------------------------

If the data converge to the underlying power-law relationship, then the expected CDF at some point *x* × *Δ**x* following an observed value of *y*(*x*) will be,
\begin{equation}
    y = A \left(\frac{x}{x_0}\right)^{-\alpha+1},\ {\rm and}\, 
        y^\prime = A \left(\frac{x\cdot\Delta x}{x_0}\right)^{-\alpha+1},
\end{equation}
so that
\begin{equation}
    \frac{y^\prime}{y} = \Delta x^{-\alpha+1}.
\end{equation}
The equivalent PDF will be,
\begin{equation}
    y = \frac{a-1}{x_0} \left(\frac{x}{x_0}\right)^{-\alpha},\ {\rm and}\, 
        y^\prime = \frac{a-1}{x_0} \left(\frac{x\cdot\Delta x}{x_0}\right)^{-\alpha},
\end{equation}
so that the ratio of consecutive values in this case is,
\begin{equation}
    \frac{y^\prime}{y} = \Delta x^{-\alpha}.
\end{equation}
Raw Values
----------

The raw values look like this:

``` r
dx <- 1.1
x <- cumprod (rep (dx, 50))
x <- x [x > 1 & x < 100]
a <- 2.2
x0 <- 2
n <- 2
y0 <- (x / x0) ^ (-a+1)
y1 <- 1 - paretoconv (x=x, a=a, x0=x0, n=n, cdf=TRUE)
plot (x, y0, "l", col="blue", log="xy", ylim=range (c (y0, y1)))
lines (x, y1, col="red")
points (x, y1, pch=1, col="red")
```

![](https://github.com/mpadge/paretoconv/blob/master/fig/wiki1-raw-values-1.png)

This shows that the lines themselves don't necessarily converge, yet the gradients do, and thus convergence tests can only be performed on the gradients.

Gradients
---------

The ratios of consectuve points for `dx=1.1` are:

``` r
indx <- 2:length (x)
y0r <- y0 [indx] / y0 [indx - 1]
y1r <- y1 [indx] / y1 [indx - 1]
plot (x [indx - 1], y0r, "l", col="blue", ylim=range (c (y0r, y1r)), xlab="x")
lines (x [indx - 1], y1r, col="red")
```

![](https://github.com/mpadge/paretoconv/blob/master/fig/wiki1-ratios-1.png)

Both lines converge to the expected ratio of:

``` r
dx ^ (-a + 1)
## [1] 0.8919259
```

``` r
tail (y0r, n=1); tail (y1r, n=1)
## [1] 0.8919259
## [1] 0.8905734
```

------------------------------------------------------------------------

Convergence Loop
----------------

A convergence loop looks something like this:

``` r
convloop <- function (dx=1.1, tol=0.05, xs=10, a=2.2, x0=2, n=2, quiet=TRUE)
{
    yold <- 1 - paretoconv (x=xs, a=a, x0=x0, n=n, cdf=TRUE)
    y <- 1 - paretoconv (x=xs * dx, a=a, x0=x0, n=n, cdf=TRUE)
    dy0 <- dx ^ (-a + 1) # expected value
    dy <- abs (y / yold - dy0) / (dx - 1)

    count <- 1
    while (dy > tol)
    {
        yold <- y
        xs <- xs * dx
        y <- 1 - paretoconv (x=xs, a=a, x0=x0, n=n, cdf=TRUE)
        dy <- abs (y / yold - dy0) / (dx - 1)
        if (!quiet)
            message ("[", count, "] (x,dy)=(", xs, ", ", dy, ")")
        count <- count + 1
    }
    # return x as mean of final two values
    list (count=count, x=mean (c (xs, xs / dx)))
}
```

Different values of `dx` give the following results

``` r
convloop (dx=1.01)
## $count
## [1] 92
## 
## $x
## [1] 24.60876
```

``` r
convloop (dx=1.05)
## $count
## [1] 20
## 
## $x
## [1] 24.66785
```

``` r
convloop (dx=1.1)
## $count
## [1] 11
## 
## $x
## [1] 24.75845
```

``` r
convloop (dx=1.5)
## $count
## [1] 3
## 
## $x
## [1] 18.75
```

``` r
convloop (dx=2)
## $count
## [1] 1
## 
## $x
## [1] 7.5
```

This suggests in turn that one could start with large values of `dx` and decrease them until `count` exceeds, say, 5 or so.

Convergence of dx
-----------------

A second try:

``` r
convloop2 <- function (tol=0.05, xs=10, a=2.2, x0=2, n=2, quiet=TRUE)
{
    count <- 1
    dx <- 2
    while (count < 5)
    {
        res <- convloop (dx=dx, tol=tol, xs=xs, a=a, x0=x0, n=n, quiet=quiet)
        count <- res$count
        if (!quiet)
            message ("dx=", dx, " (count, x) = (", count, ", ", res$x, ")")
        dx <- sqrt (dx)
    }
    list (dx=dx^2, n=count, x=res$x)
}
convloop2 ()
## $dx
## [1] 1.189207
## 
## $n
## [1] 6
## 
## $x
## [1] 21.89207
```

The relationship between `tol` and resultant estimates of `x` is as expected:

``` r
tol <- 5:1 / 100
res <- sapply (tol, function (i) {
                   pt0 <- proc.time ()[3]
                   c (convloop2 (tol=i)$x, proc.time ()[3] - pt0)
                })
res <- data.frame (x=res[1,], t=res[2,])
```

Then plot the results

``` r
cols <- c ("blue", "red")
plot.new ()
par (mar=c(5.1,4.1,2.1,4.1))
plot (tol, res$x, "l", ylab="x", col=cols [1])
points (tol, res$x, pch=1, col=cols [1])
par ("new"=TRUE)
plot (tol, res$t, "l", col=cols [2], xaxt="n", yaxt="n", xlab="", ylab="")
points (tol, res$t, pch=1, col=cols [2])
Axis (res$t, side=4)
legend ("topright", lwd=1, col=cols, bty="n", legend=c("x", "time"))
```

![](https://github.com/mpadge/paretoconv/blob/master/fig/wiki1-convplot-1.png)

------------------------------------------------------------------------

American civil war example
--------------------------

Values for the American civil war data analysed by Colin Gillespie are:

``` r
x0 <- 20
alpha <- 2.21
```

The point of convergence can be estimated with the above routines like this:

``` r
pt0 <- proc.time ()
convloop2 (tol=0.05, xs=50, a=alpha, x0=x0, n=1, quiet=FALSE)
proc.time () - pt0
## $dx
## [1] 1.414214
##
## $n
## [1] 7
##
## $x
## [1] 341.4214
##   user  system elapsed 
## 71.516   0.000  71.516 
```

And that takes quite a long time (over one minute) but gets there in the end.

What is needed, however, is for convergence to measured for each given value of `x`. This is illustrated here by dividing the values for native American casualties by 10 to enable calculation in a reasonable time (and it is presumed for convenience that `x` values are sorted):

``` r
data ("native_american", package="poweRlaw")
x <- sort (unique (native_american$Cas)) / 10
x0 <- 20 / 10
alpha <- 2.21
fity <- function (x, x0, alpha, tol=0.01, quiet=TRUE)
{
    y <- rep (NA, length (x))
    y [1:2] <- 1 - paretoconv (x=x [1:2], a=alpha, x0=x0, n=1, cdf=TRUE)
    dy0 <- (y [2] / y [1])
    dy <- abs (y [2] / y [1] - (x [2] / x [1]) ^ (-a + 1)) / (x [2] / x [1] - 1)
    i <- 3
    yold <- y [2]
    while (dy > tol)
    {
        y [i] <- 1 - paretoconv (x=x [i], a=alpha, x0=x0, n=1, cdf=TRUE)
        dy0 <- (y [i] / y [i-1])
        dx <- x [i] / x [i-1] - 1
        dy <- abs (y [i] / y [i-1] - (x [i] / x [i-1]) ^ (-a + 1)) / dx
        if (!quiet)
            message ("[", i, "] (x, dx, dy) = (", x [i], ", ", dx, ", ", dy, ")")
        i <- i + 1
    }
    x [i-1]
}
```

``` r
pt0 <- proc.time ()
fity (x=x, x0=x0, alpha=alpha)
proc.time () - pt0
## [1] 100
##   user  system elapsed
## 16.696   0.000  16.698
```

Replacing `paretoconv` values with asymptotic equivalents
---------------------------------------------------------

What is finally needed is to replace the values for `x>=100` with asymptotic estimates.

``` r
x3 <- tail (x, n=3)
y3 <- 1 - paretoconv (x=x3 [1], a=alpha, x0=x0, n=1, cdf=TRUE)
y3 <- c (y3, NA, NA)
y3 [2] <- y3 [1] * (x3 [2] / x3 [1]) ^ (-alpha + 1)
y3 [3] <- y3 [2] * (x3 [3] / x3 [2]) ^ (-alpha + 1)
```

Calculate the preceding values with `paretoconv`

``` r
x1 <- x [1:(length (x) - 3)]
y1 <- 1 - paretoconv (x=x1, a=alpha, x0=x0, n=1, cdf=TRUE)
```

And finally combine the two on the same plot:

``` r
plot (x1, y1, "l", col="blue", log="xy", 
      xlim=range (c (x1, x3)), ylim=range (c (y1, y3)))
lines (x3, y3, col="red", lwd=2, lty=2)
```

![](https://github.com/mpadge/paretoconv/blob/master/fig/wiki1-native-american-tail-plot-1.png)

Equivalent values for the PDF look like this:

``` r
x3 <- tail (x, n=3)
y3 <- paretoconv (x=x3 [1], a=alpha, x0=x0, n=1, cdf=FALSE)
y3 <- c (y3, NA, NA)
y3 [2] <- y3 [1] * (x3 [2] / x3 [1]) ^ (-alpha)
y3 [3] <- y3 [2] * (x3 [3] / x3 [2]) ^ (-alpha)
```

Calculate the preceding values with `paretoconv`

``` r
x1 <- x [1:(length (x) - 3)]
y1 <- paretoconv (x=x1, a=alpha, x0=x0, n=1, cdf=FALSE)
```

And finally combine the two on the same plot:

``` r
plot (x1, y1, "l", col="blue", log="xy", 
      xlim=range (c (x1, x3)), ylim=range (c (y1, y3)))
lines (x3, y3, col="red", lwd=2, lty=2)
```

![](aaa_files/figure-markdown_github/native-american-tail-plot-pdf-1.png)
