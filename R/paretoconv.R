#' Convolution of Pareto distributions 
#'
#' Convolutes multiple Pareto distributions following 
#' \itemize{
#' \item For integer shape parameters: 'The Distribution of Sums of Certain
#' I.I.D.  Pareto Variates' by Colin Ramsay (Communications in Statistics -
#' Theory and Methods 35:395-405, 2006);
#' \item For non-integer shape parameters: 'The Distribution of Sums of I.I.D.
#' Pareto Random Variables with Arbitrary Shape Parameter' by Colin Ramsay
#' (Communications in Statistics - Theory and Methods 37:2177-2184, 2008).
#' }
#'
#' @param x value of independent variable (may be a vector)
#' @param a The primary shape parameter of the Pareto distribution (single value
#' only)
#' @param n Number of convolutions (single value only)
#' @param x0 Lower cut-off point of classical heavy-tailed distribution
#' (generally obtained emprically with the poweRlaw package).
#' @param cdf If TRUE, returns the cumulative distribution function, otherwise
#' returns the probability density function.
#' @param asymp If TRUE and \code{x} is a vector of 10 or more elements,
#' asymptotic convergence is checked for large \code{x}, with values beyond
#' convergence replaced with directly power-law values
#' @param quiet If FALSE, issue progress messages
#'
#' @note The Pareto distribution may be defined as f(x)=(a/b)(b/x)^(a-1), where
#' a and b are the primary and secondary shape parameters, respectively. It
#' presumed here without loss of generality that b=1 and thus f(x)=a x^(1-a).
#' Convolution of multiple distritions (that is, \code{n>0}) are NOT normalised,
#' so CDFs to not sum to unity, and PDFs do not integrate to unity.
#'
#' @return Value for the CDF or PDF from the convolution of two Pareto
#' distributions of shape a at the value x.
#'
#' @export
#' @examples 
#' paretoconv (x=1:10, a=2, n=0)
paretoconv <- function (x, a, n, x0=1, cdf=FALSE, asymp=TRUE, quiet=TRUE)
{
    if (missing (x)) stop ('x must be supplied')
    if (missing (a)) stop ('a must be supplied')
    if (missing (n)) stop ('n must be supplied')
    if (!is.numeric (x)) stop ('x must be numeric')
    if (!is.numeric (a)) stop ('a must be numeric')
    if (!is.numeric (n)) stop ('n must be numeric')
    if (length (a) > 1) stop ('a must be a single value only')
    if (length (n) > 1) stop ('n must be a single value only')

    if (any (x < 0)) stop ('x must be positive')
    if (n%%1 != 0) stop ('n must be an integer')
    if (n < 0) stop ('n must be a non-negative integer')

    fn <- "" # supress R CMD check warnings
    if (n == 0)
    {
        if (cdf)
            y <- (x0 / x) ^ (a - 1)
        else
            y <- (a / x0) * (x / x0) ^ (-a)
    } else
    {
        a <- a - 1 # Ramsay transforms to this value
        if (a%%1 == 0) 
        {
            if (cdf)
                fn <- "ramsay_int_cdf"
            else
                fn <- "ramsay_int_pdf"
        } else
        {
            if (cdf)
                fn <- "ramsay_nonint_cdf"
            else
                fn <- "ramsay_nonint_pdf"
        }
        # The value of 10 is arbitrary, but serves to implement aymptotic
        # approximations for large x in longer vectors.
        if (length (x) < 10 | !asymp)
            y <- sapply (x, function (i) 
                         do.call (fn, list (x=i, a=a, n=n, x0=x0)))
        else
        {
            # aymptotic approximation
            y <- asympt (x=x, fn=fn, a=a, n=n, x0=x0, quiet=quiet)
        }
    }

    return (y)
}

#' Asymptotic approximation for large x
#'
#' @param x Vector of values of independent variable
#' @param fn Name of function (int/nonint, pdf/cdf) as set in \code{paretoconv}
#' @param a The primary shape parameter of the Pareto distribution (single value
#' only)
#' @param n Number of convolutions (single value only)
#' @param x0 Lower cut-off point of classical heavy-tailed distribution
#' (generally obtained emprically with the poweRlaw package).
#' @param tol Convergence tolerance (values around 0.01 are typically
#' sufficient; values < 0.001 may take a long time to converge).
#' @param quiet If FALSE, issue progress messages
#'
#' @note The Pareto distribution may be defined as f(x)=(a/b)(b/x)^(a-1), where
#' a and b are the primary and secondary shape parameters, respectively. It
#' presumed here without loss of generality that b=1 and thus f(x)=a x^(1-a).
#'
#' @return Value for the CDF or PDF from the convolution of two Pareto
#' distributions of shape a at the value x.
asympt <- function (x, fn, a, n, x0, tol=0.05, quiet=quiet)
{
    indx <- sort (x, index.return=TRUE)$ix
    x <- x [indx]
    y <- rep (NA, length (x))
    y [1] <- do.call (fn, list (x=x[1], a=a, x0=x0, n=n))
    dy <- 1
    alph <- a + 1 # because 1 has already been subtracted above
    if (grepl ("cdf", fn))
        alph <- a
    msg <- TRUE
    for (i in 2:length (x))
    {
        if (dy > tol)
            y [i] <- do.call (fn, list (x=x[i], a=a, x0=x0, n=n))
        else # asymptotic approximation
            y [i] <- y [i-1] * (x [i] / x [i-1]) ^ (-alph)

        dx <- x [i] / x [i-1]
        if (x [i-1] == 0) dx <- 1 # in case x [1] == 0
        if (y [i] == 0)
            dy <- 0
        else
            dy <- abs (y [i] / y [i-1] - dx ^ (-alph)) / (dx - 1)

        if (!quiet & msg)
            message ("\rnot converged to asymptotic tail at x [", i-1, 
                     "/", length (x), "] = ", x [i], "; dy = ", 
                     formatC (dy, format="f", digits=3), " > ", tol, appendLF=FALSE)
        if (dy < tol)
        {
            if (!quiet & msg)
                message ("\rconverged to asymptotic tail at x [", i-1, 
                         "/", length (x), "] = ", x [i], "; dy = ", 
                         formatC (dy, format="f", digits=3), " < ", tol, "             ")
            msg <- FALSE
        }
    }
    if (!quiet & msg) message ("") # line feed
    y [indx]
}

#' Get limit of asymptotic convergence
#'
#' @param a The primary shape parameter of the Pareto distribution (single value
#' only)
#' @param n Number of convolutions (single value only)
#' @param x0 Lower cut-off point of classical heavy-tailed distribution
#' (generally obtained emprically with the poweRlaw package).
#' @param tol Convergence tolerance (values around 0.05 are typically
#' sufficient; values < 0.01 will generally not work).
#' @param quiet If FALSE, issue progress messages
#'
#' @return Single integer value of limit beyond which PDF may be approximated by
#' its asymptotic power law
#' @export
get_asymp_limit <- function (a, n, x0, tol=0.05, quiet=FALSE)
{
    # NOTE: The point of asymptotic convergence is approached by multiplying
    # each value of `x` by `(1 + dy)`, so that `x` increases by less as the
    # convergence point is approached. This is necessary because increasing `x`
    # too rapidly leads to `paretoconv()` returning `y=0`, and so `dy` becomes
    # undefined.
    x <- 1:2
    dy <- 1
    while (dy > tol)
    {
        y <- paretoconv (x, a=alpha, x0=x0, n=n, cdf=FALSE)
        dx <- x [2] / x [1]
        dy <- abs (y [2] / y [1] - dx ^ (-a)) / (dx - 1)
        if (any (y == 0)) dy <- 0
        if (!quiet & dy > tol)
            message ('\rdy [x = ', round (x [1]), '] = ', dy, ' > ', tol,
                     appendLF=FALSE)
        x [1] <- x [1] * (1 + dy)
        x [2] <- x [1] + 1
    }
    if (!quiet)
        message ('\rdy [x = ', round (x [1]), '] = ', dy, ' < ', tol)
    return (x [1])
}
