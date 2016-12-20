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
#'
#' @note The Pareto distribution may be defined as f(x)=(a/b)(b/x)^(a-1), where
#' a and b are the primary and secondary shape parameters, respectively. It
#' presumed here without loss of generality that b=1 and thus f(x)=a x^(1-a).
#'
#' @return Value for the CDF or PDF from the convolution of two Pareto
#' distributions of shape a at the value x.
#'
#' @export
#' @examples 
#' paretoconv (1:10, 1, 1)
paretoconv <- function (x, a, n, x0=1, cdf=FALSE)
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
        y <- sapply (x, function (i) do.call (fn, list (x=i, a=a, n=n, x0=x0)))
    }

    return (y)
}
