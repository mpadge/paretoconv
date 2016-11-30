#' Convolution of Pareto distributions for non-integer alpha
#'
#' Convolutes multiple Pareto distributions following 'The Distribution of Sums
#' of I.I.D. Pareto Random Variables with Arbitrary Shape Parameter' by Colin
#' Ramsay (Communications in Statistics - Theory and Methods 37:2177-2184,
#' 2008).
#'
#' @param x value of independent variable
#' @param a The primary shape parameter of the Pareto distribution - alpha in
#' Ramsay's notation. 
#' @param n Number of convolutions
#'
#' @note This version is for non-integer values of a. Integer values can be
#' calculated with ramsay_int
#'
#' @return Value for the convolution of two Pareto distributions of shape a at
#' the value x.
ramsay_nonint <- function (x, a, n)
{
    # NOTE: compare ramsay_nonint
    if (a%%1 == 0) stop ("a can not be an integer")

    bet <- 1 # beta value fixed for my purposes always at 1

    # Kummer's M function, the confluent hypergeometric _1F_1
    kM <- function (a, b, x) fAsianOptions::kummerM (a, b, x)
    # NOTE: cot (z) = 1 / tan (z)
    Rxa <- function (x, a) 
        kM (1, 1-a, -x) - pi * x ^ a * exp (-x) / (tan (pi * a) * gamma (a))
    Ixa <- function (x, a) pi * x ^ a * exp (-x) / gamma (a)
    chi <- function (x, a, n)
    {
        s1 <- function (x, a, n, r)
        {
            Rn <- Rxa (x, a) ^ (n - 2 * r - 1)
            In <- Ixa (x, a) ^ (2 * r + 1)
            (-1) ^ r * choose (n, 2 * r + 1) * Rn * In / pi
        }
        nseq <- 0:floor ((n - 1) / 2)
        sum (sapply (nseq, function (i) s1 (x, a, n, i)))
    }
    # In this integrand, y is ramsay's x and x is his t, so his F_n(t) is here
    # F_n(x), and the integral is over y-values
    integrand <- function (y, x, a, n) 
        (1 - exp (-x * y / bet)) * Re (chi (y, a, n)) / y
    cubature::adaptIntegrate (integrand, lowerLimit=0, upperLimit=1e3,
                              x=x, a=a, n=n)$integral
}
