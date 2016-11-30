#' Ramsay's (2006) chi function of Eq. (8)
#'
#' @param x value of independent variable
#' @param a The primary shape parameter of the Pareto distribution - alpha in
#' Ramsay's notation. 
#' @param n Number of convolutions
#'
#' @return Single value of chi
chi <- function (x, a, n)
{
    # Kummer's M function, the confluent hypergeometric _1F_1
    kM <- function (a, b, x) fAsianOptions::kummerM (a, b, x) # Eq. (13)
    # NOTE: cot (z) = 1 / tan (z)
    # Rxa is Eq. (12)
    Rxa <- function (x, a) 
        kM (1, 1-a, -x) - pi * x ^ a * exp (-x) / (tan (pi * a) * gamma (a))
    Ixa <- function (x, a) pi * x ^ a * exp (-x) / gamma (a) # Eq. (6)

    # Then the value of chi from Eq. (8)
    s1 <- function (x, a, n, r)
    {
        Rn <- Rxa (x, a) ^ (n - 2 * r - 1)
        In <- Ixa (x, a) ^ (2 * r + 1)
        (-1) ^ r * choose (n, 2 * r + 1) * Rn * In / pi
    }
    nseq <- 0:floor ((n - 1) / 2)
    sum (sapply (nseq, function (i) s1 (x, a, n, i)))
}

#' CDF for Convolution of Pareto distributions for non-integer alpha
#'
#' Convolutes multiple Pareto distributions following 'The Distribution of Sums
#' of I.I.D. Pareto Random Variables with Arbitrary Shape Parameter' by Colin
#' Ramsay (Communications in Statistics - Theory and Methods 37:2177-2184,
#' 2008).
#'
#' This evaluates Eq. (9) via Eqs. (6), (8), and (12).
#'
#' @param x value of independent variable
#' @param a The primary shape parameter of the Pareto distribution - alpha in
#' Ramsay's notation. 
#' @param n Number of convolutions
#'
#' @note This version is for non-integer values of a. Integer values can be
#' calculated with ramsay_int
#'
#' @return Value for the CDF of the convolution of two Pareto distributions of
#' shape a at the value x.
ramsay_nonint_cdf <- function (x, a, n)
{
    bet <- set_beta () # Always = 1

    # In this integrand, y is ramsay's x and x is his t, so his F_n(t) is here
    # F_n(x), and the integral is over y-values
    integrand <- function (y, x, a, n) 
        (1 - exp (-x * y / bet)) * Re (chi (y, a, n)) / y
    cubature::adaptIntegrate (integrand, lowerLimit=0, upperLimit=1e3,
                              x=x, a=a, n=n)$integral
}

#' PDF for Convolution of Pareto distributions for non-integer alpha
#'
#' Convolutes multiple Pareto distributions following 'The Distribution of Sums
#' of I.I.D. Pareto Random Variables with Arbitrary Shape Parameter' by Colin
#' Ramsay (Communications in Statistics - Theory and Methods 37:2177-2184,
#' 2008).
#'
#' This evaluates Eq. (7) via Eqs. (6), (8), and (12).
#'
#' @param x value of independent variable
#' @param a The primary shape parameter of the Pareto distribution - alpha in
#' Ramsay's notation. 
#' @param n Number of convolutions
#'
#' @note This version is for non-integer values of a. Integer values can be
#' calculated with ramsay_int
#'
#' @return Value for the PDF of the convolution of two Pareto distributions of
#' shape a at the value x.
ramsay_nonint_pdf <- function (x, a, n)
{
    bet <- set_beta () # Always = 1

    # In this integrand, y is ramsay's x and x is his t, so his F_n(t) is here
    # F_n(x), and the integral is over y-values
    integrand <- function (y, x, a, n)
        exp (-x * y / bet) * Re (chi (y, a, n))
    cubature::adaptIntegrate (integrand, lowerLimit=0, upperLimit=1e3,
                              x=x, a=a, n=n)$integral / bet
}
