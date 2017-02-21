#' Ramsay's (2006) chi function of Eq. (8)
#'
#' @param z value of independent variable
#' @param a The primary shape parameter of the Pareto distribution - alpha in
#' Ramsay's notation. 
#' @param n Number of convolutions
#'
#' @return Single value of chi
chi <- function (z, a, n)
{
    # Kummer's M function, the confluent hypergeometric _1F_1. NOTE:
    # fAsianOptions reverses the order of the arguments from Abramowitz &
    # Stegun!
    kM <- function (a, b, z) fAsianOptions::kummerM (z, a, b) # Eq. (13)
    # NOTE: cot (z) = 1 / tan (z)
    RI <- function (z, a, n, r)
    {
        Iza <- pi * z ^ a * exp (-z) / gamma (a) # Eq. (6)
        Rza <- kM (1, 1-a, -z) - Iza / tan (pi * a) # Eq. (12)
        Iza ^ (2 * r + 1) * Rza ^ (n - 2 * r - 1) # part of Eq. (8)
    }
    # Then the value of chi from Eq. (8)
    s1 <- function (z, a, n, r)
        (-1) ^ r * choose (n, 2 * r + 1) * RI (z, a, n, r) / pi
    nseq <- 0:floor ((n - 1) / 2)
    sum (sapply (nseq, function (i) s1 (z, a, n, i)))
}


#' CDF for Convolution of Pareto distributions for non-integer alpha
#'
#' Calculates complementary Cumulative Distribution Function (cdf) from
#' convolution of multiple Pareto distributions following 'The Distribution of
#' Sums of I.I.D. Pareto Random Variables with Arbitrary Shape Parameter' by
#' Colin Ramsay (Communications in Statistics - Theory and Methods 37:2177-2184,
#' 2008).
#'
#' This evaluates Eq. (9) via Eqs. (6), (8), and (12).
#'
#' @param x value of independent variable
#' @param a The primary shape parameter of the Pareto distribution - alpha in
#' Ramsay's notation. 
#' @param n Number of convolutions
#' @param x0 Lower cut-off point of classical heavy-tailed distribution
#' (generally obtained emprically with the poweRlaw package).
#'
#' @note This version is for non-integer values of a. Integer values can be
#' calculated with ramsay_int
#'
#' @return Value for the CDF of the convolution of two Pareto distributions of
#' shape a at the value x.
ramsay_nonint_cdf <- function (x, a, n, x0)
{
    # In this integrand, z is ramsay's x and x is his t, so his F_n(t) is here
    # F_n(x), and the integral is over z-values
    integrand <- function (z, x, a, x0, n) 
    {
        if (z == 0)
            0
        else
            (1 - exp (-x * z / x0)) * Re (chi (z, a, n)) / z
    }
    if (x == 0)
        1
    else
        1 - calc_integral (integrand, x, a, x0, n)
}

#' PDF for Convolution of Pareto distributions for non-integer alpha
#'
#' Calculates complementary Cumulative Distribution Function (cdf) from
#' convolution of multiple Pareto distributions following 'The Distribution of
#' Sums of I.I.D. Pareto Random Variables with Arbitrary Shape Parameter' by
#' Colin Ramsay (Communications in Statistics - Theory and Methods 37:2177-2184,
#' 2008).
#'
#' This evaluates Eq. (7) via Eqs. (6), (8), and (12).
#'
#' @param x value of independent variable
#' @param a The primary shape parameter of the Pareto distribution - alpha in
#' Ramsay's notation. 
#' @param n Number of convolutions
#' @param x0 Lower cut-off point of classical heavy-tailed distribution
#' (generally obtained emprically with the poweRlaw package).
#'
#' @note This version is for non-integer values of a. Integer values can be
#' calculated with ramsay_int
#'
#' @return Value for the PDF of the convolution of two Pareto distributions of
#' shape a at the value x.
ramsay_nonint_pdf <- function (x, a, n, x0)
{
    # In this integrand, z is ramsay's x and x is his t, so his F_n(t) is here
    # F_n(x), and the integral is over z-values
    integrand <- function (z, x, a, x0, n)
        exp (-x * z / x0) * Re (chi (z, a, n))
    if (x == 0)
        0
    else
        calc_integral (integrand, x, a, x0, n)
}
