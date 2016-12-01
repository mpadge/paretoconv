#' Ramsay's (2006) chi function of Eq. (8)
#'
#' @param x value of independent variable
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

#' Calculates integral for CDF or PDF function, f
#'
#' @param f The integrand function to be integrated
#' @param x value of independent variable
#' @param a The primary shape parameter of the Pareto distribution 
#' @param n Number of convolutions
#' @param incr Proportional factor by which to increase upper limit at each step
#' @param rough If TRUE, calculates a rough upper limit from the integrand
#' rather than the integral, taking the first point where both absolute values
#' and differences decrease below tol.
#'
#' @note This is used for both non-integer and integer versions.
#'
#' @note This procedure generally involves very few loops, so \code{rough}
#' doesn't really make much difference to speed, yet yields inferior estimates
#' of the integral.
#'
#' @return Single value of integral
calc_integral <- function (f, x, a, n, incr=0.1, rough=FALSE)
{
    # Setting appropriate upper limits is very important because integrals in
    # these non-integer cases often diverge for high upper limits. An
    # appropriate upper limit is first found here by finding two consecutive
    # values of the integrand (f) for which the latter lies within tol of zero
    # and for which the difference between them is also < tol.
    #
    # Starting value is taken as 10 * x, increasing by incr each step.
    tol <- 1e-5
    chng <- 1
    upper <- 10 * x
    if (rough)
        val_old <- do.call (f, list (z=upper, x=x, a=a, n=n))
    else
        val_old <- cubature::adaptIntegrate (f, lowerLimit=0, upperLimit=upper,
                                             x=x, a=a, n=n)$integral
    maxiter <- 1e3
    niter <- 1
    while (chng > tol && niter < maxiter)
    {
        upper <- upper * (1 + incr)
        if (rough)
            val <- do.call (f, list (z=upper, x=x, a=a, n=n))
        else
            val <- cubature::adaptIntegrate (f, lowerLimit=0, upperLimit=upper,
                                             x=x, a=a, n=n)$integral
        if (abs (val) < tol)
            chng <- 0 # will stop loop at first instance, so will never / 0
        else
            chng <- abs (val - val_old) / val_old
        val_old <- val
        niter <- niter + 1
    }
    if (niter >= maxiter)
        warning ("Integral did not converge!")
    if (rough)
        val <- cubature::adaptIntegrate (f, lowerLimit=0, upperLimit=upper, 
                                         x=x, a=a, n=n)$integral
    return (val)
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

    # In this integrand, z is ramsay's x and x is his t, so his F_n(t) is here
    # F_n(x), and the integral is over z-values
    integrand <- function (z, x, a, n) 
        (1 - exp (-x * z / bet)) * Re (chi (z, a, n)) / z
    calc_integral (integrand, x, a, n)
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

    # In this integrand, z is ramsay's x and x is his t, so his F_n(t) is here
    # F_n(x), and the integral is over z-values
    integrand <- function (z, x, a, n)
        exp (-x * z / bet) * Re (chi (z, a, n))
    calc_integral (integrand, x, a, n)
}
