#' Ramsay's (2006) modified phi function of Eq. (16)
#'
#' @param x value of independent variable
#' @param a The primary shape parameter of the Pareto distribution - alpha in
#' Ramsay's notation. 
#' @param n Number of convolutions
#'
#' @return Single value of phi
phi16 <- function (x, a, n)
{
    gam <- -digamma (1) # Euler's constant

    m <- a + 1 
    # Eq (13): phi_{m,n}(\nu) \sim (Ei_{m+1}(\nu)); here m=a+1
    if (m < 2)
        t1 <- gam + log (x)
    else
        t1 <- gam + log (x) - sum(1/1:(m - 1))
    t1 <- exp (-x) * x^(m - 1) * t1 / factorial (m-1)
    # Note that R's factorial goes up to and includes 170
    qq <- 0:170
    qq <- qq [!qq %in% (m-1)]
    # TODO: Compare timing using cumprod instead of factorial
    # Then note that x^qq generates Inf for large x, but in which case
    # exp(-x)/q! will always be smaller, so Inf values may safely be
    # presumed to be zero 
    t2 <- x ^ qq * exp (-x) / ((qq - m + 1) * factorial (qq))
    t2 <- sum (t2 [which (is.finite (t2) & !is.nan (t2))])
    eim <- t1 + t2 

    # Then the sum for phi_{a,n}(nu)
    rr <- 0:floor((n-1)/2)
    summand <- (-pi^2)^rr * choose (n, 2*rr+1) * eim ^ (n - 2 * rr - 1) * 
        (x ^ a * exp (-x) / factorial (a)) ^ (2 * rr + 1)
    (-1) ^ (n + 1) * a ^ n * sum (summand)
}

#' CDF for Convolution of Pareto distributions for integer alpha
#'
#' Calculates Cumulative Distribution Function (cdf) from convolution of
#' multiple Pareto distributions following 'The Distribution of Sums of Certain
#' I.I.D. Pareto Variates' by Colin Ramsay (Communications in Statistics -
#' Theory and Methods 35:395-405, 2006).
#'
#' This is a combination of Eqs. (9), (13), and (14), with (13) actually
#' evaluated in the form of (16).
#'
#' @param x value of independent variable
#' @param a The primary shape parameter of the Pareto distribution - alpha in
#' Ramsay's notation. 
#' @param n Number of convolutions
#'
#' @note This version is for integer values of a. Non-Integer values can be
#' calculated with ramsay_nonint
#'
#' @return Value for the CDF of the convolution of two Pareto distributions of
#' shape a at the value x.
ramsay_int_cdf <- function (x, a, n)
{
    bet <- set_beta () # Always = 1

    integrand <- function (y, x, a, n) 
    {
        if (y != 0)
            ret <- (1 - exp (-x * y / bet)) * phi16 (y, a, n) / y
        else
            ret <- 0
        return (ret)
    }
    upper <- 1e2
    count <- 0
    while (integrand (upper, max (x), a, n) != 0)
    {
        upper <- upper * 10
        count <- count + 1
        if (count > 6)
            stop ('Integrand not convergent')
    }
    cubature::adaptIntegrate (integrand, lowerLimit=0, upperLimit=upper,
                              x=x, a=a, n=n)$integral
}

#' PDF for Convolution of Pareto distributions for integer alpha
#'
#' Calculates Probability Density Function (pdf) from convolution of
#' multiple Pareto distributions following 'The Distribution of Sums of Certain
#' I.I.D. Pareto Variates' by Colin Ramsay (Communications in Statistics -
#' Theory and Methods 35:395-405, 2006).
#'
#' This is a combination of Eqs. (9), (13), and (12), with (13) actually
#' evaluated in the form of (16).
#'
#' @param x value of independent variable
#' @param a The primary shape parameter of the Pareto distribution - alpha in
#' Ramsay's notation. 
#' @param n Number of convolutions
#'
#' @note This version is for integer values of a. Non-Integer values can be
#' calculated with ramsay_nonint
#'
#' @return Value for the PDF of the convolution of two Pareto distributions of
#' shape a at the value x.
ramsay_int_pdf <- function (x, a, n)
{
    bet <- set_beta () # Always = 1

    integrand <- function (y, x, a, n) 
    {
        if (y == 0)
            ret <- 0
        else
            ret <- exp (-x * y / (n * bet)) * phi16 (y / n, a, n)
        return (ret)
    }
    upper <- 1e2
    count <- 0
    while (integrand (upper, max (x), a, n) != 0)
    {
        upper <- upper * 10
        count <- count + 1
        if (count > 6)
            stop ('Integrand not convergent')
    }
    # Integrals sometimes don't seem to converge even through the error is well
    # below tol. A work-around is to set maxEval, and only remove it if error is
    # > tol
    int <- cubature::adaptIntegrate (integrand, lowerLimit=0, upperLimit=upper,
                                     x=x, a=a, n=n, maxEval=1e3)
    if (int$error > 1e-5) # default adaptIntegrate tolerance
        int <- cubature::adaptIntegrate (integrand, lowerLimit=0,
                                         upperLimit=upper, x=x, a=a, n=n)

    int$integral / (n * bet)
}
