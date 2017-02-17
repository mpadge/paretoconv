#' Ramsay's (2006) modified phi function of Eq. (16)
#'
#' @param z value of independent variable
#' @param a The primary shape parameter of the Pareto distribution - alpha in
#' Ramsay's notation. 
#' @param n Number of convolutions
#'
#' @return Single value of phi
phi16 <- function (z, a, n)
{
    gam <- -digamma (1) # Euler's constant

    m <- a + 1 
    # Eq (13): phi_{m,n}(\nu) \sim (Ei_{m+1}(\nu)); here m=a+1
    if (m < 2)
        t1 <- gam + log (z)
    else
        t1 <- gam + log (z) - sum(1/1:(m - 1))
    t1 <- exp (-z) * z^(m - 1) * t1 / factorial (m-1)
    # Note that R's factorial goes up to and includes 170
    qq <- 0:170
    qq <- qq [!qq %in% (m-1)]
    # TODO: Compare timing using cumprod instead of factorial
    # Then note that z^qq generates Inf for large z, but in which case
    # exp(-z)/q! will always be smaller, so Inf values may safely be
    # presumed to be zero 
    t2 <- z ^ qq * exp (-z) / ((qq - m + 1) * factorial (qq))
    t2 <- sum (t2 [which (is.finite (t2) & !is.nan (t2))])
    eim <- t1 + t2 

    # Then the sum for phi_{a,n}(nu)
    rr <- 0:floor((n-1)/2)
    summand <- (-pi^2)^rr * choose (n, 2*rr+1) * eim ^ (n - 2 * rr - 1) * 
        (z ^ a * exp (-z) / factorial (a)) ^ (2 * rr + 1)
    (-1) ^ (n + 1) * a ^ n * sum (summand)
}

#' CDF for Convolution of Pareto distributions for integer alpha
#'
#' Calculates complementary Cumulative Distribution Function (cdf) from
#' convolution of multiple Pareto distributions following 'The Distribution of
#' Sums of Certain I.I.D. Pareto Variates' by Colin Ramsay (Communications in
#' Statistics - Theory and Methods 35:395-405, 2006).
#'
#' This is a combination of Eqs. (9), (13), and (14), with (13) actually
#' evaluated in the form of (16).
#'
#' @param x value of independent variable
#' @param a The primary shape parameter of the Pareto distribution - alpha in
#' Ramsay's notation. 
#' @param n Number of convolutions
#' @param x0 Lower cut-off point of classical heavy-tailed distribution
#' (generally obtained emprically with the poweRlaw package).
#'
#' @note This version is for integer values of a. Non-Integer values can be
#' calculated with ramsay_nonint
#'
#' @return Value for the CDF of the convolution of two Pareto distributions of
#' shape a at the value x.
ramsay_int_cdf <- function (x, a, n, x0)
{
    integrand <- function (z, x, a, n) 
    {
        if (z != 0)
            ret <- (1 - exp (-x * z / x0)) * phi16 (z, a, n) / z
        else
            ret <- 0
        return (ret)
    }
    if (x == 0)
        1
    else
        1 - calc_integral (integrand, x, a, n, incr=10, rough=TRUE)
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
#' @param x0 Lower cut-off point of classical heavy-tailed distribution
#' (generally obtained emprically with the poweRlaw package).
#'
#' @note This version is for integer values of a. Non-Integer values can be
#' calculated with ramsay_nonint
#'
#' @return Value for the PDF of the convolution of two Pareto distributions of
#' shape a at the value x.
ramsay_int_pdf <- function (x, a, n, x0)
{
    integrand <- function (z, x, a, n) 
    {
        if (z == 0)
            ret <- 0
        else
            ret <- exp (-x * z / (n * x0)) * phi16 (z / n, a, n)
        return (ret)
    }
    if (x == 0)
        0
    else
        calc_integral (integrand, x, a, n, incr=10, rough=TRUE) / (n * x0)
}
