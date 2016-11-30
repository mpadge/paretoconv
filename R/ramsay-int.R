#' Convolution of Pareto distributions for non-integer alpha
#'
#' Convolutes multiple Pareto distributions following 'The Distribution of Sums
#' of Certain I.I.D. Pareto Variates' by Colin Ramsay (Communications in
#' Statistics - Theory and Methods 35:395-405, 2006).
#'
#' @param x value of independent variable
#' @param a The primary shape parameter of the Pareto distribution - alpha in
#' Ramsay's notation. 
#' @param n Number of convolutions
#'
#' @note This version is for integer values of a. Non-Integer values can be
#' calculated with ramsay_nonint
#'
#' @return Value for the convolution of two Pareto distributions of shape a at
#' the value x.
#'
#' @export
#' @examples 
#' ramsay_int (1, 1, 1)
ramsay_int <- function (x, a, n)
{
    if (a%%1 != 0) stop ("a must be an integer")
    # euler's constant is -digamma(1)

    bet <- 1 # beta value fixed for my purposes always at 1
    gam <- -digamma (1)

    # Start with Ei_m of Eq.(9), but noting that this is used in Eq.(13), and
    # ultimately in the c.d.f. of Eq.(14), in which phi(nu) is multiplied by
    # exp(-n nu), so this multiple may be transferred into Eq.(9) to keep the
    # sums stable. The exponential term in Eq.(14) then becomes,
    # e^(-n nu) * e^nu = e^(-nu (n - 1)). The term in (13) is nevertheless
    # (Ei_(m+1) (nu)) ^ (n - 2r + 1), so actually requiring multiplication by
    # e^(-n nu) * (e^nu) ^ (n - 2r + 1) = e^(-n nu) * e^(n nu - 2r nu + nu)
    # = e^(-2r nu + nu) = e^(-nu (2r - 1)).
    #
    # The primary summation in Eq.(9) can then be calculated as
    # sum_{r=0}^\infty exp(-x) x ^ r / ((r - m + 1) * r!)

    phi <- function (x, a, n)
    {
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
    integrand <- function (y, x, a, n) 
    {
        if (y != 0)
            ret <- (1 - exp (-x * y / bet)) * phi (y, a, n) / y
        else
            ret <- 0
        return (ret)
    }
    cubature::adaptIntegrate (integrand, lowerLimit=0, upperLimit=1e3,
                              x=x, a=a, n=n)$integral
}
