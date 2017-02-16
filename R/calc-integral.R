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
#' @note This function calculates the integrals defined by Ramsay and used to
#' evaoluate the PDF/CDF functions (for both non-integer and integer versions).
#' The integrals of the PDF are calculated with the separate function,
#' \cdoe{pdf_integral}.
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
