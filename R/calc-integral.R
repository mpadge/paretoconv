#' Calculates integral for CDF or PDF function, f
#'
#' @param f The integrand function to be integrated
#' @param x value of independent variable
#' @param a The primary shape parameter of the Pareto distribution 
#' @param x0 Lower cut-off point of classical heavy-tailed distribution
#' @param n Number of convolutions
#' @param incr Proportional factor by which to increase upper limit at each step
#'
#' @note This function calculates the integrals defined by Ramsay and used to
#' evaoluate the PDF/CDF functions (for both non-integer and integer versions).
#' The integrals of the PDF are calculated with the separate function,
#' \code{pdf_integral}.
#'
#' @return Single value of integral
#' @noRd
calc_integral <- function (f, x, a, x0, n, incr=0.1)
{
    # This parameter is for hcubature and is important for very low values of x
    # for which hcubature often fails to converge at all, presumably because the
    # tolerance is relative and low x -> *VERY* low integrals, and so unreliable
    # measures of relative tolerance. The resultant values are very, very low
    # regardless, so stopping at maxEval is fine in these cases, and 1e4 seems
    # to generate perfectly fine values yet not affect "normal" operation at
    # all.
    maxEval <- 0 # default of no limit
    if (x < x0 / 4)
        maxEval <- 1e4

    # Setting appropriate upper limits is very important because integrals in
    # these non-integer cases often diverge for high upper limits. An
    # appropriate upper limit is first found here by finding two consecutive
    # values of the integrand (f) which are both within tol of zero
    tol <- 1e-10
    val <- chng <- upper <- 1
    nlim <- 1e6
    count <- 1
    err <- FALSE
    while (chng > tol)
    {
        upper <- upper * (1 + incr)
        val_old <- val
        val <- do.call (f, list (z=upper, x=x, a=a, x0=x0, n=n))
        chng <- abs (val) + abs (val_old) 
        count <- count + 1
        if (count > nlim)
        {
            err <- TRUE
            break
        }
    }

    if (err)
    {
        warning ("Integral did not converge")
        ret <- NA
    } else
    {
        upper <- ceiling (upper)
        ret <- cubature::hcubature (f, lowerLimit=0, upperLimit=upper,
                             x=x, a=a, x0=x0, n=n, maxEval=1e4)$integral
    }
    return (ret)
}
