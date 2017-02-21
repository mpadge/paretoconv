#' Kolmogorow-Smirnov distance between convoluted Pareto distribution and
#' observed data
#'
#' @param m A \code{poweRlaw::displ} object containing data to be modelled
#' @param x0 Initial guess for lower limit of Pareto distribution
#' @param n Initial guess for number of convolutions
#' @param quiet If FALSE, progress information is dumped to screen
#'
#' @return A single value of the KS statistic
ks_dist <- function (m, x0, n, quiet=TRUE)
{
    if (!is (m, 'displ')) stop ('m must be a displ object')

    xvals <- sort (unique (m$getDat ()))
    y0 <- as.vector (table (sort (m$getDat ())))
    cdf <- rev (cumsum (rev (y0)))
    cdf <- cdf / max (cdf)
    y <- paretoconv (xvals, a=m$getPars (), x0=x0, n=n, cdf=TRUE, quiet=quiet)
    yfit <- y * lm (y ~ cdf + 0)$coefficients
    max (abs (yfit - cdf))
}

#' Locally optimal fitting of convoluted Pareto distribution to observed data
#'
#' @param m A \code{poweRlaw::displ} object containing data to be modelled
#' @param x0 Initial guess for lower limit of Pareto distribution
#' @param n Initial guess for number of convolutions
#' @param quiet If FALSE, progress information is dumped to screen
#'
#' @note This function only finds local optima - it is up to the user to ensure
#' that the given start values are near the global optimum. Values of \code{n=0}
#' are NOT searched.
#'
#' @return Position of the local optimum as quantified by \code{x0} and
#' \code{n}, along with associated Kolmogorow-Smirnov statistic quantifying
#' maximal distance from convoluted Pareto Cumulative Distribution Function
#' (CDF) and empirical CDF of model \code{m}.
#'
#' @export
pareto_optimise <- function (m, x0=1, n=1, quiet=TRUE)
{
    if (!is (m, 'displ'))
        stop ('m must be a poweRlaw::displ object')

    x0vec <- rep ((x0 - 1):(x0 + 1), 3)
    nvec <- rep ((n - 1):(n + 1), each=3)
    indx <- which (x0vec > 0 & nvec > 0)
    x0vec <- x0vec [indx]
    nvec <- nvec [indx]
    if (!quiet)
        message ('calculating initial KS statistics ...')
    y <- sapply (seq (indx), function (i) ks_dist (m, x0=x0vec [i], n=nvec [i]))
    i <- which.min (y)
    x0vec_old <- x0vec
    nvec_old <- nvec
    count <- 1
    while (which.min (y) != 5) # while min is not mid-point
    {
        x0 <- x0vec [which.min (y)]
        n <- nvec [which.min (y)]
        if (!quiet)
            message ('iteration#', count, ' -> (', x0, ', ', n, ')')
        x0vec_old <- x0vec
        nvec_old <- nvec
        x0vec <- rep ((x0 - 1):(x0 + 1), 3)
        nvec <- rep ((n - 1):(n + 1), each=3)
        indx <- which (x0vec %in% x0vec_old & nvec %in% nvec_old)
        indx_old <- which (x0vec_old %in% x0vec & nvec_old %in% nvec)
        y_old <- y
        y <- rep (NA, 9)
        y [indx] <- y_old [indx_old]
        y [which (is.na (y))] <- sapply (which (is.na (y)), function (i)
                                         ks_dist (m, x0=x0vec [i], n=nvec [i]))
        count <- count + 2
    }
    ret <- c (x0, n, min (y))
    names (ret) <- c ('x0', 'n', 'KS')
    return (ret)
}
