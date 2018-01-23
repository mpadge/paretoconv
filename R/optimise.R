#' Kolmogorow-Smirnov distance between convoluted Pareto distribution and
#' observed data
#'
#' @param m A \code{poweRlaw::displ} object containing data to be modelled
#' @param x0 Initial guess for lower limit of Pareto distribution
#' @param n Initial guess for number of convolutions
#' @param quiet If FALSE, progress information is dumped to screen
#'
#' @return A single value of the KS statistic
#' @noRd
ks_dist <- function (m, x0, n, quiet=TRUE)
{
    if (!is (m, 'displ')) stop ('m must be a displ object')

    xvals <- sort (unique (m$getDat ()))
    y0 <- as.vector (table (sort (m$getDat ())))
    cdf <- rev (cumsum (rev (y0)))
    cdf <- cdf / max (cdf)
    y <- paretoconv (xvals, a=m$getPars (), x0=x0, n=n, cdf=TRUE, quiet=quiet)

    # Then the model y is shifted until the KS statistic is minimised
    const <- 10
    n <- 0
    while (abs (const - 1) > 1e-10)
    {
        i1 <- which.min (cdf - y)
        i2 <- which.max (cdf - y)
        const <- (cdf [i1] + cdf [i2]) / (y [i1] + y [i2])
        y <- y * const
        n <- n + 1
        if (n > 100) # this should never happen
            stop ('KS statistic calculation did not converge')
    }
    max (abs (y - cdf))
}


#' Locally optimal fitting of convoluted Pareto distribution to observed data
#'
#' @param m A \code{poweRlaw::displ} object containing data to be modelled
#' @param x0 Initial guess for lower limit of Pareto distribution
#' @param n Initial guess for number of convolutions
#' @param check_non_conv If TRUE first checks whether a non-convoluted model may
#' be optimal
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
pareto_optimise <- function (m, x0=1, n=1, check_non_conv=TRUE, quiet=TRUE)
{
    if (!is (m, 'displ'))
        stop ('m must be a poweRlaw::displ object')

    # First check whether n=0 is an optimum. This is a special case because x0
    # has no effect for n=0
    y_n0 <- ks_dist (m, x0=m$getXmin (), n=0)
    # TODO: Improve this check, which currently only uses x=1:4
    if (!quiet)
        message ('checking whether non-convoluted version is optimal',
                 appendLF=FALSE)
    if (check_non_conv)
        y_n1 <- sapply (1:4, function (i) ks_dist (m, x0=i, n=1))
    else
        y_n1 <- rep (NA, 4)
    if (all (y_n1 > y_n0) & all (diff (y_n1) > 0))
    {
        if (!quiet)
            message (': YES')
        ret <- c (x0, 0, y_n0)
    } else
    {
        if (!quiet)
            message (': NO')
        x0vec <- rep ((x0 - 1):(x0 + 1), 3)
        nvec <- rep ((n - 1):(n + 1), each=3)
        indx <- which (x0vec > 0 & nvec >= 0)
        x0vec <- x0vec [indx]
        nvec <- nvec [indx]
        if (!quiet)
            message ('calculating initial KS statistics ...')
        y <- rep (NA, length (indx))
        if (check_non_conv)
        {
            indx <- which (x0vec %in% 1:4 & nvec == 1)
            y [indx] <- y_n1 [which (1:4 %in% x0vec)]
        }
        y [is.na (y)] <- sapply (which (is.na (y)), function (i) 
                                 ks_dist (m, x0=x0vec [i], n=nvec [i]))
        count <- 1
        at_min <- FALSE
        if (x0vec [which.min (y)] == x0 & nvec [which.min (y)] == n)
            at_min <- TRUE
        while (!at_min)
        {
            x0 <- x0vec [which.min (y)]
            n <- nvec [which.min (y)]
            if (!quiet)
                message ('iteration#', count, ' -> (x0, n) = (', x0, ', ', n, ')')
            x0vec_old <- x0vec
            nvec_old <- nvec
            x0vec <- rep ((x0 - 1):(x0 + 1), 3)
            nvec <- rep ((n - 1):(n + 1), each=3)
            indx <- which (x0vec > 0 & nvec >= 0)
            x0vec <- x0vec [indx]
            nvec <- nvec [indx]
            indx <- which (x0vec %in% x0vec_old & nvec %in% nvec_old)
            indx_old <- which (x0vec_old %in% x0vec & nvec_old %in% nvec)
            y_old <- y
            y <- rep (NA, 9)
            y [indx] <- y_old [indx_old]
            indx <- which (nvec >= 0)
            nvec <- nvec [indx]
            x0vec <- x0vec [indx]
            y <- y [indx]
            y [which (is.na (y))] <- sapply (which (is.na (y)), function (i)
                                             ks_dist (m, x0=x0vec [i], n=nvec [i]))
            count <- count + 1
            if (x0vec [which.min (y)] == x0 & nvec [which.min (y)] == n)
                at_min <- TRUE
            if (count > 1e4)
                stop ('pareto_optimise failed to converge')
        }
        ret <- c (x0, n, min (y))
    }
    names (ret) <- c ('x0', 'n', 'KS')
    return (ret)
}

#' Generate one simulated paretoconv model
#'
#' @param m A \code{poweRlaw::displ} object containing data to be modelled
#' @param x0 Initial guess for lower limit of Pareto distribution
#' @param n Initial guess for number of convolutions
#' @param check_non_conv If TRUE first checks whether a non-convoluted model may
#' be optimal
#' @param quiet If FALSE, display progress messages on screen
#'
#' @return Position of the local optimum as quantified by \code{x0} and
#' \code{n}, along with associated Kolmogorow-Smirnov statistic quantifying
#' maximal distance from convoluted Pareto Cumulative Distribution Function
#' (CDF) and empirical CDF of model \code{m}.
#' @noRd
sim_mod1 <- function (m, x0, n, check_non_conv=TRUE, quiet=TRUE)
{
    if (!is (m, 'displ'))
        stop ('m must be a poweRlaw::displ object')

    if (missing (x0) | missing (n))
    {
        dat <- pareto_optimise (m, x0=1, n=1, check_non_conv=check_non_conv,
                                quiet=quiet)
        x0 <- dat [1]
        n <- dat [2]
    }
    a <- m$getPars ()
    nx <- length (m$getDat ())
    xt <- rep (0, nx)
    for (i in 0:n)
        xt <- xt + rplconv (n=nx, x0=x0, alpha=a)
    xt <- floor (xt / (n + 1))
    m2 <- poweRlaw::displ$new (xt)
    m2$setXmin (poweRlaw::estimate_xmin (m2))
    pareto_optimise (m2, x0=x0, n=n, quiet=quiet)
}

#' Generate a series of simulated paretoconv models
#'
#' @param m A \code{poweRlaw::displ} object containing data to be modelled
#' @param x0 Initial guess for lower limit of Pareto distribution
#' @param n Initial guess for number of convolutions
#' @param times Minimum number of times most successful model should be
#' generated
#' @param quiet If FALSE, display progress messages on screen
#'
#' @return Series of models specified by\code{x0} and \code{n}, along with
#' numbers of times each of those models represented the optimal model
#' @noRd
sim_mod_series <- function (m, x0, n, times=4, quiet=TRUE)
{
    mods <- counts <- NULL
    enough <- FALSE
    if (!quiet)
        message ('Generating models until same model appears ',
                 times, ' times')
    while (!enough)
    {
        mod <- sim_mod1 (m, x0=x0, n=n, check_non_conv=FALSE, quiet=quiet)
        i <- which (mod [1] == mods [,1] & mod [2] == mods [,2])
        if (length (i) == 0)
        {
            mods <- rbind (mods, mod [1:2])
            counts <- c (counts, 1)
        } else
            counts [i] <- counts [i] + 1
        if (max (counts) >= times)
            enough <- TRUE
    }
    cbind (mods, counts)
}

#' Estimate Kolmogorov-Smirnov (KS) statistic for a single synthetic paretoconv
#' model
#'
#' @param m A \code{poweRlaw::displ} object containing data to be modelled
#' @param x0 Initial guess for lower limit of Pareto distribution
#' @param n Initial guess for number of convolutions
#' @param x0mod Equivalent value of \code{x0} for synthetic model
#' @param nmod Equivalent value of \code{n} for synthetic model
#'
#' @return Single value of KS statistic
#' @noRd
ks1 <- function (m, x0, n, x0mod, nmod)
{
    a <- m$getPars ()
    nx <- length (m$getDat ())
    xt <- rep (0, nx)
    for (i in 0:n)
        xt <- xt + rplconv (n=nx, x0=x0, alpha=a)
    xt <- floor (xt / (n + 1))
    # Then the code from `ks_dist()`, include re-estimating a from simulated
    # data
    m2 <- poweRlaw::displ$new (xt)
    m2$setXmin (poweRlaw::estimate_xmin (m2))
    ks_dist (m2, x0=x0mod, n=nmod)
}

#' Estimate probability of observed KS statistic from synthetic models
#'
#' @param m A \code{poweRlaw::displ} object containing data to be modelled
#' @param x0 Initial guess for lower limit of Pareto distribution
#' @param n Initial guess for number of convolutions
#' @param neach Number of test statistics to be generated for each different
#' model
#' @param quiet If FALSE, display progress messages on screen
#'
#' @return Single value representing probability of observing given value of
#' \code{ks}
#'
#' @export
pparetoconv <- function (m, x0, n, neach=10, quiet=TRUE)
{
    ks0 <- ks_dist (m=m, x0=x0, n=n)
    if (!quiet)
        message ('Generating simulated models')
    mods <- sim_mod_series (m=m, x0=x0, n=n, quiet=quiet)
    ksvals <- NULL
    if (!quiet)
        message ('Generating synthetic series from ', nrow (mods),
                 ' simulated models')
    for (i in seq (nrow (mods)))
    {
        ni <- neach * mods [i, 3]
        if (!quiet)
            message ('mod [', mods [i, 1], ', ', mods [i, 2], '] ', 
                     appendLF=FALSE)
        for (j in seq (ni))
        {
            ksvals <- c (ksvals, 
                         ks1 (m, x0=x0, n=n, x0mod=mods [i, 1], nmod=mods [i, 2]))
            if (!quiet)
                message ('.', appendLF=FALSE)
        }
        if (!quiet)
            message ('')
    }
    if (!quiet)
        message ('Generated ', length (ksvals), ' synthetic KS statistics')
    dk <- density (ksvals, n=2^16)
    length (dk$y [dk$y > ks0]) / length (dk$y)
}
