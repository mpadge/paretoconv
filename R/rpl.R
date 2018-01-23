#' Generate random deviates for power laws
#'
#' @param n Number of deviates
#' @param x0 Traditional cutoff parameter of Pareto distribution
#' @param alpha Primary shape parameter of Pareto distribution
#' @param discrete_max Upper limit for generation (see \code{poweRlaw::rpldis}
#' for details)
#'
#' @note This is directly cribbed from \code{poweRlaw::rpldis}, modified only to
#' generate random deviates for all values including 0<x<xmin.
#'
#' @return Vector of \code{n} random deviates
#'
#' @export
rplconv <- function (n, x0, alpha, discrete_max=10000)
{
    x0 <- floor (x0)
    u <- runif (n)

    constant <- VGAM::zeta (alpha)
    if (x0 > 1)
        constant <- constant - sum ( (1:(x0 - 1)) ^ (-alpha))
    cdf <- c (0, 1 - (constant -
                      cumsum ( (x0:discrete_max) ^ (-alpha))) / constant)

    ## Due to numerical instability
    ## Not enough precision to exactly calculate the CDF
    dups <- duplicated (cdf, fromLast = TRUE)
    if (any (dups)) cdf <- cdf [1:which.min (!dups)]

    ## Simulate using look up method
    # Gillespie's code:
    #rngs <- as.numeric (cut (u, cdf)) + x0 - 1
    # but adding x0 only transposes the values; to start from 1, just
    rngs <- as.numeric (cut (u, cdf))

    ## Fill in blanks using Clauset approximation
    is_na <- is.na (rngs)
    if (any (is_na))
        rngs [is_na] <- floor ( (x0 - 0.5) *
                               (1 - u [is_na]) ^ (-1 / (alpha - 1)) + 0.5)
    rngs
}
