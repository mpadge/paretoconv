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

    # colin ramsay 2008
}
