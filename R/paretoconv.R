#' Convolution of Pareto distributions 
#'
#' Convolutes multiple Pareto distributions following 
#' \itemize{
#' \item For integer shape parameters: 'The Distribution of Sums of Certain
#' I.I.D.  Pareto Variates' by Colin Ramsay (Communications in Statistics -
#' Theory and Methods 35:395-405, 2006);
#' \item For non-integer shape parameters: 'The Distribution of Sums of I.I.D.
#' Pareto Random Variables with Arbitrary Shape Parameter' by Colin Ramsay
#' (Communications in Statistics - Theory and Methods 37:2177-2184, 2008).
#' }
#'
#' @param x value of independent variable
#' @param a The primary shape parameter of the Pareto distribution - alpha in
#' Ramsay's notation.
#' @param n Number of convolutions
#'
#' @note The Pareto distribution may be defined as f(x)=(a/b)(b/x)^(a-1), where
#' a and b are the primary and secondary shape parameters, respectively. It
#' presumed here without loss of generality that b=1 and thus f(x)=a x^(1-a).
#'
#' @return Value for the convolution of two Pareto distributions of shape a at
#' the value x.
#'
#' @export
#' @examples 
#' paretoconv (1, 1, 1)
paretoconv <- function (x, a, n)
{
    if (!is.numeric (x)) stop ('x must be numeric')
    if (!is.numeric (a)) stop ('a must be numeric')
    if (!is.numeric (n)) stop ('n must be numeric')

    if (n%%1 != 0) stop ('n must be an integer')
    if (x < 0) stop ('x must be positive')

    if (a%%1 == 0) 
        ret <- ramsay_int (x, a, n)
    else
        ret <- ramsay_nonint (x, a, n)
    return (ret)
}
