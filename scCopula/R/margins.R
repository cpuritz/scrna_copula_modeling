###############################################################################

#' Empirical Cumulative Distribution Function
#'
#' @description This is simply stats::ecdf scaled by n/(n+1) to avoid
#' the ECDF evaluating to 1.
#'
#' @param x Numeric vector.
#'
#' @returns A function implementing the ECDF of x.
#'
#' @export
empcdf <- function(x) {
    x <- sort(x)
    n <- length(x)
    if (n < 1)  {
        stop("'x' must have 1 or more non-missing values.")
    }
    vals <- unique(x)
    y <- cumsum(tabulate(match(x, vals))) / (n + 1)
    fun <- stats::approxfun(
        x = vals,
        y = y,
        method = "constant",
        yleft = 0,
        yright = y[length(y)],
        f = 0,
        ties = "ordered"
    )
    return(fun)
}

###############################################################################
