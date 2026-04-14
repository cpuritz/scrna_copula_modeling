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

#' Parametric CDF for count data
#'
#' @description Constructs either a negative binomial or zero-inflated negative
#' binomial CDF for count data. The family is automatically selected by
#' minimizing AIC.
#'
#' @param x Numeric vector of count data.
#'
#' @returns A function implementing the CDF of x.
#'
#' @export
par_margin <- function(x) {
    # Fit negative binomial
    fit_nb <- MASS::glm.nb(x ~ 1)

    # Fit zero-inflated negative binomial
    fit_zinb <- pscl::zeroinfl(x ~ 1 | 1, dist = "negbin")

    # Choose model with lower AIC
    aic_nb <- stats::AIC(fit_nb)
    aic_zinb <- stats::AIC(fit_zinb)

    if (aic_nb < aic_zinb) {
        # Extract NB parameters
        mu <- exp(fit_nb$coefficients[["(Intercept)"]])
        size <- fit_nb$theta

        # Build CDF
        pX <- function(q) {
            stats::pnbinom(q, size = size, mu = mu)
        }
    } else {
        # Extract ZINB parameters
        mu <- exp(fit_zinb$coefficients$count[["(Intercept)"]])
        size <- fit_zinb$theta
        pi <- stats::plogis(fit_zinb$coefficients$zero[["(Intercept)"]])

        # Build CDF
        pX <- function(q) {
            pi + (1 - pi) * stats::pnbinom(q, size = size, mu = mu)
        }
    }
    return(pX)
}

###############################################################################
