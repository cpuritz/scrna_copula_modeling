###############################################################################

#' Fit a t copula to a SingleCellExperiment
#'
#' @param sce A SingleCellExperiment.
#' @param margins Type of marginal distribution. Either "empirical" or "par". If
#' the former, margins are modeled empirically. If the latter, margins are
#' modeled as negative binomial/zero-inflated negative binomial (selected by
#' AIC).
#' @param Sigma Optional precomputed estimate of scale matrix.
#'
#' @returns A \code{tcop} object.
#'
#' @export
fitT <- function(sce,
                 margins = c("empirical", "par"),
                 Sigma = NULL) {
    if (!inherits(sce, "SingleCellExperiment")) {
        stop("'sce' must be a SingleCellExperiment.")
    }
    ngene <- dim(sce)[1]
    margins <- match.arg(margins)

    # Count matrix
    X <- as.matrix(Matrix::t(SingleCellExperiment::counts(sce)))

    # Construct marginal distribution functions
    mfun <- ifelse(margins == "empirical", empcdf, par_margin)
    margins <- apply(X, 2, mfun, simplify = FALSE)

    # Evaluate margins at samples
    FX <- sapply(1:ngene, function(i) { margins[[i]](X[, i]) })
    # Evaluate margins at left limits of samples
    FXm <- sapply(1:ngene, function(i) { margins[[i]](X[, i] - 1) })

    # Push values away from boundaries of unit cube
    eps <- 1e-12
    FX[FX > 1 - eps] <- 1 - eps
    FX[FX < eps] <- eps
    FXm[FXm > 1 - eps] <- 1 - eps
    FXm[FXm < eps] <- eps

    # Fit t copula via maximum likelihood
    res <- t_mle(
        FX = FX,
        FXm = FXm,
        Sigma = Sigma
    )

    # Construct t copula
    cop <- copula::tCopula(
        param = copula::P2p(res$Sigma),
        df = res$df,
        dim = ngene,
        dispstr = "un"
    )

    # Construct S3 object
    obj <- list(copula = cop,
                nobs = dim(X)[1],
                dim = dim(X)[2],
                npars = copula::nParam(cop, freeOnly = TRUE),
                names = colnames(X))
    class(obj) <- "tcop"
    return(obj)
}

###############################################################################

#' Maximum Likelihood Estimation of t Copula Parameters
#'
#' Maximum likelihood estimation of the scale matrix and degrees of freedom for
#' a t copula model of count data.
#'
#' @param FX Margins evaluated at samples.
#' @param FXm Margins evaluated at left limits of samples.
#' @param Sigma Optional precomputed estimate of scale matrix.
#' @param ... Additional arguments to pass to \code{stats::optim}.
#'
#' @returns Maximum likelihood estimation of scale matrix and degrees of
#' freedom.
t_mle <- function(FX,
                  FXm,
                  Sigma = NULL,
                  ...) {
    n <- dim(FX)[1]
    d <- dim(FX)[2]

    # If needed, estimate scale matrix under Gaussian assumption
    if (is.null(Sigma)) {
        Sigma <- gaussian_mle(FX, FXm, ...)
    }

    # Compute log-likelihood for fixed df
    loglik <- function(nu) {
        # Transform margins to standard t with df = nu
        TXm <- stats::qt(FXm, df = nu)
        TX <- stats::qt(FX, df = nu)

        # Compute average log-likelihood
        ll <- tryCatch(
            mean(log(sapply(1:n, function(i) {
                TruncatedNormal::mvTcdf(
                    l = TXm[i, ],
                    u = TX[i, ],
                    Sig = Sigma,
                    df = nu,
                    n = 1e3
                )$prob
            }))),
            error = function(e) {
                # Error often thrown for poorly conditioned matrices
                -Inf
            }
        )

        if (is.infinite(ll)) {
            ll <- sign(ll) * 1e20
        }
        # Return negative log-likelihood
        return(-ll)
    }

    # Same initial value as in copula package
    nu0 <- 4
    mle <- stats::optim(
        par = nu0,
        fn = loglik,
        lower = 2,      # lower bound at 2 to avoid undefined covariance
        upper = 300,    # large but finite upper bound
        method = "L-BFGS-B",
        ...
    )
    return(list(Sigma = Sigma, df = mle$value))
}

###############################################################################

#' Print method for S3 class tcop
#'
#' @param x Object of class \code{tcop}.
#' @param ... Additional arguments.
#'
#' @export
#' @method print tcop
print.tcop <- function(x, ...) {
    l1 <- paste0(x$dim, "-dimensional t copula fit")
    l2 <- paste(paste("nobs =", x$nobs),
                paste("npars =", x$npars),
                sep = "   ")
    cat(paste(l1, l2, sep = "\n"))
}

###############################################################################
