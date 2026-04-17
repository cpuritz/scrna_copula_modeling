###############################################################################

#' Fit a Gaussian copula to a SingleCellExperiment
#'
#' @param sce A SingleCellExperiment.
#' @param margins Type of marginal distribution. Either "empirical" or "par". If
#' the former, margins are modeled empirically. If the latter, margins are
#' modeled as negative binomial/zero-inflated negative binomial (selected by
#' AIC).
#' @param estimate How should the correlation matrix be estimated. One of
#' \code{"sample"}, \code{"jitter"}, or \code{"mle"}. The first uses the
#' sample correlation matrix of the pseudo-observations. The second uses the
#' sample correlation matrix of jittered pseudo-observations. The third performs
#' maximum likelihood estimation using the count data.
#' @param cores Number of cores.
#'
#' @returns A \code{gaussiancop} object.
#'
#' @export
fitGaussian <- function(sce,
                        margins = c("empirical", "par"),
                        estimate = c("sample", "jitter", "mle"),
                        cores = 1L) {
    if (!inherits(sce, "SingleCellExperiment")) {
        stop("'sce' must be a SingleCellExperiment.")
    }
    ngene <- dim(sce)[1]
    margins <- match.arg(margins)
    estimate <- match.arg(estimate)
    cores <- as.integer(cores)

    # Count matrix
    X <- as.matrix(Matrix::t(SingleCellExperiment::counts(sce)))

    # Construct marginal distribution functions
    mfun <- ifelse(margins == "empirical", empcdf, par_margin)
    margins <- apply(X, 2, mfun, simplify = FALSE)

    # Evaluate margins at samples
    FX <- sapply(1:ngene, function(i) { margins[[i]](X[, i]) })
    # Evaluate margins at left limits of samples
    FXm <- sapply(1:ngene, function(i) { margins[[i]](X[, i] - 1) })

    eps <- 1e-12
    if (estimate == "sample") {
        # Push values away from boundaries of unit cube
        FX[FX > 1 - eps] <- 1 - eps
        FX[FX < eps] <- eps

        # Sample correlation matrix
        Sigma <- stats::cor(stats::qnorm(FX))
    } else if (estimate == "jitter") {
        # Sample correlation matrix of jittered data
        V <- matrix(stats::runif(prod(dim(FX))), nrow = nrow(FX))
        U <- FXm + (FX - FXm) * V

        # Push values away from boundaries of unit cube
        U[U > 1 - eps] <- 1 - eps
        U[U < eps] <- eps

        # Sample correlation matrix
        Sigma <- stats::cor(stats::qnorm(U))
    } else {
        # Push values away from boundaries of unit cube
        FX[FX > 1 - eps] <- 1 - eps
        FX[FX < eps] <- eps
        FXm[FXm > 1 - eps] <- 1 - eps
        FXm[FXm < eps] <- eps

        Sigma <- gaussian_mle(FX, FXm, cores = cores)
    }

    # Construct copula
    cop <- copula::normalCopula(
        param = copula::P2p(Sigma),
        dim = ngene,
        dispstr = "un"
    )

    # Construct S3 object
    obj <- list(copula = cop,
                nobs = dim(X)[1],
                dim = ngene,
                npars = copula::nParam(cop, freeOnly = TRUE),
                names = colnames(X))
    class(obj) <- "gaussiancop"
    return(obj)
}

###############################################################################

#' Maximum Likelihood Estimation of Gaussian Correlation Matrix for Count Data
#'
#' Maximum likelihood estimation of the correlation matrix for a Gaussian
#' copula model of count data. Optimization is parallelized.
#'
#' @param FX Margins evaluated at samples.
#' @param FXm Margins evaluated at left limits of samples.
#' @param cores Number of cores to use for parallel computation.
#'
#' @returns Maximum likelihood estimation of correlation matrix.
gaussian_mle <- function(FX,
                         FXm,
                         cores) {
    n <- dim(FX)[1]
    d <- dim(FX)[2]
    NX <- stats::qnorm(FX)
    NXm <- stats::qnorm(FXm)

    # Unconstrained vector parametrization of correlation matrices
    loglik <- function(par) {
        # Build correlation matrix from vector input
        R <- vec2cor(par)
        if (any(eigen(R)$values <= 0)) {
            # Although the matrix is PD, if the smallest eigenvalue is within
            # machine precision of zero, it can be numerically negative and
            # this can cause a segfault in mvNcdf.
            ll <- -Inf
        } else {
            # Compute average log-likelihood
            ll <- tryCatch(
                mean(log(sapply(seq(dim(NX)[1]), function(i) {
                    TruncatedNormal::mvNcdf(
                        l = NXm[i, ],
                        u = NX[i, ],
                        Sig = R,
                        n = 1e3
                    )$prob
                }))),
                error = function(e) {
                    # Error often thrown for poorly conditioned matrices
                    -Inf
                }
            )
        }

        if (is.infinite(ll)) {
            ll <- sign(ll) * 1e20
        }
        # Return negative log-likelihood
        return(-ll)
    }

    # Sample correlation matrix as initial condition
    par0 <- cor2vec(stats::cor(NX, method = "pearson"))

    # Unconstrained minimization of negative log-likelihood
    if (cores > 1L) {
        # Create cluster
        cl <- parallel::makeCluster(cores)
        on.exit(parallel::stopCluster(cl))
        parallel::clusterExport(
            cl = cl,
            varlist = c("vec2cor", "NX", "NXm"),
            envir = environment()
        )

        opt <- optimParallel::optimParallel(
            par = par0,
            fn = loglik,
            parallel = list(cl = cl, forward = TRUE, loginfo = FALSE)
        )
    } else {
        opt <- stats::optim(
            par = par0,
            fn = loglik,
            method = "L-BFGS-B"
        )
    }

    # Map back to correlation matrix
    return(vec2cor(opt$par))
}

###############################################################################

#' Print method for S3 class gaussiancop
#'
#' @param x Object of class \code{gaussiancop}.
#' @param ... Additional arguments.
#'
#' @export
#' @method print gaussiancop
print.gaussiancop <- function(x, ...) {
    l1 <- paste0(x$dim, "-dimensional Gaussian copula fit")
    l2 <- paste(paste("nobs =", x$nobs),
                paste("npars =", x$npars),
                sep = "   ")
    cat(paste(l1, l2, sep = "\n"))
}

###############################################################################
