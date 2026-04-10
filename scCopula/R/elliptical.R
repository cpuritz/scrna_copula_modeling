###############################################################################

#' Fit a Gaussian copula to a SingleCellExperiment
#'
#' @param sce A SingleCellExperiment.
#' @param margins One, or a list of, marginal distribution functions. If only
#' one function is passed, all margins are modeled identically. Entries can also
#' equal "empirical", in which case the corresponding margin is modeled
#' empirically.
#' @param mle Logical indicating whether maximum likelihood estimation should
#' be performed to estimate the correlation matrix. If \code{FALSE}, the sample
#' correlation matrix is used.
#' @param jitter Logical indicating whether data should be jittered.
#' @param likelihood Whether to compute likelihood.
#' @param cores Number of cores.
#'
#' @returns A \code{gaussiancop} object.
#'
#' @export
fitGaussian <- function(sce,
                        margins,
                        mle = TRUE,
                        jitter = FALSE,
                        likelihood = TRUE,
                        cores = 1L) {
    if (!inherits(sce, "SingleCellExperiment")) {
        stop("'sce' must be a SingleCellExperiment.")
    }
    ngene <- dim(sce)[1]
    if (length(margins) == 1) {
        margins <- rep(list(margins), ngene)
    }
    if (length(margins) != ngene) {
        stop("Length of 'margins' be one or the number of genes in 'sce'.")
    }
    if (!is.logical(mle)) {
        stop("'mle' must be 'logical'.")
    }
    cores <- as.integer(cores)

    # Count matrix
    X <- as.matrix(Matrix::t(SingleCellExperiment::counts(sce)))

    # Construct empirical distribution functions if necessary
    for (i in 1:ngene) {
        if (is.character(margins[[i]]) && margins[[i]] == "empirical") {
            margins[[i]] <- empcdf(X[, i])
        }
    }

    # Check that all margins are now functions
    if (!all(sapply(margins, is.function))) {
        stop("'margins' must be a list of functions.")
    }

    # Evaluate margins at samples
    FX <- sapply(1:ngene, function(i) { margins[[i]](X[, i]) })
    # Evaluate margins at left limits of samples
    FXm <- sapply(1:ngene, function(i) { margins[[i]](X[, i] - 1) })

    if (jitter) {
        # Sample correlation matrix of jittered data
        V <- matrix(stats::runif(prod(dim(FX))), nrow = nrow(FX))
        U <- FXm + (FX - FXm) * V
        Sigma <- stats::cor(stats::qnorm(U))

        # Construct copula
        cop <- copula::normalCopula(
            param = copula::P2p(Sigma),
            dim = ngene,
            dispstr = "un"
        )

        if (likelihood) {
            # Compute copula log-likelihood
            ll <- loglikCopulaCts(
                X = U,
                copula = cop,
                margins = margins
            )
        } else {
            ll <- NULL
        }
    } else {
        if (mle) {
            # Maximum likelihood estimation of correlation matrix
            if (cores > 1L) {
                Sigma <- gaussian_mle_par(FX, FXm, cores = cores)
            } else {
                Sigma <- gaussian_mle(FX, FXm)
            }
        } else {
            # Sample correlation matrix
            Sigma <- stats::cor(stats::qnorm(FX))
        }

        # Construct copula
        cop <- copula::normalCopula(
            param = copula::P2p(Sigma),
            dim = ngene,
            dispstr = "un"
        )
        if (likelihood) {
            # Compute copula log-likelihood
            ll <- loglikEllipCount(
                X = X,
                copula = cop,
                margins = margins
            )
        } else {
            ll <- NULL
        }
    }

    # Compute AIC and BIC
    n <- dim(X)[1]
    d <- dim(X)[2]
    k <- copula::nParam(cop, freeOnly = TRUE)
    if (likelihood) {
        aic <- 2 * k - 2 * ll
        bic <- log(n) * k - 2 * ll
    } else {
        aic <- NULL
        bic <- NULL
    }

    # Construct S3 object
    obj <- list(copula = cop,
                nobs = n,
                dim = d,
                npars = k,
                loglik = ll,
                aic = aic,
                bic = bic,
                names = colnames(X))
    class(obj) <- "gaussiancop"
    return(obj)
}

###############################################################################

#' Fit a t copula to a SingleCellExperiment
#'
#' @param sce A SingleCellExperiment.
#' @param margins One, or a list of, marginal distribution functions. If only
#' one function is passed, all margins are modeled identically. Entries can also
#' equal "empirical", in which case the corresponding margin is modeled
#' empirically.
#' @param Sigma Optional precomputed estimate of scale matrix.
#' @param likelihood Whether to compute likelihood.
#'
#' @returns A \code{tcop} object.
#'
#' @export
fitT <- function(sce,
                 margins,
                 Sigma = NULL,
                 likelihood = TRUE) {
    if (!inherits(sce, "SingleCellExperiment")) {
        stop("'sce' must be a SingleCellExperiment.")
    }
    ngene <- dim(sce)[1]
    if (length(margins) == 1) {
        margins <- rep(list(margins), ngene)
    }
    if (length(margins) != ngene) {
        stop("Length of 'margins' be one or the number of genes in 'sce'.")
    }

    # Count matrix
    X <- as.matrix(Matrix::t(SingleCellExperiment::counts(sce)))

    # Construct empirical distribution functions if necessary
    for (i in 1:ngene) {
        if (is.character(margins[[i]]) && margins[[i]] == "empirical") {
            margins[[i]] <- empcdf(X[, i])
        }
    }

    # Check that all margins are now functions
    if (!all(sapply(margins, is.function))) {
        stop("'margins' must be a list of functions.")
    }

    # Evaluate margins at samples
    FX <- sapply(1:ngene, function(i) { margins[[i]](X[, i]) })
    # Evaluate margins at left limits of samples
    FXm <- sapply(1:ngene, function(i) { margins[[i]](X[, i] - 1) })

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

    if (likelihood) {
        # Compute copula log-likelihood
        ll <- loglikEllipCount(
            X = X,
            copula = cop,
            margins = margins
        )
    } else {
        ll <- NULL
    }

    # Compute AIC and BIC
    n <- dim(X)[1]
    d <- dim(X)[2]
    k <- copula::nParam(cop, freeOnly = TRUE)
    if (likelihood) {
        aic <- 2 * k - 2 * ll
        bic <- log(n) * k - 2 * ll
    } else {
        aic <- NULL
        bic <- NULL
    }

    # Construct S3 object
    obj <- list(copula = cop,
                nobs = n,
                dim = d,
                npars = k,
                loglik = ll,
                aic = aic,
                bic = bic,
                names = colnames(X))
    class(obj) <- "tcop"
    return(obj)
}

###############################################################################

#' Elliptical Copula Likelihood for Count Data
#'
#' @description Compute log-likelihood for an elliptical copula model of count
#' data.
#'
#' @param X Matrix of data.
#' @param copula A \code{normalCopula} or \code{tCopula} object.
#' @param margins Vector of marginal distribution functions.
#' @param cores Number of cores to use for parallel computation.
#'
#' @returns A scalar value.
#'
#' @export
loglikEllipCount <- function(X,
                             copula,
                             margins,
                             cores = 1L) {
    if (is.data.frame(X)) {
        X <- as.matrix(X)
    }
    if (is.vector(X)) {
        X <- matrix(X, nrow = 1)
    }
    n <- dim(X)[1]
    d <- dim(X)[2]

    if (!inherits(copula, "normalCopula") && !inherits(copula, "tCopula")) {
        stop("'copula' must be a 'normalCopula' or 'tCopula' object.")
    }
    if (length(margins) != d) {
        stop("Length of 'margins' must equal number of columns in 'X'.")
    }
    cores <- as.integer(cores)

    # Evaluate margins at samples
    FX <- sapply(1:d, function(k) margins[[k]](X[, k]))
    # Evaluate margins at left limits of samples
    FXm <- sapply(1:d, function(k) margins[[k]](X[, k] - 1))

    # Numerical estimation of full density using QMC integration
    mu <- rep(0, d)
    Sigma <- copula::getSigma(copula)
    if (inherits(copula, "tCopula")) {
        # df must be an integer
        df <- copula@parameters[copula@param.names == "df"]
        if (df - floor(df) > 0) {
            message("Rounding 'df' to nearest integer.")
            df <- round(df)
        }
        # Convert margins to standard t
        TX <- stats::qt(FX, df = df)
        TXm <- stats::qt(FXm, df = df)
        # Compute likelihood at each sample
        P <- unlist(parallel::mclapply(1:n, function(i) {
            mvtnorm::pmvt(
                lower = TXm[i, ],
                upper = TX[i, ],
                corr = Sigma,
                df = df
            )
        }))
    } else {
        # Convert margins to standard normal
        NX <- stats::qnorm(FX)
        NXm <- stats::qnorm(FXm)
        # Compute likelihood at each sample
        P <- unlist(parallel::mclapply(1:n, function(i) {
            TruncatedNormal::mvNcdf(
                l = NXm[i, ],
                u = NX[i, ],
                Sig = Sigma,
                n = 1e3
            )$prob
        }))
    }

    # Compute marginal likelihoods
    # sum_{i=1..n} prod_{j=1..d} (F_j(X_ij) - F_j(X_ij - 1))
    M <- rep(1, length = n)
    for (i in 1:n) {
        for (j in 1:d) {
            Fj <- margins[[j]]
            M[i] <- M[i] * (Fj(X[i, j]) - Fj(X[i, j] - 1))
        }
    }

    # Compute copula log-likelihood. Copula density defined as
    # full density / marginal density.
    ll <- sum(log(P)) - sum(log(M))
    return(ll)
}

###############################################################################

#' Maximum Likelihood Estimation of Gaussian Correlation Matrix for Count Data
#'
#' Maximum likelihood estimation of the correlation matrix for a Gaussian
#' copula model of count data.
#'
#' @param FX Margins evaluated at samples.
#' @param FXm Margins evaluated at left limits of samples.
#' @param ... Additional arguments to pass to \code{stats::optim}.
#'
#' @returns Maximum likelihood estimation of correlation matrix.
gaussian_mle <- function(FX,
                         FXm,
                         ...) {
    n <- dim(FX)[1]
    d <- dim(FX)[2]
    NX <- stats::qnorm(FX)
    NXm <- stats::qnorm(FXm)

    # Unconstrained vector parametrization of correlation matrices
    loglik <- function(v) {
        # Build correlation matrix from vector input
        Sigma <- vec2cor(v)
        if (any(eigen(Sigma)$values <= 0)) {
            # Although the matrix is PD, if the smallest eigenvalue is within
            # machine precision of zero, it can be numerically negative and
            # this can cause a segfault in mvNcdf
            ll <- -Inf
        } else {
            # Compute average log-likelihood
            ll <- tryCatch(
                mean(log(sapply(1:n, function(i) {
                    TruncatedNormal::mvNcdf(
                        l = NXm[i, ],
                        u = NX[i, ],
                        Sig = Sigma,
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

    # Sample correlation matrix as starting point
    Sigma0 <- stats::cor(NX)
    # Parametrize in terms of partial correlations
    v0 <- cor2vec(Sigma0)

    # Unconstrained minimization of negative log-likelihood
    mle <- stats::optim(
        par = v0,
        fn = loglik,
        method = "L-BFGS-B",
        ...
    )

    # Map back to correlation matrix
    Sigma_mle <- vec2cor(mle$par)
    return(Sigma_mle)
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
gaussian_mle_par <- function(FX,
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
            # this can cause a segfault in mvNcdf
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

    # Create cluster
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(
        cl = cl,
        varlist = c("vec2cor", "NX", "NXm"),
        envir = environment()
    )

    # Sample correlation matrix as initial condition
    par0 <- cor2vec(stats::cor(NX, method = "pearson"))

    # Unconstrained minimization of negative log-likelihood
    opt <- optimParallel::optimParallel(
        par = par0,
        fn = loglik,
        parallel = list(cl = cl, forward = TRUE)
    )

    # Map back to correlation matrix
    return(vec2cor(opt$par))
}

###############################################################################

#' Parametrize correlation matrix as unconstrained vector
#'
#' @param R A d x d correlation matrix.
#'
#' @returns A vector of length choose(d, 2).
cor2vec <- function(R) {
    d <- dim(R)[1]

    # Compute Cholesky factor
    L <- t(chol(R))

    # Parametrize the Cholesky factor to have a real and positive diagonal
    H <- matrix(0, nrow = d, ncol = d)
    diag(H) <- 1
    H[, 1] <- L[, 1]
    for (i in 3:d) {
        H[i, 2:(i - 1)] <- L[i, 2:(i - 1)] / sqrt(1 - cumsum(L[i, 1:(i - 2)]^2))
    }

    # Map to an unconstrained vector in R^(d choose 2)
    vH <- H[lower.tri(H)]
    scale <- 1.2904
    return(scale * log(acos(vH) / (pi - acos(vH))))
}

###############################################################################

#' Inverse of cor2vec
#'
#' @param v A vector of length choose(d, 2).
#'
#' @returns A d x d correlation matrix.
vec2cor <- function(v) {
    d <- as.integer(round((1 + sqrt(1 + 8 * length(v))) / 2))
    scale <- 1.2904

    # Map unconstrained vector to d x d matrix with entries in (-1, 1)
    H <- matrix(0, nrow = d, ncol = d)
    diag(H) <- 1
    H[lower.tri(H)] <- cos(pi / (1 + exp(-v / scale)))

    # Map back to Cholesky factor space
    for (i in 2:d) {
        H[i, 2:i] <- H[i, 2:i] * sqrt(cumprod(1 - H[i, 1:(i - 1)]^2))
    }

    # Return correlation matrix
    return(tcrossprod(H, H))
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
                paste("logLik =", sprintf("%.2f", x$loglik)),
                paste("npars =", x$npars),
                paste("AIC =", sprintf("%.2f", x$aic)),
                paste("BIC =", sprintf("%.2f", x$bic)),
                sep = "   ")
    cat(paste(l1, l2, sep = "\n"))
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
                paste("logLik =", sprintf("%.2f", x$loglik)),
                paste("npars =", x$npars),
                paste("AIC =", sprintf("%.2f", x$aic)),
                paste("BIC =", sprintf("%.2f", x$bic)),
                sep = "   ")
    cat(paste(l1, l2, sep = "\n"))
}

###############################################################################
