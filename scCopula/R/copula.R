###############################################################################

#' Copula Likelihood for Continuous Data
#'
#' @description Compute log-likelihood for a copula model of continuous data.
#'
#' @param X Matrix of data.
#' @param copula A \code{copula} object.
#' @param margins A list of marginal distribution functions. If \code{NULL},
#' samples are first converted to pseudo-observations.
#'
#' @returns A scalar value.
#'
#' @export
loglikCopulaCts <- function(X,
                            copula,
                            margins = NULL) {
    if (is.data.frame(X)) {
        X <- as.matrix(X)
    }
    if (is.vector(X)) {
        X <- matrix(X, nrow = 1)
    }
    d <- dim(X)[2]

    if (!inherits(copula, "copula")) {
        stop("'copula' must be a 'copula' object.")
    }
    if (!is.null(margins) && length(margins) != d) {
        stop("'margins' must have the same length as the number of columns in",
             " 'X'.")
    }

    # Convert margins to Unif(0,1)
    if (is.null(margins)) {
        # Compute pseudo-observations
        U <- copula::pobs(X)
    } else {
        # Use specified marginal distribution functions
        U <- sapply(1:d, function(k) { margins[[k]](X[, k]) })
    }

    # Compute log-likelihood
    ll <- copula::loglikCopula(u = U, copula = copula)

    return(ll)
}

###############################################################################

#' Copula Likelihood for Count Data
#'
#' @description Compute log-likelihood for a copula model of count data.
#'
#' @details This function takes a finite-difference approach and requires 2^d
#' evaluations of the copula distribution function.
#'
#' @param X Matrix of data.
#' @param copula A \code{copula} object.
#' @param margins Vector of marginal distribution functions.
#' @param cores Number of cores to use for parallel computation.
#'
#' @returns A scalar value.
#'
#' @export
loglikCopulaCount <- function(X,
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

    if (!inherits(copula, "copula")) {
        stop("'copula' must be a 'copula' object.")
    }
    if (length(margins) != d) {
        stop("'margins' must have the same length as the number of columns in",
             " 'X'.")
    }
    cores <- as.integer(cores)

    # All length d binary vectors
    bv <- as.matrix(expand.grid(replicate(d, 0:1, simplify = FALSE)))
    # Alternating +/- 1
    coeff <- (-1)^(apply(bv, 1, sum))
    # Convert to logical
    bv <- (bv == 1)

    # Pre-evaluate margins at samples
    FX <- sapply(1:d, function(k) { margins[[k]](X[, k]) })
    # Pre-evaluate left limits of margins at samples
    FXm <- sapply(1:d, function(k) { margins[[k]](X[, k] - 1) })

    # Compute full density
    # P(X=X_i)=sum_{j1,..,jd=0,1}(-1)^(j1+..+jd) C(F1(X_i1-j1),..,Fd(X_id-jd)))
    P <- unlist(parallel::mclapply(1:n, function(i) {
        s <- 0
        # Evaluate density for each sample (finite difference with 2^d terms)
        for (j in 1:2^d) {
            FY <- FX[i, ]
            FY[bv[j, ]] <- FXm[i, bv[j, ]]
            s <- s + coeff[j] * copula::pCopula(FY, copula)
        }
        return(s)
    }))

    # Compute marginal densities
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

#' Sample cells
#'
#' @description Sample cells from the model specified by the given copula and
#' marginal quantile functions.
#'
#' @param N Number of cells to sample.
#' @param Qx List of marginal quantile functions.
#' @param copula A copula object. Currently supported are \code{vinecop},
#' \code{gaussiancop}, \code{tcop}, and \code{indepcop}.
#'
#' @returns A SingleCellExperiment.
#'
#' @export
sampleCells <- function(N,
                        Qx,
                        copula) {
    if (!is.list(Qx)) {
        stop("'Qx' must be a list.")
    }
    if (length(Qx) != length(copula$names)) {
        stop("Length of 'Qx' must equal the number of genes in 'copula'.")
    }
    if (N < 1) {
        stop("'N' must be at least 1.")
    }

    # Draw sample
    if (inherits(copula, "vinecop")) {
        U <- rvinecopulib::rvinecop(N, copula)
    } else if (inherits(copula, "gaussiancop") ||
               inherits(copula, "tcop") ||
               inherits(copula, "indepcop")) {
        U <- copula::rCopula(N, copula$copula)
    } else {
        stop("'", class(copula)[1], "' is not a supported copula class.")
    }

    # Convert margins
    genes <- copula$names
    X <- lapply(seq_along(genes), function(i) { Qx[[i]](U[, i]) })
    X <- do.call(cbind, X)

    # Construct new SingleCellExperiment
    rownames(X) <- paste0("sim_cell", seq(N))
    colnames(X) <- genes
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = t(X)))
    SingleCellExperiment::logcounts(sce) <-
        log1p(SingleCellExperiment::counts(sce))
    return(sce)
}

###############################################################################

#' Fit an independence copula to a SingleCellExperiment
#'
#' @param sce A SingleCellExperiment.
#'
#' @returns An \code{indepcop} object.
#'
#' @export
fitIndep <- function(sce) {
    if (!inherits(sce, "SingleCellExperiment")) {
        stop("'sce' must be a SingleCellExperiment.")
    }
    # Construct S3 object
    d <- dim(sce)[1]
    obj <- list(copula = copula::indepCopula(dim = d),
                nobs = dim(sce)[2],
                dim = d,
                names = rownames(sce))
    class(obj) <- "indepcop"
    return(obj)
}

###############################################################################

#' Print method for S3 class indepcop
#'
#' @param x Object of class \code{indepcop}.
#' @param ... Additional arguments.
#'
#' @export
#' @method print indepcop
print.indepcop <- function(x, ...) {
    l1 <- paste0(x$dim, "-dimensional independence copula fit")
    l2 <- paste("nobs =", x$nobs)
    cat(paste(l1, l2, sep = "\n"))
}

###############################################################################
