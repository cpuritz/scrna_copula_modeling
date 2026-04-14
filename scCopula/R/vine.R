###############################################################################

#' Fit a vine copula to a SingleCellExperiment
#'
#' @param sce A SingleCellExperiment.
#' @param margins One, or a list of, marginal distribution functions. If only
#' one function is passed, all margins are modeled identically. Entries can also
#' equal "empirical" or "nb", in which case the corresponding margin is modeled
#' empirically or as a negative binomial, respectively.
#' @param jitter Logical indicating whether data should be jittered.
#' @param family_set Family set to use for vine copula construction.
#' @param cores Number of cores for parallel computation.
#'
#' @returns A vine copula object.
#'
#' @export
fitVine <- function(sce,
                    margins,
                    jitter = FALSE,
                    family_set = c("indep", "gaussian", "clayton",
                                   "gumbel", "frank", "joe"),
                    cores = 1L) {
    if (!inherits(sce, "SingleCellExperiment")) {
        stop("'sce' must be a SingleCellExperiment.")
    }
    ngene <- dim(sce)[1]
    if (length(margins) == 1) {
        margins <- rep(list(margins), ngene)
    }
    if (length(margins) != ngene) {
        stop("Length of 'margins' must equal 1 or the number of genes in",
             " 'sce'.")
    }
    family_set <- match.arg(family_set, several.ok = TRUE)
    cores <- as.integer(cores)

    # Count matrix
    X <- as.matrix(Matrix::t(SingleCellExperiment::counts(sce)))

    # Construct marginal distribution functions if necessary
    for (i in 1:ngene) {
        if (is.character(margins[[i]])) {
            if (margins[[i]] == "empirical") {
                margins[[i]] <- empcdf(X[, i])
            } else if (margins[[i]] == "nb") {
                par <- fitdistrplus::fitdist(
                    data = X[, i],
                    distr = "nbinom",
                    method = "mle"
                )$estimate
                margins[[i]] <- function(q) {
                    do.call(stats::pnbinom, c(as.list(par), list(q = q)))
                }
            }
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
        V <- matrix(stats::runif(prod(dim(FX))), nrow = nrow(FX))
        U <- FXm + (FX - FXm) * V
        var_types <- rep("c", ngene)
    } else {
        U <- cbind(FX, FXm)
        var_types <- rep("d", ngene)
    }
    colnames(U) <- NULL

    # Fit vine copula
    vine <- rvinecopulib::vinecop(
        data = U,
        var_types = var_types,
        family_set = family_set,
        par_method = "mle",
        selcrit = "aic",
        cores = cores
    )

    # Store gene names
    vine$names <- colnames(X)
    return(vine)
}

###############################################################################
