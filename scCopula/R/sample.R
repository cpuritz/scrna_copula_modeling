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
    logcounts <- log1p(SingleCellExperiment::counts(sce))
    SingleCellExperiment::logcounts(sce) <- logcounts
    return(sce)
}

###############################################################################
