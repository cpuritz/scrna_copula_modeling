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
