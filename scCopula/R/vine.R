###############################################################################

#' Fit a vine copula to a SingleCellExperiment
#'
#' @param sce A SingleCellExperiment.
#' @param margins Type of marginal distribution. Either "empirical" or "par". If
#' the former, margins are modeled empirically. If the latter, margins are
#' modeled as negative binomial/zero-inflated negative binomial (selected by
#' AIC).
#' @param jitter Logical indicating whether data should be jittered.
#' @param family_set Family set to use for vine copula construction.
#' @param cores Number of cores for parallel computation.
#'
#' @returns A vine copula object.
#'
#' @export
fitVine <- function(sce,
                    margins = c("empirical", "par"),
                    jitter = FALSE,
                    family_set = c("indep", "gaussian", "clayton",
                                   "gumbel", "frank", "joe"),
                    cores = 1L) {
    if (!inherits(sce, "SingleCellExperiment")) {
        stop("'sce' must be a SingleCellExperiment.")
    }
    ngene <- dim(sce)[1]
    margins <- match.arg(margins)
    family_set <- match.arg(family_set, several.ok = TRUE)
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
