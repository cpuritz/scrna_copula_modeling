###############################################################################

#' Pairwise mutual information matrix
#'
#' @description Compute upper triangular matrix of pairwise mutual information.
#'
#' @param X A matrix.
#'
#' @returns A matrix.
mi_mat <- function(X) {
    n <- dim(X)[1]
    d <- dim(X)[2]
    mat <- matrix(nrow = d, ncol = d)
    mi <- function(a, b) {
        table_xy <- table(a, b)
        pxy <- table_xy / n
        px <- rowSums(table_xy) / n
        py <- colSums(table_xy) / n
        sum(stats::na.omit(as.vector(pxy * log(pxy / px %*% t(py)))))
    }
    for (i in 1:(d - 1)) {
        for (j in (i + 1):d) {
            mat[[i, j]] <- mi(X[, i], X[, j])
        }
    }
    return(mat)
}

###############################################################################

#' Pairwise distance correlation matrix
#'
#' @description Compute upper triangular matrix of pairwise distance
#' correlations.
#'
#' @param X A matrix.
#'
#' @returns A matrix.
dcor_mat <- function(X) {
    d <- dim(X)[2]
    mat <- matrix(nrow = d, ncol = d)
    for (i in 1:(d - 1)) {
        for (j in (i + 1):d) {
            mat[[i, j]] <- energy::dcor(X[, i], X[, j])
        }
    }
    return(mat)
}

###############################################################################

#' Count Matrix Differences
#'
#' @description Compute various statistics comparing two count matrices.
#'
#' @param sce1 A SingleCellExperiment.
#' @param sce2 A SingleCellExperiment.
#' @param stats Which statistics to compute.
#'
#' @returns A vector of statistics.
#'
#' @export
compareCounts <- function(sce1,
                          sce2,
                          stats = c("pearson",
                                    "spearman",
                                    "kendall",
                                    "mi",
                                    "bicor",
                                    "dcor")) {
    stats <- match.arg(stats, several.ok = TRUE)

    X1 <- Matrix::t(as.matrix(SingleCellExperiment::counts(sce1)))
    X2 <- Matrix::t(as.matrix(SingleCellExperiment::counts(sce2)))

    erf <- function(A, B, fun) {
        fA <- fun(A)
        fB <- fun(B)
        vfA <- fA[upper.tri(fA)]
        vfB <- fB[upper.tri(fB)]
        return(sqrt(mean((vfA - vfB)^2)))
    }

    err <- c()
    for (i in seq_along(stats)) {
        cx <- stats[i]
        if (cx %in% c("pearson", "spearman", "kendall")) {
            fun <- function(Y) { stats::cor(Y, method = cx) }
        } else if (cx == "mi") {
            fun <- mi_mat
        } else if (cx == "bicor") {
            fun <- WGCNA::bicor
        } else if (cx == "dcor") {
            fun <- dcor_mat
        }
        err[i] <- erf(X1, X2, fun)
    }
    res <- data.frame(stats = stats, err = err)

    return(res)
}

###############################################################################
