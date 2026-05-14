###############################################################################

#' Pairwise normalized mutual information matrix
#'
#' @description Compute upper triangular matrix of pairwise normalized mutual
#' information.
#'
#' @param X A matrix.
#'
#' @returns A matrix.
nmi_mat <- function(X) {
    n <- dim(X)[1]
    d <- dim(X)[2]

    entropy <- function(x) {
        p <- table(x) / length(x)
        -sum(p * log(p))
    }

    mi <- function(x, y) {
        table_xy <- table(x, y)
        pxy <- table_xy / length(x)
        px <- rowSums(pxy)
        py <- colSums(pxy)
        terms <- as.vector(pxy * log(pxy / (px %*% t(py))))
        sum(terms[is.finite(terms)])
    }

    nmi <- function(x, y) {
        hx <- entropy(x)
        hy <- entropy(y)
        if (hx == 0 || hy == 0) {
            return(0)
        }
        mi(x, y) / sqrt(hx * hy)
    }

    mat <- matrix(0, nrow = d, ncol = d)
    for (i in 1:(d - 1)) {
        for (j in (i + 1):d) {
            mat[[i, j]] <- mat[[j, i]] <- nmi(X[, i], X[, j])
        }
    }
    diag(mat) <- 1

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
            mat[[i, j]] <- mat[[j, i]] <- energy::dcor(X[, i], X[, j])
        }
    }
    diag(mat) <- 1
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
                                    "nmi",
                                    "bicor",
                                    "dcor")) {
    stats <- match.arg(stats, several.ok = TRUE)

    X1 <- Matrix::t(as.matrix(SingleCellExperiment::counts(sce1)))
    X2 <- Matrix::t(as.matrix(SingleCellExperiment::counts(sce2)))

    erf <- function(A, B, fun) {
        fA <- fun(A)
        fB <- fun(B)
        return(sqrt(sum((fA - fB)^2)) / sqrt(sum(fB^2)))
    }

    err <- c()
    for (i in seq_along(stats)) {
        cx <- stats[i]
        if (cx %in% c("pearson", "spearman", "kendall")) {
            fun <- function(Y) { stats::cor(Y, method = cx) }
        } else if (cx == "nmi") {
            fun <- nmi_mat
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
