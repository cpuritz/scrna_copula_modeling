###############################################################################

#' Pairwise mutual information matrix
#'
#' @description Compute matrix of pairwise mutual information.
#'
#' @param X A matrix.
#'
#' @returns A matrix.
mut_inf <- function(X) {
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
    for (i in seq(d)) {
        mat[[i, i]] <- mi(X[, i], X[, i])
    }
    for (i in seq(1, d - 1)) {
        for (j in seq(i + 1, d)) {
            mat[[i, j]] <- mat[[j, i]] <- mi(X[, i], X[, j])
        }
    }
    return(mat)
}

###############################################################################

#' Pairwise distance correlation matrix
#'
#' @description Compute matrix of pairwise distance correlations.
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
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
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
#' @param assay The assay to use.
#' @param stats Which statistics to compute.
#'
#' @returns A vector of statistics.
#'
#' @export
compareCounts <- function(sce1,
                          sce2,
                          assay = "counts",
                          stats = c("pearson",
                                    "spearman",
                                    "kendall",
                                    "mi",
                                    "bicor",
                                    "dcor",
                                    "rhoA")) {
    if (!assay %in% names(SummarizedExperiment::assays(sce1))) {
        stop("Assay ", assay, " not found in sce1.")
    }
    if (!assay %in% names(SummarizedExperiment::assays(sce2))) {
        stop("Assay ", assay, " not found in sce2.")
    }
    stats <- match.arg(stats, several.ok = TRUE)

    X1 <- Matrix::t(as.matrix(SummarizedExperiment::assay(sce1, assay)))
    X2 <- Matrix::t(as.matrix(SummarizedExperiment::assay(sce2, assay)))

    # Compute L1, L2, and Linf matrix norms
    errs <- function(A, B, name) {
        d <- abs(A - B)
        data.frame(
            stat = rep(name, 3),
            norm = c("one", "two", "inf"),
            err = c(sum(d), sqrt(sum(d^2)), max(d))
        )
    }

    res <- NULL
    for (cx in stats) {
        if (cx %in% c("pearson", "spearman", "kendall")) {
            fun <- function(Y) { stats::cor(Y, method = cx) }
        } else if (cx == "mi") {
            fun <- mut_inf
        } else if (cx == "bicor") {
            fun <- WGCNA::bicor
        } else if (cx == "dcor") {
            fun <- dcor_mat
        } else if (cx == "rhoA") {
            fun <- rhoA_mat
        }
        res <- rbind(res, errs(fun(X1), fun(X2), cx))
    }

    return(res)
}

###############################################################################

#' Pairwise zero-inflated Spearman's rho
#'
#' @description Compute matrix of pairwise zero-inflated Spearman's rho.
#'
#' @param X A matrix.
#'
#' @returns A matrix.
rhoA_mat <- function(X) {
    d <- dim(X)[2]
    mat <- matrix(nrow = d, ncol = d)
    for (i in 1:d) {
        for (j in i:d) {
            mat[[i, j]] <- mat[[j, i]] <- spm_Arends(X[, i], X[, j])
        }
    }
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
    diag(mat) <- 1
    return(mat)
}

###############################################################################

spm_Arends <- function(A, B) {
	# Modified from https://github.com/JasperArends/SpmZID_Rcode

    N <- length(A)

    p00 <- sum(A == 0 & B == 0) / N  # P(X = 0, Y = 0)
    p01 <- sum(A == 0 & B > 0) / N   # P(X = 0, Y > 0)
    p10 <- sum(A > 0 & B == 0) / N   # P(X > 0, Y = 0)
    p11 <- sum(A > 0 & B > 0) / N    # P(x > 0, Y > 0)

    ######################################
    # Margins
    p0m <- p00 + p01  # P(X = 0)
    p1m <- p10 + p11  # P(X > 0)
    pm0 <- p00 + p10  # P(Y = 0)
    pm1 <- p01 + p11  # P(Y > 0)


    ######################################
    # Split data
    X10 <- which(A > 0 & B == 0)
    X01 <- which(A == 0 & B > 0)
    X11 <- which(A > 0 & B > 0)

    n11 <- length(X11)
    n10 <- length(X10)
    n01 <- length(X01)

    # Initialize probabilities
    p1s <- 0  # P(X10 > X11)
    p1c <- 0  # P(X10 = X11)
    p2s <- 0  # P(Y01 > Y11)
    p2c <- 0  # P(Y01 = Y11)

    for (ix in X10) {
        p1s <- p1s + sum(A[ix] > A[X11])
        p1c <- p1c + sum(A[ix] == A[X11])
    }
    p1s <- p1s / max(n10 * n11, 1)
    p1c <- p1c / max(n10 * n11, 1)

    for (ix in X01) {
        p2s <- p2s + sum(B[ix] > B[X11])
        p2c <- p2c + sum(B[ix] == B[X11])
    }
    p2s <- p2s / max(n01 * n11, 1)
    p2c <- p2c / max(n01 * n11, 1)

    D <- p11 * (p10 * (1 - 2 * p1s - p1c) + p01 * (1 - 2 * p2s - p2c))

    ######################################
    # Estimate P(C | C1111) - P(D | C1111)
    # Probabilities are stored in pC1, pC2, pC3 and pC4 such that pC1[1]
    # and pC1[2] correspond to concordance, discordance and ties for the
    # first case, where X3 > 0 and Y2 > 0 (similarly for pC2, pC3 and pC4).
    pC2 <- c(0, 0) # X3 = 0, Y2 > 0
    pC3 <- c(0, 0) # X3 > 0, Y2 = 0
    pC4 <- c(0, 0) # X3 = 0, Y2 = 0

    nX <- c(0, 0, 0, 0)
    nY <- c(0, 0, 0, 0)
    nXdb <- c(0, 0)
    for (ix in X11) {
        X1 <- A[ix]
        Y1 <- B[ix]

        nX[1] <- sum(X1 > A[X11])
        nX[2] <- sum(X1 > A[X10])
        nX[3] <- sum(X1 < A[X11])
        nX[4] <- sum(X1 < A[X10])
        nY[1] <- sum(Y1 > B[X11])
        nY[2] <- sum(Y1 > B[X01])
        nY[3] <- sum(Y1 < B[X11])
        nY[4] <- sum(Y1 < B[X01])

        pC2[1] <- pC2[1] + nX[1] * nY[2] + nX[3] * nY[4]
        pC2[2] <- pC2[2] + nX[1] * nY[4] + nX[3] * nY[2]

        pC3[1] <- pC3[1] + nX[2] * nY[1] + nX[4] * nY[3]
        pC3[2] <- pC3[2] + nX[2] * nY[3] + nX[4] * nY[1]

        pC4[1] <- pC4[1] + nX[2] * nY[2] + nX[4] * nY[4]
        pC4[2] <- pC4[2] + nX[2] * nY[4] + nX[4] * nY[2]
    }

    # Concordance - discordance
    pC <- c((pC2[1] - pC2[2]),
            (pC3[1] - pC3[2]),
            (pC4[1] - pC4[2]))
    # Convert to probabilities
    cnt <- c(as.double(n11) * (n11 - 1) * n01,
             as.double(n11) * n10 * (n11 - 1),
             as.double(n11) * n10 * n01)
    pC <- pC / pmax(cnt, 1) # Avoid division by 0

    rho11 <- stats::cor(A[X11], B[X11], method = "spearman") / 3
    if (is.na(rho11)) {
        rho11 <- 0
    }

    pPos <- p11 * (p11^2 * rho11 + p01 * p11 * pC[1] +
                       p11 * p10 * pC[2] + p01 * p10 * pC[3])

    # Spearman's rho estimate
    rho_est <- 3 * pPos + 3 * (p00 * p11 - p10 * p01) + 3 * D
    return(rho_est)
}

###############################################################################
