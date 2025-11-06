library(dplyr)
library(scCopula)
library(parallel)
library(SingleCellExperiment)
library(fasano.franceschini.test)
set.seed(0, kind = "L'Ecuyer-CMRG")

cores <- as.integer(commandArgs(trailingOnly = TRUE)[1])
ix <- as.integer(commandArgs(trailingOnly = TRUE)[2])

refs <- list.dirs("Results/02copulas", recursive = FALSE, full.names = FALSE)
ref <- refs[ix]

families <- c("ind", "norm", "norm_jitter", "vine", "vine_jitter", "nmle", "t")
ntrial <- 50L
nsample <- 20L

# Project X and Y onto the PCA space of X
proj <- function(X, Y) {
    pr <- prcomp(X)
    rot <- pr$rotation
    Y <- Y[, colnames(X)]
    Yc <- scale(Y, center = colMeans(X), scale = FALSE)
    return(list(X_pc = pr$x, Y_pc = Yc %*% rot))
}

sce <- readRDS(paste0("Data/References/Subsets/", ref, ".rds"))
shuffles <- readRDS(paste0("Results/01shuffles/", ref, ".rds"))
N <- floor(dim(sce)[2] / 2)
    
sims <- expand.grid(family = families, trial = seq(ntrial))
res <- mclapply(seq(dim(sims)[1]), function(i) {
    family <- sims[i, "family"]
    trial <- sims[i, "trial"]
    sce_test <- sce[, shuffles[trial, (N + 1):(2 * N)]]
    X1 <- as.matrix(t(counts(sce_test)))

    res_i <- NULL
    for (j in seq(nsample)) {
        fname <- paste0("Results/03samples/", ref, "/", family, "_", trial, "-", j, ".rds")
        sce_sim <- readRDS(fname)
        X2 <- as.matrix(t(counts(sce_sim)))
        pc <- proj(X = X1, Y = X2)

        seeds <- 10 * j + seq(0, 9)
        test2  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:2],  S2 = pc$Y_pc[, 1:2],  verbose = FALSE, seed = seeds[1])
        test3  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:3],  S2 = pc$Y_pc[, 1:3],  verbose = FALSE, seed = seeds[2])
        test4  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:4],  S2 = pc$Y_pc[, 1:4],  verbose = FALSE, seed = seeds[3])
        test5  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:5],  S2 = pc$Y_pc[, 1:5],  verbose = FALSE, seed = seeds[4])
        test6  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:6],  S2 = pc$Y_pc[, 1:6],  verbose = FALSE, seed = seeds[5])
        test7  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:7],  S2 = pc$Y_pc[, 1:7],  verbose = FALSE, seed = seeds[6])
        test8  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:8],  S2 = pc$Y_pc[, 1:8],  verbose = FALSE, seed = seeds[7])
        test9  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:9],  S2 = pc$Y_pc[, 1:9],  verbose = FALSE, seed = seeds[8])
        test10 <- fasano.franceschini.test(S1 = pc$X_pc[, 1:10], S2 = pc$Y_pc[, 1:10], verbose = FALSE, seed = seeds[9])

        res_i <- rbind(res_i, data.frame(
            ref = ref,
            family = family,
            trial = trial,
            sample = j,
            ngene = dim(sce_sim)[1],
            ncell = dim(sce_sim)[2],
            pval2 = as.numeric(test2$p.value),
            pval3 = as.numeric(test3$p.value),
            pval4 = as.numeric(test4$p.value),
            pval5 = as.numeric(test5$p.value)
            pval6 = as.numeric(test6$p.value),
            pval7 = as.numeric(test7$p.value),
            pval8 = as.numeric(test8$p.value)
            pval9 = as.numeric(test9$p.value),
            pval10 = as.numeric(test10$p.value)
        ))
    }
    return(res_i)
}, mc.cores = cores)
res <- do.call(rbind, res)
rownames(res) <- NULL
saveRDS(res, file = paste0("Results/04pca/pcap_", ix, ".rds"))
