suppressMessages(library(dplyr))
suppressMessages(library(scCopula))
suppressMessages(library(parallel))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(fasano.franceschini.test))
set.seed(0, kind = "L'Ecuyer-CMRG")

setwd("~/copgen")

cores <- as.integer(commandArgs(trailingOnly = TRUE)[1])
ix <- as.integer(commandArgs(trailingOnly = TRUE)[2])

refs <- list.dirs("Results/02copulas", recursive = FALSE, full.names = FALSE)
ref <- refs[ix]

families <- c("norm", "vine", "nmle", "t")
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

    res_i <- NULL
    for (j in seq(nsample)) {
        fname <- paste0("Results/03samples/", ref, "/", family, "_", trial, "-", j, ".rds")
        sce_sim <- readRDS(fname)
        X1 <- as.matrix(t(counts(sce_test)))
        X2 <- as.matrix(t(counts(sce_sim)))
        pc <- proj(X = X1, Y = X2)
        ff2  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:2],  S2 = pc$Y_pc[, 1:2], nPermute = 0)$stat  / (2 * N)^(3/2)
        ff3  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:3],  S2 = pc$Y_pc[, 1:3], nPermute = 0)$stat  / (2 * N)^(3/2)
        ff4  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:4],  S2 = pc$Y_pc[, 1:4], nPermute = 0)$stat  / (2 * N)^(3/2)
        ff5  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:5],  S2 = pc$Y_pc[, 1:5], nPermute = 0)$stat  / (2 * N)^(3/2)
        ff6  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:6],  S2 = pc$Y_pc[, 1:6], nPermute = 0)$stat  / (2 * N)^(3/2)
        ff7  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:7],  S2 = pc$Y_pc[, 1:7], nPermute = 0)$stat  / (2 * N)^(3/2)
        ff8  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:8],  S2 = pc$Y_pc[, 1:8], nPermute = 0)$stat  / (2 * N)^(3/2)
        ff9  <- fasano.franceschini.test(S1 = pc$X_pc[, 1:9],  S2 = pc$Y_pc[, 1:9], nPermute = 0)$stat  / (2 * N)^(3/2)
        ff10 <- fasano.franceschini.test(S1 = pc$X_pc[, 1:10], S2 = pc$Y_pc[, 1:10], nPermute = 0)$stat / (2 * N)^(3/2)
        ffall <- fasano.franceschini.test(S1 = X1, S2 = X2, nPermute = 0)$stat / (2 * N)^(3/2)
        res_i <- rbind(res_i, data.frame(
            ref = ref,
            family = family,
            trial = trial,
            sample = j,
            ngene = dim(sce_sim)[1],
            ncell = dim(sce_sim)[2],
            ff2 = ff2,
            ff3 = ff3,
            ff4 = ff4,
            ff5 = ff5,
            ff6 = ff6,
            ff7 = ff7,
            ff8 = ff8,
            ff9 = ff9,
            ff10 = ff10,
            ffall = ffall
        ))
    }
    return(res_i)
}, mc.cores = cores)
res <- do.call(rbind, res)
rownames(res) <- NULL
saveRDS(res, file = paste0("Results/03samples/pcaall_", ix, ".rds"))
message("Done!")
