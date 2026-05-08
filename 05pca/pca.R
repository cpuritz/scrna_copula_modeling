library(dplyr)
library(scCopula)
library(parallel)
library(SingleCellExperiment)
library(fasano.franceschini.test)
set.seed(0, kind = "L'Ecuyer-CMRG")

args <- commandArgs(trailingOnly = TRUE)
cores <- as.integer(args[1])
ix <- as.integer(args[2])
npc <- as.integer(args[3])

refs <- list.dirs("Results/02copulas", recursive = FALSE, full.names = FALSE)
ref <- refs[ix]

message("Running for reference ", ix, ", npc = ", npc)

families <- c("ind", "norm", "norm_jitter", "vine", "vine_jitter", "nmle", "t", "zinbwave", "sparsim")
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
        sce_sim <- readRDS(paste0("Results/03samples/", ref, "/", family, "_", trial, "-", j, "_par.rds"))
        X2 <- as.matrix(t(counts(sce_sim)))
        pc <- proj(X = X1, Y = X2)

        p <- fasano.franceschini.test(
            S1 = pc$X_pc[, 1:npc], 
            S2 = pc$Y_pc[, 1:npc],
            verbose = FALSE,
            seed = 10 * j + npc
        )$p.value

        res_i <- rbind(res_i, data.frame(
            ref = ref,
            family = family,
            trial = trial,
            sample = j,
            npc = npc,
            p = as.numeric(p)    
        ))
    }
    return(res_i)
}, mc.cores = cores)
res <- do.call(rbind, res)
rownames(res) <- NULL
saveRDS(res, file = paste0("Results/05pca/pca_", ix, "_", npc, ".rds"))
