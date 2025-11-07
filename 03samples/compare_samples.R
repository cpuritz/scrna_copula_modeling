library(dplyr)
library(scCopula)
library(parallel)
library(SingleCellExperiment)
set.seed(0, kind = "L'Ecuyer-CMRG")

cores <- as.integer(commandArgs(trailingOnly = TRUE)[1])
ix <- as.integer(commandArgs(trailingOnly = TRUE)[2])
stat <- commandArgs(trailingOnly = TRUE)[3]

refs <- list.dirs("Results/02copulas", recursive = FALSE, full.names = FALSE)
ref <- refs[ix]
message("Computing stats for ", stat, ", ", ref)

families <- c("ind", "norm", "norm_jitter", "vine", "vine_jitter", "nmle", "t")
ntrial <- 50L
nsample <- 20L

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
        comp <- compareCounts(sce_sim, sce_test, stats = stat)
        comp$ref <- ref
        comp$family <- family
        comp$trial <- trial
        comp$sample <- j
        comp$ngene <- dim(sce_sim)[1]
        comp$ncell <- dim(sce_sim)[2]
        res_i <- rbind(res_i, comp)
    }
    return(res_i)
}, mc.cores = cores)
res <- do.call(rbind, res)
rownames(res) <- NULL
saveRDS(res, file = paste0("Results/03samples/comparisons_", stat, "_", ix, ".rds"))
