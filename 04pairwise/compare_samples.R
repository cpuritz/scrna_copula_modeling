library(dplyr)
library(scCopula)
library(parallel)
library(SingleCellExperiment)
set.seed(0, kind = "L'Ecuyer-CMRG")

args <- commandArgs(trailingOnly = TRUE)
cores <- as.integer(args[1])
ix <- as.integer(args[2])
stat <- args[3]
family <- args[4]

refs <- list.dirs("Results/02copulas", recursive = FALSE, full.names = FALSE)
ref <- refs[ix]

message("Computing stats for ", stat, ", ", ref, ", ", family)
ntrial <- 50L
nsample <- 20L

sce <- readRDS(paste0("Data/References/Subsets/", ref, ".rds"))
shuffles <- readRDS(paste0("Results/01shuffles/", ref, ".rds"))
N <- floor(dim(sce)[2] / 2)
            
sims <- expand.grid(family = family, trial = seq(ntrial))
res <- mclapply(seq(dim(sims)[1]), function(i) {
    family <- sims[i, "family"]
    trial <- sims[i, "trial"]
    sce_test <- sce[, shuffles[trial, (N + 1):(2 * N)]]

    res_i <- NULL
    for (j in seq(nsample)) {
        sce_sim <- readRDS(paste0("Results/03samples/", ref, "/", family, "_", trial, "-", j, "_par.rds"))
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
saveRDS(res, file = paste0("Results/04pairwise/par_comparisons_", stat, "_", ix, "_", family, ".rds"))
