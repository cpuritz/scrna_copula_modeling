suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scCopula))
suppressMessages(library(parallel))
set.seed(0, kind = "L'Ecuyer-CMRG")

setwd("~/copgen")

cores <- as.integer(commandArgs(trailingOnly = TRUE)[1])

refs <- list.dirs("Results/02copulas", recursive = FALSE, full.names = FALSE)
families <- c("norm", "norm_jitter", "vine", "vine_jitter", "nmle", "t")
ntrial <- 50L
nsample <- 20L

for (ref_ix in seq_along(refs)) {
    ref <- refs[ref_ix]
    message("Building ", ref, " ", ref_ix, "/", length(refs))
    sce <- readRDS(paste0("Data/References/Subsets/", ref, ".rds"))
    shuffles <- readRDS(paste0("Results/01shuffles/", ref, ".rds"))
    N <- floor(dim(sce)[2] / 2)
    
    sims <- expand.grid(family = families, trial = seq(ntrial))
    sink <- mclapply(seq(dim(sims)[1]), function(i) {
        family <- sims[i, "family"]
        trial <- sims[i, "trial"]

        ifname <- paste0("Results/02copulas/", ref, "/", family, "_", trial, ".rds") 
        if (file.exists(ifname)) {
            # Read in copula
            cop <- readRDS(ifname)
        
            # Train and test sets
            sce_train <- sce[, shuffles[trial, 1:N]]
            sce_test <- sce[, shuffles[trial, (N + 1):(2 * N)]]
        
            # Construct empirical quantile functions
            Qx <- apply(as.matrix(counts(sce_train)), 1, function(x) {
                function(p) { quantile(x, probs = p, type = 1) }
            })
        
            # Draw samples
            for (j in seq(nsample)) {
                ofname <- paste0("Results/03samples/", ref, "/", family, "_", trial, "-", j, ".rds")
                sce_sim <- sampleCells(N, Qx, cop)
                saveRDS(sce_sim, file = ofname)
            }
        }
        return(NULL)
    }, mc.cores = cores)
}
