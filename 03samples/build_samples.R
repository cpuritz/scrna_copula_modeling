library(SingleCellExperiment)
library(scCopula)
library(parallel)
library(zinbwave)
library(BiocParallel)
set.seed(0, kind = "L'Ecuyer-CMRG")

cores <- as.integer(commandArgs(trailingOnly = TRUE)[1])

refs <- list.dirs("Results/02copulas", recursive = FALSE, full.names = FALSE)
families <- c("ind", "norm", "norm_jitter", "vine", "vine_jitter", "nmle", "t")
simulators <- c("zinbwave")
methods <- c(families, simulators)
ntrial <- 50L
nsample <- 20L

ofname <- function(w, x, y, z) {
    paste0("Results/03samples/", w, "/", x, "_", y, "-", z, ".rds")
}

for (ref_ix in seq_along(refs)) {
    ref <- refs[ref_ix]
    sce <- readRDS(paste0("Data/References/Subsets/", ref, ".rds"))
    shuffles <- readRDS(paste0("Results/01shuffles/", ref, ".rds"))
    N <- floor(dim(sce)[2] / 2)
    
    sims <- expand.grid(method = methods, trial = seq(ntrial))
    sink <- mclapply(seq(dim(sims)[1]), function(i) {
        method <- sims[i, "method"]
        trial <- sims[i, "trial"]
        
        if (method %in% families) {
            # Read in copula
            cop <- readRDS(paste0("Results/02copulas/", ref, "/", method, "_",
                                  trial, ".rds") )
            
            # Train and test sets
            sce_train <- sce[, shuffles[trial, 1:N]]
            sce_test <- sce[, shuffles[trial, (N + 1):(2 * N)]]
            
            # Construct empirical quantile functions
            Qx <- apply(as.matrix(counts(sce_train)), 1, function(x) {
                function(p) { quantile(x, probs = p, type = 1) }
            })
            
            get_sample <- function() {
                sampleCells(N, Qx, cop)
            }
        } else if (method == "zinbwave") {
            # Train and test sets
            sce_train <- sce[, shuffles[trial, 1:N]]
            sce_test <- sce[, shuffles[trial, (N + 1):(2 * N)]]
            
            # Fit ZINB-WAVE model
            zinb_model <- zinbFit(
                sce_train, 
                verbose = FALSE,
                BPPARAM = SerialParam()
            )
            
            get_sample <- function() {
                SingleCellExperiment(
                    list(counts = zinbSim(zinb_model)$counts), 
                    colData = colData(sce_train)
                )
            }
        }
        
        # Draw samples
        for (j in seq(nsample)) {
            sce_sim <- get_sample()
            saveRDS(sce_sim, file = ofname(ref, method, trial, j))
        }
        
        return(NULL)
    }, mc.cores = cores)
}
