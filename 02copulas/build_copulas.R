suppressMessages(library(SingleCellExperiment))
suppressMessages(library(parallel))
suppressMessages(library(scCopula))
suppressMessages(library(copula))
set.seed(1, kind = "L'Ecuyer-CMRG")

setwd("~/copgen")

# Read inputs
ref <- commandArgs(trailingOnly = TRUE)[1]
family <- commandArgs(trailingOnly = TRUE)[2]
cores <- as.integer(commandArgs(trailingOnly = TRUE)[3])

ref <- unlist(strsplit(ref, ".rds"))

# Number of trials to run
ntrial <- 50L

# Read in dataset
sce <- readRDS(paste0("Data/References/Subsets/", ref, ".rds"))

# Read in train-test shuffles
shuffles <- readRDS(paste0("Results/01shuffles/", ref, ".rds"))

# Use half of cells for training and half for testing
N <- floor(dim(sce)[2] / 2)

sink <- mclapply(1:ntrial, function(i) {
    sce_train <- sce[, shuffles[i, 1:N]]
    if (family == "vine") {
        cop <- fitVine(sce = sce_train, margins = "empirical")
    } else if (family == "vine_jitter") {
        cop <- fitVine(sce = sce_train, margins = "empirical", jitter = TRUE)
    } else if (family == "norm") {
        cop <- fitGaussian(sce = sce_train, margins = "empirical", mle = FALSE, likelihood = FALSE)
    } else if (family == "norm_jitter") {
        cop <- fitGaussian(sce = sce_train, margins = "empirical", mle = FALSE, jitter = TRUE, likelihood = FALSE)
    } else if (family == "nmle") {
        cop <- fitGaussian(sce = sce_train, margins = "empirical", mle = TRUE, likelihood = FALSE)
    } else if (family == "t") {
        fname <- paste0("Results/02copulas/", ref, "/nmle_", i, ".rds")
        if (!file.exists(fname)) {
            message("Gaussian copula not fit")
            return(NULL)
        }
        cop_norm <- readRDS(fname)
        Sigma <- getSigma(cop_norm$copula)
        cop <- fitT(sce = sce_train, margins = "empirical", Sigma = Sigma, likelihood = FALSE)
    }
    saveRDS(cop, paste0("Results/02copulas/", ref, "/", family, "_", i, ".rds"))
    
    return(NULL)
}, mc.cores = cores)
message("Done!")
