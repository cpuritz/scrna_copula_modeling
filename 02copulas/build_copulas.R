library(SingleCellExperiment)
library(parallel)
library(scCopula)
library(copula)
set.seed(1, kind = "L'Ecuyer-CMRG")

# Read inputs
ref_ix <- as.integer(commandArgs(trailingOnly = TRUE)[1])
family <- commandArgs(trailingOnly = TRUE)[2]
cores <- as.integer(commandArgs(trailingOnly = TRUE)[3])

refs_all <- list.files("Data/References/Subsets")
refs <- refs_all[grepl("hsa|HVG", refs_all)]
ref <- unlist(strsplit(refs[ref_ix], ".rds"))

# Number of trials to run
ntrial <- 50L

# Read in dataset
sce <- readRDS(paste0("Data/References/Subsets/", ref, ".rds"))

# Read in train-test shuffles
shuffles <- readRDS(paste0("Results/01shuffles/", ref, ".rds"))

# Use half of cells for training and half for testing
N <- floor(dim(sce)[2] / 2)

margins <- "par"
sink <- mclapply(seq(ntrial), function(i) {
    sce_train <- sce[, shuffles[i, 1:N]]
    if (family == "vine") {
        cop <- fitVine(sce = sce_train, margins = margins, jitter = FALSE)
    } else if (family == "vine_jitter") {
        cop <- fitVine(sce = sce_train, margins = margins, jitter = TRUE)
    } else if (family == "norm") {
        cop <- fitGaussian(sce = sce_train, margins = margins, estimate = "sample")
    } else if (family == "norm_jitter") {
        cop <- fitGaussian(sce = sce_train, margins = margins, estimate = "jitter")
    } else if (family == "nmle") {
        cop <- fitGaussian(sce = sce_train, margins = margins, estimate = "mle")
    } else if (family == "t") {
        fname <- paste0("Results/02copulas/", ref, "/nmle_", i, "_par.rds")
        if (!file.exists(fname)) {
            message("Gaussian copula not fit")
            return(NULL)
        }
        cop_norm <- readRDS(fname)
        Sigma <- getSigma(cop_norm$copula)
        cop <- fitT(sce = sce_train, margins = margins, Sigma = Sigma)
    } else if (family == "ind") {
        cop <- fitIndep(sce_train)
    }
    saveRDS(cop, paste0("Results/02copulas/", ref, "/", family, "_", i, "_par.rds"))
    return(NULL)
}, mc.cores = cores)
