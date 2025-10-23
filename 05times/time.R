suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scCopula))
suppressMessages(library(parallel))
set.seed(0)
setwd("~/copgen")

family <- commandArgs(trailingOnly = TRUE)[1]
task <- as.integer(commandArgs(trailingOnly = TRUE)[2])
cores <- as.integer(commandArgs(trailingOnly = TRUE)[3])

sce_all <- readRDS("Data/References/Time/bailey24-tram1.rds")

genes <- c(10, 20, 30, 40, 50)
cells <- c(100, 1000, 2000, 3000, 4000, 5000)
params <- expand.grid(genes, cells)

ng <- params[task, 1]
nc <- params[task, 2]

message("Timing for family = ", family, ", genes = ", ng, ", cells = ", nc)
nrun <- 100L
times <- parallel::mclapply(seq(nrun), function(i) {
    # Randomly subsample
    cells <- sample(dim(sce_all)[2], nc)
    sce <- sce_all[, cells]
    sce <- sce[rowSums(counts(sce)) > 10, ]
    genes <- sample(dim(sce)[1], ng)
    sce <- sce[genes, ]
    
    X <- as.matrix(t(counts(sce)))
    margins <- lapply(seq(ng), function(i) { empcdf(X[, i]) })
    
    if (family == "vine") {
        t0 <- Sys.time()
        res <- fitVine(sce, margins, jitter = FALSE)
        t1 <- Sys.time()
    } else if (family == "vine_jitter") {
        t0 <- Sys.time()
        res <- fitVine(sce, margins, jitter = TRUE)
        t1 <- Sys.time()
    } else if (family == "norm") {
        t0 <- Sys.time()
        res <- fitGaussian(sce, margins, jitter = FALSE, mle = FALSE, likelihood = FALSE)
        t1 <- Sys.time()
    } else if (family == "norm_jitter") {
        t0 <- Sys.time()
        res <- fitGaussian(sce, margins, jitter = TRUE, mle = FALSE, likelihood = FALSE)
        t1 <- Sys.time()
    } else if (family == "nmle") {
        t0 <- Sys.time()
        res <- fitGaussian(sce, margins, mle = TRUE, likelihood = FALSE)
        t1 <- Sys.time()
    } else if (family == "t") {
        t0 <- Sys.time()
        res <- fitT(sce, margins, Sigma = NULL, likelihood = FALSE)
        t1 <- Sys.time()
    }
    return(as.numeric(difftime(t1, t0, units = "secs")))
}, mc.cores = cores)

df <- data.frame(
    family = rep(family, nrun),
    genes = rep(ng, nrun),
    cells = rep(nc, nrun),
    time = unlist(times)
)
saveRDS(df, file = paste0("Results/04times/", family, "_", ng, "_", nc, ".rds"))
message("Done!")
