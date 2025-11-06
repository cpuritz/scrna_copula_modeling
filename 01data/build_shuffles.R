library(SingleCellExperiment)
set.seed(1)

refs <- list.files("Data/References/Subsets")
ntrial <- 50L

for (ref in refs) {
    sce <- readRDS(paste("Data/References/Subsets", ref, sep = '/'))
    shuffles <- t(replicate(ntrial, sample(dim(sce)[2])))
    saveRDS(shuffles, file = paste("Results/01shuffles", ref, sep = '/'))
}
