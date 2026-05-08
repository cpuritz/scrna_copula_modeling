library(SingleCellExperiment)
set.seed(1)

refs_all <- list.files("Data/References/Subsets")

refs <- refs_all[grepl("hsa|HVG", refs_all)]
for (ref in refs) {
    sce <- readRDS(paste("Data/References/Subsets", ref, sep = '/'))
    shuffles <- t(replicate(50L, sample(dim(sce)[2])))
    saveRDS(shuffles, file = paste("Results/01shuffles", ref, sep = '/'))
}

refs <- refs_all[grepl("hdwgcna", refs_all)]
for (ref in refs) {
    sce <- readRDS(paste("Data/References/Subsets", ref, sep = '/'))
    shuffles <- t(replicate(10L, sample(dim(sce)[2])))
    saveRDS(shuffles, file = paste("Results/01shuffles", ref, sep = '/'))
}
