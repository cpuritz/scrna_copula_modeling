library(SingleCellExperiment)
set.seed(1)

refs_all <- list.files("Data/References/Subsets")

ntrial <- 50L
refs <- refs_all[grepl("hsa|HVG", refs_all)]
for (ref in refs) {
    sce <- readRDS(paste("Data/References/Subsets", ref, sep = '/'))
    shuffles <- t(replicate(ntrial, sample(dim(sce)[2])))
    saveRDS(shuffles, file = paste("Results/01shuffles", ref, sep = '/'))
}

ntrial <- 10L
refs <- refs_all[grepl("hdwgcna", refs_all)]
for (ref in refs) {
    sce <- readRDS(paste("Data/References/Subsets", ref, sep = '/'))
    shuffles <- t(replicate(ntrial, sample(dim(sce)[2])))
    saveRDS(shuffles, file = paste("Results/01shuffles", ref, sep = '/'))
}
