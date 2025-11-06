library(Seurat)
library(SingleCellExperiment)
library(scran)
library(dplyr)
set.seed(0)
setwd("~/Documents/Graduate School/Copula Paper")

subs <- list(
    c("king21", "tonsil germinal center B cell", "BCP5"),
    c("king21", "naive B cell", "BCP6")
)
subs <- data.frame(do.call(rbind, subs))
colnames(subs) <- c("ref", "cell_type", "donor")
subs$nhvg <- floor(exp(runif(1000, min = log(500), max = log(1000))))

for (i in seq_len(dim(subs)[1])) {
    message(i, "/", dim(subs)[1])
    ref <- subs[i, "ref"]
    ct <- subs[i, "cell_type"]
    id <- subs[i, "donor"]

    obj <- readRDS(paste0("Data/Raw/", ref, ".rds"))
    obj <- obj[, obj@meta.data$cell_type == ct & obj@meta.data$donor_id == id]
    sce <- as.SingleCellExperiment(obj)
    sce <- sce[colSums(counts(sce) > 0) > 0.05 * dim(sce)[2], ]
    hvgs <- getTopHVGs(modelGeneVar(sce), n = as.integer(subs[i, "nhvg"]))
    sce <- sce[hvgs, ]
    
    fname <- paste0("Data/References/hdWGCNA/",
                    gsub(' ', '-', paste(ref, ct, id, sep = '_')), ".rds")
    saveRDS(sce, file = fname)
}
