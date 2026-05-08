library(ggplot2)
library(SingleCellExperiment)
library(scales)
library(scran)
set.seed(0)

setwd("~/Documents/Graduate School/Copula Paper")

# Read in gene sets
hsa04260 <- readRDS("Data/genesets/hsa04260.rds")
hsa04915 <- readRDS("Data/genesets/hsa04915.rds")
hsa05224 <- readRDS("Data/genesets/hsa05224.rds")
hsa04911 <- readRDS("Data/genesets/hsa04911.rds")

res <- list(
    c("bailey24-rpra04-moam", 17, "HVG"),
    c("bailey24-rpra09-tram", 12, "HVG"),
    c("bailey24-rpra24-tram", 15, "HVG"),
    c("bailey24-rpra29-dc2",  35, "HVG"),
    c("bailey24-rpra20-cd8t", 40, "HVG"),
    c("bailey24-rpra30-cd4t", 32, "HVG"),
    c("reed24-donor6-myoepithelial",  30, "HVG"),
    c("reed24-donor18-myoepithelial", 36, "HVG"),
    c("reed24-donor19-myoepithelial", 45, "HVG"),
    c("reed24-donor26-myoepithelial", 50, "HVG"),
    c("reed24-donor6-myoepithelial",  NA, "hsa05224"),
    c("reed24-donor18-myoepithelial", NA, "hsa05224"),
    c("reed24-donor19-myoepithelial", NA, "hsa05224"),
    c("reed24-donor26-myoepithelial", NA, "hsa05224"),
    c("jones24-donor3-stromal", 25, "HVG"),
    c("jones24-donor4-stromal", 20, "HVG"),
    c("jones24-donor5-stromal", 15, "HVG"),
    c("jones24-donor3-stromal", NA, "hsa04915"),
    c("jones24-donor4-stromal", NA, "hsa04915"),
    c("jones24-donor5-stromal", NA, "hsa04915"),
    c("lukassen20-9JQK55ng-monocyte", 32, "HVG"),
    c("lukassen20-A9LCTZng-monocyte", 10, "HVG"),
    c("lukassen20-QZY9VQng-monocyte", 14, "HVG"),
    c("orozco20-109373-mueller", 20, "HVG"),
    c("orozco20-109829-mueller", 12, "HVG"),
    c("orozco20-110814-mueller", 42, "HVG"),
    c("orozco20-120974-mueller", 10, "HVG"),
    c("litvinukova20-d3-myocyte", NA, "hsa04260"),
    c("litvinukova20-d6-myocyte", NA, "hsa04260"),
    c("litvinukova20-h6-myocyte", NA, "hsa04260"),
    c("litvinukova20-h7-myocyte", NA, "hsa04260"),
    c("king21-BCP5-naiveB", 16, "HVG"),
    c("king21-BCP6-naiveB", 20, "HVG"),
    c("king21-BCP8-naiveB", 24, "HVG"),
    c("james20-290b-memB", 10, "HVG"),
    c("james20-302c-memB", 36, "HVG"),
    c("james20-390c-memB", 48, "HVG"),
    c("james20-290b-iga",  15, "HVG"),
    c("james20-298c-iga",  26, "HVG"),
    c("gu24-pooled-mp", 10, "HVG"),
    c("gu24-pooled-dc", 46, "HVG"),
    c("gu24-pooled-plasma", 28, "HVG"),
    c("tritschler22-pig1-panB", 10, "HVG"),
    c("tritschler22-pig2-panA", 30, "HVG"),
    c("tritschler22-pig2-panB", 20, "HVG"),
    c("tritschler22-pig1-panB", NA, "hsa04911"),
    c("tritschler22-pig2-panA", NA, "hsa04911"),
    c("tritschler22-pig2-panB", NA, "hsa04911")
)
res <- as.data.frame(do.call(rbind, res))
colnames(res) <- c("subset", "ngene", "type")
res$ngene <- as.numeric(res$ngene)
res$ncell <- rep(NA, dim(res)[1])

dist <- NULL
for (i in seq(dim(res)[1])) {
    x <- res[i, ]
    sce <- readRDS(paste0("Data/References/SingleCellExperiment/", x$subset, ".rds"))
    if (x$type == "HVG") {
        # Compute HVGs
        hvgs <- getTopHVGs(modelGeneVar(sce))
        # Only keep HVGs expressed in at least 2% of cells
        genes_keep <- names(which(rowMeans(counts(sce) > 0) > 0.02))
        hvgs <- intersect(hvgs, genes_keep)
        # Subset to top HVGs
        genes <- hvgs[seq_len(x$ngene)]
    } else {
        # Only keep genes expressed in at least 2% of cells
        genes_keep <- names(which(rowMeans(counts(sce) > 0) > 0.02))
        genes <- intersect(get(x$type), genes_keep)
    }
    sce <- sce[genes, ]
    
    # Remove cells not expressing any gene
    sce <- sce[, colSums(counts(sce)) > 0]
    
    res$ngene[i] <- dim(sce)[1]
    res$ncell[i] <- dim(sce)[2]
    
    ofname <- paste0("Data/References/Subsets/", x$subset, "-", x$type, ".rds")
    saveRDS(sce, ofname)
}