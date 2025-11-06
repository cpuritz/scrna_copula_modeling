library(ggplot2)
library(SingleCellExperiment)
library(scales)
library(scran)
set.seed(2)

setwd("~/Documents/Graduate School/Copula Paper")

files <- list.files("Data/References/SingleCellExperiment")
ngenes <- floor(10^runif(length(files), min = log10(500), max = log10(1000)))
refs <- sapply(files, function(x) {
    paste(unlist(strsplit(gsub(".rds", "", x), '-'))[1:3], collapse = '-')
})
names(refs) <- NULL
refs <- data.frame(ref = refs, ngene = ngenes)

dist <- NULL
for (i in seq_len(dim(refs)[1])) {
    ref <- refs[i, "ref"]
    sce <- readRDS(paste0("Data/References/SingleCellExperiment/", ref, ".rds"))
    # Remove genes expressed in less than 5% of cells
    sce <- sce[rowSums(counts(sce) > 0) >= 0.05 * dim(sce)[2], ]
    hvgs <- getTopHVGs(modelGeneVar(sce))
    if (length(hvgs) < refs[i, "ngene"]) {
        refs[i, "ngene"] <- length(hvgs)
    }
    hvgs <- hvgs[seq_len(refs[i, "ngene"])]
    sce <- sce[hvgs, ]
    assertthat::assert_that(dim(sce)[1] == refs[i, "ngene"])
    saveRDS(sce, paste0("Data/References/Subsets/", ref, "-hdwgcna", ".rds"))
    dist <- rbind(dist, c(ref = ref, ngene = dim(sce)[1], ncell = dim(sce)[2]))
}
dist <- as.data.frame(dist)
dist$ngene <- as.integer(dist$ngene)
dist$ncell <- as.integer(dist$ncell)

pplt <- ggplot(data = dist, aes(x = ncell, y = ngene)) +
    geom_point(size = 2.2) +
    xlab("Cells") +
    ylab("Genes") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 13),
          axis.text = element_text(color = "black", size = 13),
          axis.title = element_text(size = 13))
saveRDS(pplt, "Figures/fig_s1b.rds")
