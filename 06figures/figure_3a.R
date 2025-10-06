library(SingleCellExperiment)
library(scuttle)
library(dplyr)
library(scran)
library(scater)
library(fasano.franceschini.test)
set.seed(7)

sce <- readRDS("Data/References/Clustering/bailey24-clust.rds")
sce <- sce[, sce$cell_type == "Tregs"]
sce <- scuttle::logNormCounts(sce)
hvg <- scran::getTopHVGs(scran::modelGeneVar(sce))
hvg <- sample(hvg[1:100], 20)
N <- 1000
sce <- sce[hvg, sample(dim(sce)[2], N)]

# Project X and Y onto the PCA space of X
proj <- function(X, Y) {
    pr <- prcomp(X)
    rot <- pr$rotation
    Y <- Y[, colnames(X)]
    Yc <- scale(Y, center = colMeans(X), scale = FALSE)
    return(list(X_pc = pr$x, Y_pc = Yc %*% rot))
}

sim <- function(sce, family) {
    cell_types <- unique(sce$cell_type)
    sce_sims <- lapply(cell_types, function(ct) {
        sce_ct <- sce[, sce$cell_type == ct]
        if (family == "norm") {
            cop <- scCopula::fitGaussian(
                sce_ct,
                margins = "empirical",
                mle = FALSE,
                likelihood = FALSE
            )
        } else if (family == "vine") {
            cop <- scCopula::fitVine(
                sce_ct,
                margins = "empirical",
                cores = 4L
            )
        }
        Qx <- apply(as.matrix(counts(sce_ct)), 1, function(x) {
            function(p) { quantile(x, probs = p, type = 1) }
        })
        sce_sim <- scCopula::sampleCells(dim(sce_ct)[2], Qx, cop)
        sce_sim$cell_type <- ct
        return(sce_sim)
    })
    sce_sim <- do.call(cbind, sce_sims)
    return(t(counts(sce_sim)))
}

X <- t(counts(sce))
X_norm <- sim(sce, family = "norm")
X_vine <- sim(sce, family = "vine")

pr <- prcomp(X)
proj <- function(Y) {
    scale(Y[, colnames(X)], center = colMeans(X), scale = FALSE) %*% pr$rotation
}

pc_ref <- pr$x[, 1:2]
pc_norm <- proj(X_norm)[, 1:2]
pc_vine <- proj(X_vine)[, 1:2]

labels <- sapply(
    X = list("2" = pc_norm, "3" = pc_vine),
    FUN = function(x) {
        stat <- fasano.franceschini.test(pc_ref, x, nPermute = 0)$stat
        paste("FF =", sprintf("%.2f", stat / (2 * N)^(3/2)))
    }
)
labels["1"] <- "Reference"

df0 <- mutate(as.data.frame(pc_ref), group = "Reference")
df <- rbind(
    mutate(df0, plot = "1"),
    mutate(rbind(df0, mutate(as.data.frame(pc_norm),
                             group = "Gaussian")), plot = "2"),
    mutate(rbind(df0, mutate(as.data.frame(pc_vine),
                             group = "Vine")), plot = "3")
)
df$group <- factor(df$group, levels = c("Reference",
                                        "Gaussian",
                                        "Vine"))
pca_ex <- list(df = df, labels = labels)
colors <- setNames(c(ggsci::pal_npg()(10)[c(2, 3)], "#909090"),
                   c("Gaussian", "Vine", "Reference"))
pplt <- ggplot(pca_ex$df,
                 aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 1.0) +
    scale_color_manual(values = colors) +
    facet_wrap(~plot, scales = "free", nrow = 1,
               labeller = labeller(plot = pca_ex$labels)) +
    ylim(c(min(pca_ex$df$PC2), max(pca_ex$df$PC2))) +
    coord_cartesian(clip = "off") +
    guides(color = guide_legend(override.aes = list(size = 2.2))) +
    theme_bw() +
    theme(axis.title = element_text(size = 13),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 13),
          strip.background = element_blank(),
          strip.text = element_text(size = 13),
          legend.title = element_blank())
