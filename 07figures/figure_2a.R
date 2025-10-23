library(SingleCellExperiment)
library(scuttle)
library(dplyr)
library(scran)
library(scater)
library(fasano.franceschini.test)
set.seed(10)
setwd("~/Documents/Graduate School/Copula Paper")

sce <- readRDS("Data/References/Clustering/bailey24-clust.rds")
sce <- sce[, sce$cell_type == "Tregs"]
sce <- scuttle::logNormCounts(sce)
hvg <- scran::getTopHVGs(scran::modelGeneVar(sce), n = 200)
hvg <- sample(hvg, 30)
N <- 2000
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
        } else if (family == "ind") {
            cop <- scCopula::fitIndep(sce_ct)
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
X_ind <- sim(sce, family = "ind")
X_norm <- sim(sce, family = "norm")

pr <- prcomp(X)
proj <- function(Y) {
    scale(Y[, colnames(X)], center = colMeans(X), scale = FALSE) %*% pr$rotation
}

pc_ref <- pr$x[, 1:2]
pc_ind <- proj(X_ind)[, 1:2]
pc_norm <- proj(X_norm)[, 1:2]

labels <- sapply(
    X = list("2" = pc_ind, "3" = pc_norm),
    FUN = function(x) {
        p <- fasano.franceschini.test(pc_ref, x, threads = 4, seed = 0)$p.value
        paste("p =", sprintf("%.2f", as.numeric(p)))
    }
)
labels["1"] <- "Reference"

df0 <- mutate(as.data.frame(pc_ref), group = "Reference")
df <- rbind(
    mutate(df0, plot = "1"),
    mutate(rbind(df0, mutate(as.data.frame(pc_ind),
                             group = "Independence")), plot = "2"),
    mutate(rbind(df0, mutate(as.data.frame(pc_norm),
                             group = "Gaussian")), plot = "3")
)
df$group <- factor(df$group, levels = c("Reference",
                                        "Independence",
                                        "Gaussian"))
pca_ex <- list(df = df, labels = labels)
colors <- setNames(c(ggsci::pal_npg()(10)[c(7, 5)], "#767676"),
                   c("Independence", "Gaussian", "Reference"))
pplt <- ggplot(pca_ex$df,
                 aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 1.3) +
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
saveRDS(pca_ex, file = "Results/04pca/pca_ex.rds")
