suppressMessages(library(Seurat))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(WGCNA))
suppressMessages(library(hdWGCNA))
suppressMessages(library(dplyr))
suppressMessages(library(mclust))
set.seed(0, kind = "L'Ecuyer-CMRG")

setwd("~/copgen")
enableWGCNAThreads(nThreads = 20)

refs <- list.files("Data/References/Subsets/")
ref <- refs[grepl("hdwgcna", refs)][3]

sce <- readRDS(paste0("Data/References/Subsets/", ref))
shuffles <- readRDS(paste0("Results/01shuffles/", ref))
ref <- unlist(strsplit(ref, ".rds"))
N <- floor(dim(sce)[2] / 2)

sce_ref <- sce[, shuffles[1, (N + 1):(2 * N)]]
sce <- sce[, shuffles[1, 1:N]]

process <- function(obj) {
    fun <- function(x) {
        x %>%
            NormalizeData() %>%
            ScaleData() %>%
            RunPCA(seed = NULL)
    }
    return(suppressMessages(fun(obj)))
}

seurat_ref <- as.Seurat(sce_ref)
VariableFeatures(seurat_ref) <- rownames(seurat_ref)
seurat_ref <- process(seurat_ref)
seurat_ref <- SetupForWGCNA(seurat_ref, wgcna_name = "test")
seurat_ref <- MetacellsByGroups(
    seurat_obj = seurat_ref,
    group.by = "cell_type",
    reduction = "pca",
    ident.group = "cell_type",
    target_metacells = 500
)
seurat_ref <- suppressMessages(NormalizeMetacells(seurat_ref))
seurat_ref <- SetDatExpr(seurat_ref, assay = "RNA", layer = "data")
seurat_ref <- suppressMessages(TestSoftPowers(seurat_ref))
seurat_ref <- suppressMessages(ConstructNetwork(seurat_ref, overwrite_tom = TRUE))


sce_sim <- sce[, sample(dim(sce)[2], replace = TRUE)]
colnames(sce_sim) <- paste0("cell_", seq_len(dim(sce_sim)[2])) 
seurat_sim <- as.Seurat(sce_sim)
VariableFeatures(seurat_sim) <- rownames(seurat_sim)
seurat_sim <- process(seurat_sim)
seurat_sim <- ProjectModules(
    seurat_obj = seurat_sim,
    seurat_ref = seurat_ref,
    wgcna_name = "test",
    wgcna_name_proj = "projected",
    assay = "RNA"
)
seurat_sim <- suppressMessages(MetacellsByGroups(
    seurat_obj = seurat_sim,
    group.by = "cell_type",
    reduction = "pca",
    ident.group = "cell_type",
    target_metacells = 500
))
seurat_sim <- suppressMessages(NormalizeMetacells(seurat_sim))
seurat_sim <- SetDatExpr(seurat_sim, assay = "RNA", layer = "data")
seurat_sim <- suppressMessages(ModulePreservation(
    seurat_obj = seurat_sim,
    seurat_ref = seurat_ref,
    name = "comp",
    parallel = TRUE,
    n_permutations = 100
))
p <- GetModulePreservation(seurat_sim, "comp")$Z
p <- p[!rownames(p) %in% c("gold", "grey"), ]
z <- data.frame(module = rownames(p), Z = p$Zsummary.pres, size = p$moduleSize)

# Figure 4a
PlotDendrogram(seurat_ref)

# Figure 4b
colors <- c("blue" = "#2A6EBBFF", "brown" = "#CD202CFF",
            "turquoise" = "#00B2A9FF", "yellow" = "#F0AB00FF")
pplt <- ggplot(z, aes(x = size, y = Z, color = module)) +
    geom_point(size = 4) +
    scale_color_manual(values = colors) +
    xlab("Module size") +
    ylab(expression(Z[summary])) +
    geom_hline(yintercept = mean(z$Z), color = "#767676") +
    scale_x_log10() +
    annotate("text", x = 273, y = 1.02 * mean(z$Z), label = "bar(Z)[summary]",
             parse = TRUE, size = 4.2) +
    theme_bw() +
    theme(
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 13, color = "black"),
        legend.position = "none"
    )
