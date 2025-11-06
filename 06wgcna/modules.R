library(Seurat)
library(SingleCellExperiment)
library(WGCNA)
library(hdWGCNA)
library(dplyr)
library(scCopula)
set.seed(0, kind = "L'Ecuyer-CMRG")

args <- commandArgs(trailingOnly = TRUE)
ncores <- as.integer(args[1])
family <- as.character(args[2])
ref_ix <- as.integer(args[3])
shuffle_ix <- as.integer(args[4])

refs <- list.files("Data/References/Subsets/")
refs <- refs[grepl("hdwgcna", refs)]
ref <- refs[ref_ix]

sce <- readRDS(paste0("Data/References/Subsets/", ref))
shuffles <- readRDS(paste0("Results/01shuffles/", ref))
ref <- unlist(strsplit(ref, ".rds"))
N <- floor(dim(sce)[2] / 2)

nsample <- 20L

sce_ref <- sce[, shuffles[shuffle_ix, (N + 1):(2 * N)]]
sce <- sce[, shuffles[shuffle_ix, 1:N]]

# Basic preprocessing
process <- function(obj) {
    fun <- function(x) {
        x %>%
            NormalizeData() %>%
            ScaleData() %>%
            RunPCA(seed = NULL)
    }
    return(suppressMessages(fun(obj)))
}

# Fit copula
if (family == "ind") {
    cop <- scCopula::fitIndep(sce)
} else if (family == "norm") {
    cop <- scCopula::fitGaussian(
        sce = sce,
        margins = "empirical",
        mle = FALSE,
        jitter = FALSE,
        likelihood = FALSE
    )
} else if (family == "norm_jitter") {
    cop <- scCopula::fitGaussian(
        sce = sce,
        margins = "empirical",
        mle = FALSE,
        jitter = TRUE,
        likelihood = FALSE
    )
} else if (family == "vine_jitter") {
    cop <- scCopula::fitVine(
        sce = sce,
        margins = "empirical",
        jitter = TRUE,
        cores = ncores
    )
}

# Preprocess reference dataset
seurat_ref <- as.Seurat(sce_ref)
VariableFeatures(seurat_ref) <- rownames(seurat_ref)
seurat_ref <- process(seurat_ref)

# Generate and preprocess simulated datasets
sims <- list()
for (i in seq_len(nsample)) {
    Qx <- apply(as.matrix(counts(sce)), 1, function(x) {
        function(p) { quantile(x, probs = p, type = 1) }
    })
    sce_sim <- scCopula::sampleCells(
        N = dim(sce)[2],
        Qx = Qx,
        copula = cop
    )
    sce_sim$cell_type <- sce$cell_type
    
    seurat_sim <- as.Seurat(sce_sim)
    DefaultAssay(seurat_sim) <- "originalexp"
    seurat_sim[["RNA"]] <- seurat_sim[["originalexp"]]
    DefaultAssay(seurat_sim) <- "RNA"
    seurat_sim[["originalexp"]] <- NULL
    VariableFeatures(seurat_sim) <- rownames(seurat_sim)
    
    sims[[i]] <- process(seurat_sim)
}

# Get modules for reference dataset
enableWGCNAThreads(nThreads = ncores)
seurat_ref <- SetupForWGCNA(seurat_ref, wgcna_name = "test")
seurat_ref <- MetacellsByGroups(
    seurat_obj = seurat_ref,
    group.by = "cell_type",
    reduction = "pca",
    ident.group = "cell_type",
    target_metacells = 500
)
seurat_ref <- NormalizeMetacells(seurat_ref)
seurat_ref <- SetDatExpr(seurat_ref, assay = "RNA", layer = "data")
seurat_ref <- TestSoftPowers(seurat_ref)
seurat_ref <- ConstructNetwork(seurat_ref, overwrite_tom = TRUE)

# Module preservation
res <- c()
for (i in seq_len(nsample)) {
    seurat_sim <- sims[[i]]
    
    # Project modules from reference dataset onto simulated dataset
    seurat_sim <- ProjectModules(
        seurat_obj = seurat_sim,
        seurat_ref = seurat_ref,
        wgcna_name = "test",
        wgcna_name_proj = "projected",
        assay = "RNA"
    )
    
    # Construct metacells for simulated dataset
    seurat_sim <- MetacellsByGroups(
        seurat_obj = seurat_sim,
        group.by = "cell_type",
        reduction = "pca",
        ident.group = "cell_type",
        target_metacells = 500
    )
    seurat_sim <- NormalizeMetacells(seurat_sim)
    
    # Run module preservation analysis
    seurat_sim <- SetDatExpr(seurat_sim, assay = "RNA", layer = "data")

    res[i] <- tryCatch({
        seurat_sim <- ModulePreservation(
            seurat_obj = seurat_sim,
            seurat_ref = seurat_ref,
            name = "comp",
            parallel = TRUE,
            n_permutations = 100
        )
        p <- GetModulePreservation(seurat_sim, "comp")$Z
        # Extract Z statistics
        mean(p[!rownames(p) %in% c("gold", "grey"), "Zsummary.pres"])
    },
    error = function(e) {
        return(NA)
    })
}
res <- data.frame(
    Z = res,
    family = family,
    ref = ref,
    trial = shuffle_ix
)
fname <- paste0("Results/06wgcna/", family, "_", ref_ix, "_", shuffle_ix, ".rds")
saveRDS(res, file = fname)
