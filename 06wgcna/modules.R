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

message("Running hdWGCNA for reference ", ref, ", shuffle ", shuffle_ix, ", for family ", family, "\n\n\n")

sce <- readRDS(paste0("Data/References/Subsets/", ref))
shuffles <- readRDS(paste0("Results/01shuffles/", ref))
ref <- unlist(strsplit(ref, ".rds"))
N <- floor(dim(sce)[2] / 2)

nsample <- 20L

sce_ref <- sce[, shuffles[shuffle_ix, (N + 1):(2 * N)]]
sce <- sce[, shuffles[shuffle_ix, 1:N]]

###############

setwd(paste0("Scripts/06wgcna/", family, "/", ref_ix, "_", shuffle_ix))

# Preprocess reference dataset
message("Preprocessing reference dataset")
seurat_ref <- as.Seurat(sce_ref)
VariableFeatures(seurat_ref) <- rownames(seurat_ref)
seurat_ref <- NormalizeData(seurat_ref)
seurat_ref <- ScaleData(seurat_ref)
seurat_ref <- RunPCA(seurat_ref, seed = NULL)

# Get modules for reference dataset
message("Computing reference modules")
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

message("Done with reference modules")

###############

# Fit model
if (family == "ind") {
    model <- scCopula::fitIndep(sce)
    Qx <- apply(as.matrix(counts(sce)), 1, par_quantile, simplify = FALSE)
} else if (family == "norm") {
    model <- scCopula::fitGaussian(
        sce = sce,
        margins = "par",
        estimate = "sample"
    )
    Qx <- apply(as.matrix(counts(sce)), 1, par_quantile, simplify = FALSE)
} else if (family == "norm_jitter") {
    model <- scCopula::fitGaussian(
        sce = sce,
        margins = "par",
        estimate = "jitter"
    )
    Qx <- apply(as.matrix(counts(sce)), 1, par_quantile, simplify = FALSE)
} else if (family == "vine_jitter") {
    model <- scCopula::fitVine(
        sce = sce,
        margins = "par",
        jitter = TRUE,
        cores = ncores
    )
    Qx <- apply(as.matrix(counts(sce)), 1, par_quantile, simplify = FALSE)
} else if (family == "zinbwave") {
    model <- zinbwave::zinbFit(
        sce,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )
} else if (family == "sparsim") {
    sink(tempfile())
    model <- SPARSim::SPARSim_estimate_parameter_from_data(
        raw_data = as.matrix(counts(sce)),
        norm_data = as.matrix(logcounts(sce)),
        conditions = list(cond = seq_len(dim(sce)[2]))
    )
    sink()
}

# Generate and preprocess simulated datasets
sims <- list()
for (i in seq_len(nsample)) {
    message("Simulating ", i, "/", nsample)
    if (family == "zinbwave") {
        X <- zinbwave::zinbSim(model)$counts
        sce_sim <- SingleCellExperiment(
            list(counts = X, logcounts = log1p(X)),
            colData = colData(sce)
        )
        rownames(sce_sim) <- rownames(sce)
    } else if (family == "sparsim") {
        sink(tempfile())
        sim <- SPARSim::SPARSim_simulation(model)
        X <- sim$count_matrix
        sce_sim <- SingleCellExperiment(
            list(counts = X, logcounts = log1p(X)),
            colData = colData(sce)
        )
        rownames(sce_sim) <- rownames(sce)
        sink()
    } else {
        sce_sim <- scCopula::sampleCells(
            N = dim(sce)[2],
            Qx = Qx,
            copula = model
        )
    }
    sce_sim$cell_type <- sce$cell_type
    
    seurat_sim <- as.Seurat(sce_sim)
    DefaultAssay(seurat_sim) <- "originalexp"
    seurat_sim[["RNA"]] <- seurat_sim[["originalexp"]]
    DefaultAssay(seurat_sim) <- "RNA"
    seurat_sim[["originalexp"]] <- NULL
    VariableFeatures(seurat_sim) <- rownames(seurat_sim)
    
    seurat_sim <- NormalizeData(seurat_sim)
    seurat_sim <- ScaleData(seurat_sim)
    seurat_sim <- RunPCA(seurat_sim, seed = NULL)
    
    sims[[i]] <- seurat_sim
}
message("\n\n\n!!!!! DONE SIMULATING !!!!!\n\n\n")

# Module preservation
res <- c()
for (i in seq_len(nsample)) {
    message("\n\n\n!!!!!!! MODULES ", i, "/", nsample, " !!!!!!!\n\n\n")
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
saveRDS(
    object = data.frame(Z = res, family = family, ref = ref, trial = shuffle_ix),
    file = paste0("Results/06wgcna/", family, "_", ref_ix, "_", shuffle_ix, ".rds")
)
