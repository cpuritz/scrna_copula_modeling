suppressMessages(library(Seurat))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(WGCNA))
suppressMessages(library(hdWGCNA))
suppressMessages(library(dplyr))
suppressMessages(library(mclust))
set.seed(1, kind = "L'Ecuyer-CMRG")

setwd("~/copgen")

ncores <- as.integer(commandArgs(trailingOnly = TRUE)[1])
family <- as.character(commandArgs(trailingOnly = TRUE)[2])
message("Running hdWGCNA analysis for ", family, " with ", ncores, " cores.")

sce <- readRDS("Data/References/Subsets/bailey24-rpra30-cd8t.rds")

runv <- "v1"
if (runv == "v1") {
    # Model all cells as one cell type
    sce$cell_type <- "CD8 T cell"
    mu_formula <- corr_by <- "1"
} else {
    # Use original clusters as individual cell types
    mu_formula <- corr_by <- "cell_type"
}

process <- function(obj) {
    obj %>%
        NormalizeData() %>%
        ScaleData() %>%
        RunPCA(seed = NULL)
}

get_modules <- function(obj) {
    obj <- SetupForWGCNA(obj, wgcna_name = "test")
    obj <- MetacellsByGroups(
        seurat_obj = obj,
        group.by = "cell_type",
        reduction = "pca",
        ident.group = "cell_type"
    )
    obj <- NormalizeMetacells(obj)
    obj <- SetDatExpr(obj, assay = "RNA", layer = "data")
    obj <- TestSoftPowers(obj)
    obj <- ConstructNetwork(obj, overwrite_tom = TRUE)
    
    modules <- GetModules(obj)
    missing_genes <- rownames(obj)[!rownames(obj) %in% modules$gene_name]
    if (length(missing_genes) > 0) {
        missing_genes_df <- do.call(rbind, lapply(missing_genes, function(x) {
            c(x, "grey", "grey")
        }))
        missing_genes_df <- data.frame(missing_genes_df)
        rownames(missing_genes_df) <- missing_genes
        colnames(missing_genes_df) <- colnames(modules)
        modules <- rbind(modules, missing_genes_df)
    }
    return(modules)
}

# Construct model of reference dataset
if (family != "boot") {
    message("Constructing input data")
    input_data <- scDesign3::construct_data(
        sce = sce,
        assay_use = "counts",
        celltype = "cell_type",
        other_covariates = NULL,
        pseudotime = NULL,
        spatial = NULL,
        corr_by = corr_by
    )
    
    message("Fitting margins")
    margins <- scDesign3::fit_marginal(
        mu_formula = mu_formula,
        sigma_formula = "1",
        n_cores = ncores,
        usebam = FALSE,
        data = input_data,
        family_use = "nb"
    )
    
    message("Extracting parameters")
    para_list <- scDesign3::extract_para(
        sce = sce,
        assay_use = "counts",
        marginal_list = margins,
        family_use = "nb",
        new_covariate = input_data$newCovariate,
        n_cores = ncores,
        data = input_data$dat
    )
    
    cop <- scDesign3::fit_copula(
        sce = sce,
        assay_use = "counts",
        input_data = input_data$dat,
        marginal_list = margins,
        family_use = "nb",
        copula = ifelse(grepl("norm", family), "gaussian", "vine"),
        DT = grepl("jitter", family),
        family_set = c("indep", "gaussian", "frank", "clayton", "joe", "gumbel"),
        n_cores = ncores
    )$copula_list
}

# Process reference dataset
message("Processing reference dataset")
seurat_ref <- as.Seurat(sce)
VariableFeatures(seurat_ref) <- rownames(seurat_ref)
seurat_ref <- process(seurat_ref)

nrun <- 50L
sims <- list()
for (i in seq_len(nrun)) {
    message("Simulating ", i, "/", nrun)
    if (family == "boot") {
        uct <- unique(sce$cell_type)
        if (length(uct) > 1) {
            boot_ix <- do.call(c, sapply(uct, function(x) {
                sample(which(sce$cell_type == x), replace = TRUE)
            }))
        } else {
            boot_ix <- sample(dim(sce)[2], replace = TRUE)
        }
        sce_boot <- sce[, boot_ix]
        colnames(sce_boot) <- paste0("cell_", seq_len(dim(sce_boot)[2]))
        
        seurat_boot <- as.Seurat(sce_boot)
        VariableFeatures(seurat_boot) <- rownames(seurat_boot)
        seurat_boot <- process(seurat_boot)
        sims[[i]] <- seurat_boot
    } else {
        new_counts <- scDesign3::simu_new(
            sce = sce,
            assay_use = "counts",
            mean_mat = para_list$mean_mat,
            sigma_mat = para_list$sigma_mat,
            zero_mat = para_list$zero_mat,
            copula_list = cop,
            n_cores = ncores,
            family_use = "nb",
            input_data = input_data$dat,
            new_covariate = input_data$newCovariate,
            filtered_gene = NULL,
            important_feature = rep(TRUE, dim(sce)[1])
        )
        sce_sim <- SingleCellExperiment(list(counts = new_counts),
                                        colData = input_data$newCovariate)
        logcounts(sce_sim) <- log1p(counts(sce_sim))

        seurat_sim <- as.Seurat(sce_sim)
        DefaultAssay(seurat_sim) <- "originalexp"
        seurat_sim[["RNA"]] <- seurat_sim[["originalexp"]]
        DefaultAssay(seurat_sim) <- "RNA"
        seurat_sim[["originalexp"]] <- NULL
        
        VariableFeatures(seurat_sim) <- rownames(seurat_sim)
        seurat_sim <- process(seurat_sim)
        sims[[i]] <- seurat_sim
    }
}
message("\n\n\n!!!!! DONE SIMULATING !!!!!\n\n\n")

setwd(paste0("Scripts/06wgcna/", family))

# Get modules for reference dataset
message("Computing reference modules")
enableWGCNAThreads(nThreads = ncores)
mod_ref <- get_modules(seurat_ref)

res <- c()
for (i in seq_len(nrun)) {
    message("\n\n\n!!!!!!! MODULES ", i, "/", nrun, " !!!!!!!\n\n\n")
    # Get modules for simulated dataset
    res[i] <- tryCatch({
        mod_sim <- get_modules(sims[[i]])
        mod_sim <- mod_sim[rownames(mod_ref), ]
        mclust::adjustedRandIndex(mod_ref$module, mod_sim$module)
    },
    error = function(e) {
        return(NA)
    })
}
saveRDS(
    object = data.frame(ARI = res),
    file = paste0("~/copgen/Results/06wgcna/mod_", family, "_", vrun, ".rds")
)

message("Done!")
