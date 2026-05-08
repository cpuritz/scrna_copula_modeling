library(graphite)
library(igraph)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
set.seed(0)

setwd("Data")

# Download KEGG pathways
hsa_pw <- pathways("hsapiens", "kegg")

# Extract pathways of interest
pw_list <- list(
    hsa_breast_cancer = hsa_pw[["Breast cancer"]],
    hsa_cardiac = hsa_pw[["Cardiac muscle contraction"]],
    hsa_estrogen = hsa_pw[["Estrogen signaling pathway"]],
    hsa_insulin = hsa_pw[["Insulin secretion"]]
)

# Convert to igraph objects
pw_list <- lapply(pw_list, function(x) {
    graph_from_graphnel(pathwayGraph(x, "proteins"))
})

# Convert to undirected graphs
pw_list <- lapply(pw_list, as_undirected)

downsample <- function(x, n) {
    cl <- cluster_leiden(x, resolution = 5e-3)
    dcl <- abs(table(cl$membership) - n)
    return(cl$names[cl$membership == which.min(dcl)])
}

# Downsample with specified target sizes
hsa_breast_cancer <- downsample(pw_list$hsa_breast_cancer, 50)
hsa_estrogen <- downsample(pw_list$hsa_estrogen, 60)

# Already small enough
hsa_insulin <- names(V(pw_list$hsa_insulin))
hsa_cardiac <- names(V(pw_list$hsa_cardiac))

# Map to Ensembl identifiers
hsa_map <- AnnotationDbi::select(
    x = org.Hs.eg.db,
    keys = AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID"),
    columns = "ENSEMBL",
    keytype = "ENTREZID"
)
hsa_map$ENTREZID <- paste0("ENTREZID:", hsa_map$ENTREZID)

hsa_breast_cancer <- hsa_map$ENSEMBL[hsa_map$ENTREZID %in% hsa_breast_cancer]
hsa_estrogen <- hsa_map$ENSEMBL[hsa_map$ENTREZID %in% hsa_estrogen]
hsa_cardiac <- hsa_map$ENSEMBL[hsa_map$ENTREZID %in% hsa_cardiac]
hsa_insulin <- hsa_map$ENSEMBL[hsa_map$ENTREZID %in% hsa_insulin]

saveRDS(hsa_breast_cancer, "hsa05224.rds")
saveRDS(hsa_estrogen, "hsa04915.rds")
saveRDS(hsa_cardiac, "hsa04260.rds")
saveRDS(hsa_insulin, "hsa04911.rds")
