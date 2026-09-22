library(Seurat)

pbmc <- readRDS("results/scaling/pbmc_scaled.rds")

pbmc <- RunPCA(
    pbmc,
    features = VariableFeatures(pbmc)
)

pbmc <- FindNeighbors(
    pbmc,
    dims = 1:10
)

print(pbmc)

dir.create(
    "results/neighbors",
    recursive = TRUE,
    showWarnings = FALSE
)

saveRDS(
    pbmc,
    "results/neighbors/pbmc_neighbors.rds"
)
