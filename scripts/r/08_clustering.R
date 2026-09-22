library(Seurat)

pbmc <- readRDS("results/neighbors/pbmc_neighbors.rds")

pbmc <- FindClusters(
    pbmc,
    resolution = 0.5
)

print(pbmc)

dir.create(
    "results/clustering",
    recursive = TRUE,
    showWarnings = FALSE
)

saveRDS(
    pbmc,
    "results/clustering/pbmc_clustered.rds"
)
