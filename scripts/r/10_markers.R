library(Seurat)

pbmc <- readRDS("results/umap/pbmc_umap.rds")

markers <- FindAllMarkers(
    pbmc,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
)

dir.create(
    "results/markers",
    recursive = TRUE,
    showWarnings = FALSE
)

write.csv(
    markers,
    "results/markers/all_cluster_markers.csv",
    row.names = FALSE
)

print(markers)
