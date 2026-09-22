library(Seurat)
library(ggplot2)

pbmc <- readRDS("results/clustering/pbmc_clustered.rds")

pbmc <- RunUMAP(
    pbmc,
    dims = 1:10
)

umap_plot <- DimPlot(
    pbmc,
    reduction = "umap",
    label = TRUE
)

print(umap_plot)

dir.create(
    "results/umap",
    recursive = TRUE,
    showWarnings = FALSE
)

ggsave(
    "results/umap/umap_clusters.png",
    plot = umap_plot,
    width = 8,
    height = 6
)

saveRDS(
    pbmc,
    "results/umap/pbmc_umap.rds"
)
