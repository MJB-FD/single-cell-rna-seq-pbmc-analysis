library(Seurat)
library(ggplot2)

pbmc <- readRDS("results/scaling/pbmc_scaled.rds")

pbmc <- RunPCA(
    pbmc,
    features = VariableFeatures(pbmc)
)

print(pbmc)

elbow_plot <- ElbowPlot(
    pbmc,
    ndims = 30
)

dir.create(
    "results/pca",
    recursive = TRUE,
    showWarnings = FALSE
)

ggsave(
    "results/pca/elbow_plot.png",
    plot = elbow_plot,
    width = 8,
    height = 6
)

saveRDS(
    pbmc,
    "results/pca/pbmc_pca.rds"
)
