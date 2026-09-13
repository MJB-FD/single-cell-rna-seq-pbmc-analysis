library(Seurat)

pbmc <- readRDS("results/variable_features/pbmc_variable_features.rds")

pbmc <- ScaleData(
    pbmc,
    features = VariableFeatures(pbmc)
)

print(pbmc)


dir.create(
    "results/scaling",
    recursive = TRUE,
    showWarnings = FALSE
)

saveRDS(
    pbmc,
    "results/scaling/pbmc_scaled.rds"
)
