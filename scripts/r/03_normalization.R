library(Seurat)

# Load the QC-filtered Seurat object
pbmc <- readRDS("results/qc/pbmc_filtered.rds")

# Normalize gene expression
pbmc <- NormalizeData(
    pbmc,
    normalization.method = "LogNormalize",
    scale.factor = 10000
)

# Check the dimensions of the normalized data
print(dim(pbmc))

# Check that the RNA assay now contains a normalized data layer
print(pbmc[["RNA"]])

saveRDS(
    pbmc,
    file = "results/normalized/pbmc_normalized.rds"
)
