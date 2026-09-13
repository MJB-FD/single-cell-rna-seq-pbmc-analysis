library(Seurat)
library(ggplot2)

# Load the normalized Seurat object
pbmc <- readRDS("results/normalized/pbmc_normalized.rds")

# Identify highly variable genes
pbmc <- FindVariableFeatures(
	pbmc,
	selection.method = "vst",
	nfeatures = 2000
)


# Display the number of variable features selected
cat("Number of variable features:", length(VariableFeatures(pbmc)), "\n")


# Dispaly the first variablw features
print(head(VariableFeatures(pbmc), 20))

# Create a plot showing highly variable features
variable_plot <- VariableFeaturePlot(pbmc)

# Save the plot
dir.create("results/variable_features",
           recursive = TRUE,
           showWarnings = FALSE)

ggsave(
    "results/variable_features/highly_variable_genes.png",
    plot = variable_plot,
    width = 10,
    height = 6
)

dir.create(
    "results/variable_features",
    recursive = TRUE,
    showWarnings = FALSE
)

saveRDS(
    pbmc,
    "results/variable_features/pbmc_variable_features.rds"
)
