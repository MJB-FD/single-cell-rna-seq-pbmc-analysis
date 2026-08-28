# Load Seurat
library(Seurat)

# Location of the 10x Genomics count matrix
data_dir <- "data/raw/filtered_gene_bc_matrices/hg19"

# Read the 10x count matrix
pbmc_counts <- Read10X(data.dir = data_dir)

# Inspect the matrix
print(class(pbmc_counts))
print(dim(pbmc_counts))
print(pbmc_counts[1:5, 1:5])

# Create a Seurat object
pbmc <- CreateSeuratObject(
    counts = pbmc_counts,
    project = "PBMC3K"
)

# Inspect the Seurat object
print(pbmc)
print(dim(pbmc))
