# 11_cluster_annotation.R
#
# Purpose:
# Assign biological cell-type labels to Seurat clusters
# based on canonical markers and differential expression results.

library(Seurat)

# Load UMAP-processed Seurat object
pbmc <- readRDS("results/umap/pbmc_umap.rds")

# Map cluster numbers to biological identities
cluster_labels <- c(
  "0" = "Naive T cells",
  "1" = "Classical monocytes",
  "2" = "IL7R+ CD4 T cells",
  "3" = "B cells",
  "4" = "Cytotoxic CD8 T cells",
  "5" = "FCGR3A+ monocytes",
  "6" = "NK cells",
  "7" = "CD1C+ dendritic cells",
  "8" = "Platelets"
)

# Rename cluster identities
pbmc <- RenameIdents(pbmc, cluster_labels)

# Store cell-type labels in metadata
pbmc$cell_type <- Idents(pbmc)

# Verify annotation
print(table(Idents(pbmc)))

# Save annotated object
dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)

saveRDS(
  pbmc,
  "results/annotation/pbmc_annotated.rds"
)
