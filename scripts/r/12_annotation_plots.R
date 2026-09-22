# 12_annotation_plots.R
#
# Purpose:
# Create final visualization of annotated PBMC cell types.

library(Seurat)
library(ggplot2)

# Load annotated Seurat object
pbmc <- readRDS("results/annotation/pbmc_annotated.rds")

# Create annotated UMAP
annotated_umap <- DimPlot(
  pbmc,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) +
  NoLegend()

# Save plot
ggsave(
  filename = "results/annotation/annotated_umap.png",
  plot = annotated_umap,
  width = 10,
  height = 7,
  dpi = 300
)
