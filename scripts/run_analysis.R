# ============================================================
# PBMC 3K Single-Cell RNA-seq Analysis
# Complete reproducible workflow
#
# Input:
#   10x Genomics filtered gene-barcode matrix
#
# Output:
#   QC figures
#   Variable-feature plot
#   PCA elbow plot
#   UMAP plots
#   Cluster marker table
#   Marker DotPlot
#   Final annotated Seurat object
#
# Unlike the individual teaching scripts in scripts/r/,
# this master workflow keeps the Seurat object in memory
# and saves only the final processed object.
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(Seurat)
library(ggplot2)


# ------------------------------------------------------------
# 2. Create output directories
# ------------------------------------------------------------

dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)
dir.create("results/variable_features", recursive = TRUE, showWarnings = FALSE)
dir.create("results/pca", recursive = TRUE, showWarnings = FALSE)
dir.create("results/umap", recursive = TRUE, showWarnings = FALSE)
dir.create("results/markers", recursive = TRUE, showWarnings = FALSE)
dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 3. Read the 10x Genomics count matrix
# ------------------------------------------------------------

cat("\n[1/12] Reading 10x Genomics count matrix...\n")

data_dir <- "data/raw/filtered_gene_bc_matrices/hg19"

pbmc_counts <- Read10X(
  data.dir = data_dir
)

cat(
  "Genes:", nrow(pbmc_counts),
  "\nCells:", ncol(pbmc_counts),
  "\n"
)


# ------------------------------------------------------------
# 4. Create Seurat object
# ------------------------------------------------------------

cat("\n[2/12] Creating Seurat object...\n")

pbmc <- CreateSeuratObject(
  counts = pbmc_counts,
  project = "PBMC3K"
)

cat("Starting cells:", ncol(pbmc), "\n")


# ------------------------------------------------------------
# 5. Quality control
# ------------------------------------------------------------

cat("\n[3/12] Performing quality control...\n")

pbmc[["percent.mt"]] <- PercentageFeatureSet(
  pbmc,
  pattern = "^MT-"
)


# QC violin plots

qc_violin <- VlnPlot(
  pbmc,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)

ggsave(
  "results/qc/qc_violin_plots.png",
  plot = qc_violin,
  width = 12,
  height = 4
)


# UMI count vs genes detected

qc_scatter <- FeatureScatter(
  pbmc,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

ggsave(
  "results/qc/umi_vs_genes.png",
  plot = qc_scatter,
  width = 6,
  height = 5
)


# QC thresholds

min_features <- 200
max_features <- 2500
max_mito <- 5


# Filter cells

pbmc <- subset(
  pbmc,
  subset =
    nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    percent.mt < max_mito
)

cat("Cells after QC:", ncol(pbmc), "\n")


# ------------------------------------------------------------
# 6. Normalize expression
# ------------------------------------------------------------

cat("\n[4/12] Normalizing gene expression...\n")

pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


# ------------------------------------------------------------
# 7. Identify highly variable genes
# ------------------------------------------------------------

cat("\n[5/12] Identifying highly variable genes...\n")

pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000
)

cat(
  "Variable genes:",
  length(VariableFeatures(pbmc)),
  "\n"
)

variable_plot <- VariableFeaturePlot(pbmc)

ggsave(
  "results/variable_features/highly_variable_genes.png",
  plot = variable_plot,
  width = 10,
  height = 6
)


# ------------------------------------------------------------
# 8. Scale the variable genes
# ------------------------------------------------------------

cat("\n[6/12] Scaling expression data...\n")

pbmc <- ScaleData(
  pbmc,
  features = VariableFeatures(pbmc)
)


# ------------------------------------------------------------
# 9. Principal component analysis
# ------------------------------------------------------------

cat("\n[7/12] Running PCA...\n")

pbmc <- RunPCA(
  pbmc,
  features = VariableFeatures(pbmc)
)

elbow_plot <- ElbowPlot(
  pbmc,
  ndims = 30
)

ggsave(
  "results/pca/elbow_plot.png",
  plot = elbow_plot,
  width = 8,
  height = 6
)


# ------------------------------------------------------------
# 10. Neighbor graph and clustering
# ------------------------------------------------------------

cat("\n[8/12] Constructing neighbor graph...\n")

pbmc <- FindNeighbors(
  pbmc,
  dims = 1:10
)

cat("\n[9/12] Clustering cells...\n")

# Fixed seed improves reproducibility of the workflow
set.seed(1234)

pbmc <- FindClusters(
  pbmc,
  resolution = 0.5
)

cat("Number of clusters:", length(levels(Idents(pbmc))), "\n")
print(table(Idents(pbmc)))


# ------------------------------------------------------------
# 11. UMAP
# ------------------------------------------------------------

cat("\n[10/12] Running UMAP...\n")

set.seed(1234)

pbmc <- RunUMAP(
  pbmc,
  dims = 1:10
)

umap_plot <- DimPlot(
  pbmc,
  reduction = "umap",
  label = TRUE
)

ggsave(
  "results/umap/umap_clusters.png",
  plot = umap_plot,
  width = 8,
  height = 6
)


# ------------------------------------------------------------
# 12. Identify cluster markers
# ------------------------------------------------------------

cat("\n[11/12] Identifying cluster markers...\n")

markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  markers,
  "results/markers/all_cluster_markers.csv",
  row.names = FALSE
)

cat("Marker rows identified:", nrow(markers), "\n")


# ------------------------------------------------------------
# 13. Canonical marker DotPlot
# ------------------------------------------------------------

marker_genes <- c(
  "CD3D",
  "CD3E",
  "CD79A",
  "MS4A1",
  "CD37",
  "NKG7",
  "GNLY",
  "GZMB",
  "LYZ",
  "S100A8",
  "S100A9",
  "FCN1",
  "FCER1A",
  "CST3",
  "CD1C",
  "PPBP",
  "PF4"
)

marker_dotplot <- DotPlot(
  pbmc,
  features = marker_genes
) +
  RotatedAxis()

ggsave(
  "results/annotation/cluster_marker_dotplot.png",
  plot = marker_dotplot,
  width = 12,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 14. Cell-type annotation
# ------------------------------------------------------------

cat("\n[12/12] Assigning cell-type annotations...\n")

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

pbmc <- RenameIdents(
  pbmc,
  cluster_labels
)

pbmc$cell_type <- Idents(pbmc)

cat("\nFinal cell-type counts:\n")
print(table(pbmc$cell_type))


# ------------------------------------------------------------
# 15. Final annotated UMAP
# ------------------------------------------------------------

annotated_umap <- DimPlot(
  pbmc,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  label.size = 4
) +
  NoLegend() +
  ggtitle("PBMC Cell Type Annotation")

ggsave(
  "results/annotation/annotated_umap.png",
  plot = annotated_umap,
  width = 10,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 16. Save final Seurat object
# ------------------------------------------------------------

saveRDS(
  pbmc,
  "results/annotation/pbmc_annotated.rds"
)


# ------------------------------------------------------------
# Complete
# ------------------------------------------------------------

cat("\n============================================\n")
cat("PBMC analysis complete\n")
cat("============================================\n")

cat(
  "\nFinal cells:",
  ncol(pbmc),
  "\nFinal genes:",
  nrow(pbmc),
  "\nCell types:",
  length(levels(Idents(pbmc))),
  "\n"
)

cat(
  "\nFinal Seurat object saved to:\n",
  "results/annotation/pbmc_annotated.rds\n"
)
