# Load Seurat
library(Seurat)
library(ggplot2)

# Location of the 10x Genomics count matrix
data_dir <- "data/raw/filtered_gene_bc_matrices/hg19"

# Read the count matrix
pbmc_counts <- Read10X(data.dir = data_dir)

# Create the Seurat object
pbmc <- CreateSeuratObject(
    counts = pbmc_counts,
    project = "PBMC3K"
)

# Calculate the percentage of counts coming from mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(
    pbmc,
    pattern = "^MT-"
)

# Display the first five cells' QC metrics
print(head(pbmc[[]], 5))

# Summarize the three QC metrics in one table
qc_summary <- data.frame(
    Metric = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    Min = c(
        min(pbmc$nFeature_RNA),
        min(pbmc$nCount_RNA),
        min(pbmc$percent.mt)
    ),
    Q1 = c(
        quantile(pbmc$nFeature_RNA, 0.25),
        quantile(pbmc$nCount_RNA, 0.25),
        quantile(pbmc$percent.mt, 0.25)
    ),
    Median = c(
        median(pbmc$nFeature_RNA),
        median(pbmc$nCount_RNA),
        median(pbmc$percent.mt)
    ),
    Mean = c(
        mean(pbmc$nFeature_RNA),
        mean(pbmc$nCount_RNA),
        mean(pbmc$percent.mt)
    ),
    Q3 = c(
        quantile(pbmc$nFeature_RNA, 0.75),
        quantile(pbmc$nCount_RNA, 0.75),
        quantile(pbmc$percent.mt, 0.75)
    ),
    Max = c(
        max(pbmc$nFeature_RNA),
        max(pbmc$nCount_RNA),
        max(pbmc$percent.mt)
    )
)

print(qc_summary, row.names = FALSE)


# Create violin plots for the three QC metrics
qc_violin <- VlnPlot(
    pbmc,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
)

# Save the violin plots
ggsave(
    "results/qc/qc_violin_plots.png",
    plot = qc_violin,
    width = 12,
    height = 4
)

# Create a scatter plot comparing total UMIs with genes detected
qc_scatter <- FeatureScatter(
    pbmc,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
)

# Save the scatter plot
ggsave(
    "results/qc/umi_vs_genes.png",
    plot = qc_scatter,
    width = 6,
    height = 5
)

# Candidate QC thresholds
min_features <- 200
max_features <- 2500
max_mito <- 5

# Count cells excluded by each individual criterion
qc_exclusions <- data.frame(
    Criterion = c(
        "Too few genes",
        "Too many genes",
        "High mitochondrial percentage"
    ),
    Cells_excluded = c(
        sum(pbmc$nFeature_RNA <= min_features),
        sum(pbmc$nFeature_RNA >= max_features),
        sum(pbmc$percent.mt >= max_mito)
    )
)

print(qc_exclusions, row.names = FALSE)

# Determine which cells pass all QC criteria
pass_qc <- (
    pbmc$nFeature_RNA > min_features &
    pbmc$nFeature_RNA < max_features &
    pbmc$percent.mt < max_mito
)

print(table(pass_qc))
# Filter cells using the QC criteria
pbmc_filtered <- subset(
    pbmc,
    subset = nFeature_RNA > min_features &
             nFeature_RNA < max_features &
             percent.mt < max_mito
)

# Report the number of cells before and after filtering
cat("Cells before QC:", ncol(pbmc), "\n")
cat("Cells after QC:", ncol(pbmc_filtered), "\n")

# Create the QC results directory if it does not exist
dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

# Save the filtered Seurat object
saveRDS(
    pbmc_filtered,
    file = "results/qc/pbmc_filtered.rds"
)
