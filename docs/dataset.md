# Dataset

## PBMC 3K

This project uses the 10x Genomics PBMC 3K single-cell RNA-seq dataset.

The dataset contains peripheral blood mononuclear cells (PBMCs) from a healthy donor and was generated using the 10x Genomics Chromium single-cell platform.

For this project, analysis begins with the Cell Ranger filtered gene-barcode count matrix rather than raw FASTQ files.

## Input Files

The filtered matrix contains:

* 2,700 cell barcodes
* 32,738 genes
* 2,286,884 non-zero gene-cell entries

The Cell Ranger output consists of three files:

* `barcodes.tsv` — identifies the cell barcodes represented by matrix columns
* `genes.tsv` — identifies the genes represented by matrix rows
* `matrix.mtx` — sparse matrix containing UMI counts

The raw dataset is excluded from Git because it is publicly available and can be downloaded separately.

## Analysis Goal

The goal of this project is to identify and characterize the major immune cell populations present in the PBMC sample using single-cell gene expression profiles.
