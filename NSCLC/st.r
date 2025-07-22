#10× Visium
setwd("E:/project/nsclc/GSE189487_RAW")
list.files()

library(Seurat)
library(Matrix)
library(dplyr)

data_dir <- "TD1/"
expr <- ReadMtx(
    mtx = file.path(data_dir, "matrix.mtx.gz"),
    features = file.path(data_dir, "features.tsv.gz"),
    cells = file.path(data_dir, "barcodes.tsv.gz")
)
coords <- read.csv(file.path(data_dir, "tissue_positions_list.csv.gz"), header=FALSE)
colnames(coords) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
rownames(coords) <- coords$barcode
coords <- coords[, c("pxl_row_in_fullres", "pxl_col_in_fullres")]

seu <- CreateSeuratObject(counts = expr, assay = "Spatial")
common_cells <- intersect(colnames(seu), rownames(coords))
seu <- subset(seu, cells = common_cells)
coords <- coords[common_cells, , drop = FALSE]
seu[["image"]] <- new(Class = "SlideSeq", coordinates = coords, assay = "Spatial")

seu <- SCTransform(seu, assay = "Spatial", verbose = FALSE)

seu <- RunPCA(seu, verbose = FALSE)

seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)

seu <- RunUMAP(seu, dims = 1:30)
SpatialPlot(seu, group.by = "seurat_clusters", image.alpha = 0)

markers <- FindAllMarkers(seu,
                          only.pos = TRUE,
                          test.use = "wilcox",
                          min.pct = 0.1,
                          logfc.threshold = 0.25)

marker_genes <- markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25)

#GeoMx DSP
BiocManager::install(c("NanoStringNCTools", "limma", "edgeR"))
remotes::install_github("Nanostring-Biostats/GeoMxTools", dependencies = TRUE)

library(GeomxTools)
library(NanoStringNCTools)
library(limma)
library(edgeR)

dcc_gz <- list.files(pattern = "\\.dcc\\.gz$", full.names = TRUE)
sapply(dcc_gz, R.utils::gunzip, overwrite = TRUE)
dcc_files <- list.files(pattern = "\\.dcc$", full.names = TRUE)
raw_data <- readNanoStringGeoMxSet(dccFiles = dcc_files)
