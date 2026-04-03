library(Seurat)
library(monocle)
library(ggplot2)
library(ggsci)
library(igraph)  # 添加这一行

# 加载完整修复脚本
source("monocle_fix_complete.R")

# 现在可以正常使用
cds <- readRDS("malignant_unknown_cds_0.6_20000cells_hvg.rds")
cds <- orderCells(cds)
