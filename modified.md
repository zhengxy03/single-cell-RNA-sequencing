# single-cell-RNA-sequencing
# 1 download data
GSE144240-GSE144236
```
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144236/suppl/GSE144236%5FSCC13%5Fcounts.txt.gz
gzip -d GSE144236_SCC13_counts.txt.gz
mv GSE144236_SCC13_counts.txt scc13.gz
``` 
# 2 Seurat quality control and normalize
## 2.1 quality control
```
rm(list = ls())
if (!requireNamespace("Seurat", quietly = TRUE))
    install.packages("Seurat")
library(Seurat)

setwd("//wsl.localhost/Ubuntu/home/zxy0303/project/scRNA")

counts_matrix <- read.table("scc13.txt", header = TRUE, sep = "\t")
seurat_obj <- CreateSeuratObject(counts = counts_matrix)

seurat_obj[["percent.mito"]] <- PercentageFeatureSet(seurat_obj, pattern = "^hg19-MT-")

VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 10)

#merge different samples
#merged_seurat_obj <- merge(seurat_obj1, y = c(seurat_obj2), add.cell.ids = c("Sample1", "Sample2"))
```
## 2.2 normalized
标准化的意义：
> 数据标准化的意义: 去除测序深度带来的影响<br>
> 标准化原则：每个细胞的每个基因的count数除以该细胞总count数，然后乘以因子（10000），再进行log(n+1)转换<br>

normalization函数，默认LogNormalize的方法
```
seurat_obj <- NormalizeData(seurat_obj)
```

## find the most variable genes
FindVariableFeatures（）参数意义：
* FindVariableFeatures 函数有 3 种选择高表达变异基因的方法，可以通过 selection.method参数来选择，它们分别是： vst（默认值）， mean.var.plot 和 dispersion。 nfeatures 参数的默认值是 2000，可以改变。如果 selection.method 参数选择的是 mean.var.plot，就不需要人为规定高表达变异基因的数目，算法会自动选择合适的数目。 建议使用完 FindVariableFeatures 函数后，用 VariableFeaturePlot 对这些高表达变异基因再做个可视化，看看使用默认值 2000 会不会有问题。

```
#0.0125 <非零值均值 < 3 且标准差> 0.5
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))


variable.genes <- VariableFeatures(object = seurat_obj)
if (is.null(variable.genes)) {
  print("没有识别到高变基因")
} else {
  print(paste("识别到", length(variable.genes), "个高变基因"))
}
```
# 3 PCA analysis
ScaleData()标准化函数作用：<br>
为后续PCA降维做准备。<br>
PCA降维要求数据为正态分布，即平均值为0，方差为1。<br>
```
#回归 UMI 计数和线粒体基因百分比
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nCount_RNA", "percent.mito"))
```
```
#PCA降维
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

#选择前5个维度进行查看
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(seurat_obj, reduction = "pca")

#查看pbmc数据集中PCA降维后的细胞嵌入的前两行和前两列的数据。
head(seurat_obj[['pca']]@cell.embeddings)[1:2,1:2]
#                     PC_1      PC_2
#AAACCTGAGCCACGCT 7.965446 -3.265728
#AAACCTGGTAGCCTCG 7.356585 -6.612420


VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
```
DimPlot（）函数生成主成分分析结果图
![Dimplot](./pca1, "Dimplot")
VizDimLoadings
![vizdimloading](./pca2, "vizdimloading")

# 确定合适的主成分数（使用 ElbowPlot 可视化
确定数据集的维度<br>
> 目的：每个维度（pc）本质上代表一个“元特征”，它将相关特征集中的信息组合在一起。因此，越在顶部的主成分越可能代表数据集。然而，我们应该选择多少个主成分才认为我们选择的数据包含了绝大部分的原始数据信息呢？

方法<br>

> ElbowPlot函数，基于每个主成分所解释的方差百分比的排序，通过寻找“拐点”来判断几个维度可包含数据的大部分信息。

```
ElbowPlot(seurat_obj)
```
![elbowplot](./elbowplot, "elbowplot")

# UMAP降维聚类
基于选定的主成分（假设根据 ElbowPlot 确定为前 15 个主成分）进行 UMAP 降维并绘制 UMAP 图：
```
# 执行UMAP降维
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

# 绘制UMAP图（可根据需要调整参数，如label是否显示标签等）
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
```
![umap](./umap, "umap")
# unsupervised clustering
```
# 进行无监督聚类

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.1, dims = 1:15)


head(Idents(seurat_obj), 5)

```
# 差异表达分析
使用 FindAllMarkers 函数进行差异表达分析，获取每个聚类的标记基因（默认参数）：
```
# 进行差异表达分析
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)


```
FindAllMarkers（）参数意义：

only.pos = TRUE：只寻找上调的基因
min.pct = 0.1：某基因在细胞中表达的细胞数占相应cluster细胞数最低10%
logfc.threshold = 0.25 ：fold change倍数为0.25
该函数是决速步，执行比较耗时。
```
head(cluster_markers, n = 5)
                      p_val avg_log2FC pct.1 pct.2     p_val_adj cluster
hg19-EEF1A1   9.262398e-167  0.6032398 1.000 1.000 1.731142e-162       0
hg19-MIR205HG 4.184558e-158  0.8288175 1.000 0.999 7.820938e-154       0
hg19-RPL10    1.842909e-139  0.4595613 1.000 1.000 3.444397e-135       0
hg19-RPL3     1.866003e-132  0.5649401 1.000 1.000 3.487559e-128       0
hg19-DST      4.314248e-131  1.1875805 0.949 0.795 8.063330e-127       0
                       gene
hg19-EEF1A1     hg19-EEF1A1
hg19-MIR205HG hg19-MIR205HG
hg19-RPL10       hg19-RPL10
hg19-RPL3         hg19-RPL3
hg19-DST           hg19-DST
```

```
cluster_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
# A tibble: 6 × 7
# Groups:   cluster [3]
      p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene              
      <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>             
1 2.34e-  4       2.40 0.145 0.108 1   e+  0 0       hg19-HCAR3        
2 3.39e-  3       2.83 0.01  0.001 1   e+  0 0       hg19-RP11-346D14.1
3 0               5.66 0.657 0.025 0         1       hg19-NDC80        
4 8.38e-201       5.61 0.382 0.013 1.57e-196 1       hg19-HIST1H1B     
5 9.68e- 39       7.48 0.14  0.014 1.81e- 34 2       hg19-FLG          
6 1.76e- 38       8.69 0.077 0.002 3.29e- 34 2       hg19-SPRR2B     
```
```

#VlnPlot: 基于细胞类群的基因表达概率分布
VlnPlot(seurat_obj, features = c("hg19-EEF1A1", "hg19-RPL10"))
```
![vlnplot](./vln, "vln")
# 细胞类型注释
根据标记基因的表达情况确定细胞类型，这里通过检查每个聚类中特定细胞类型标记基因的表达来进行注释（示例代码，实际可能需要根据更详细的文献知识和数据特点进行调整）：
```
annotate_cell_type <- function(markers) {
  cell_type <- NA
  if (any(markers %in% c("hg19-KRT5", "hg19-KRT14", "hg19-KRT1", "hg19-KRT10"))) {
    cell_type <- "Epithelial"
  } else if (any(markers %in% c("hg19-LYZ", "hg19-HLA-DRB1", "hg19-HLA-DRA", "hg19-HLA-DQB2"))) {
    cell_type <- "Myeloid"
  } else if (any(markers %in% c("hg19-CD3D", "hg19-CD2", "hg19-CD7"))) {
    cell_type <- "T_cell"
  } else if (any(markers %in% c("hg19-COL1A1", "hg19-COL1A2", "hg19-LUM"))) {
    cell_type <- "Fibroblast"
  } else if (any(markers %in% c("hg19-MLANA", "hg19-DCT", "hg19-PMEL"))) {
    cell_type <- "Melanocyte"
  } else if (any(markers %in% c("hg19-TFF3", "hg19-CLDN5", "hg19-VWF"))) {
    cell_type <- "Endothelial"
  } else if (any(markers %in% c("hg19-IGLL5", "hg19-IGJ", "hg19-MS4A1", "hg19-CD79A"))) {
    cell_type <- "B/plasma"
  }
  return(cell_type)
}

cluster_cell_types <- sapply(cluster_markers$gene, annotate_cell_type)
seurat_obj$cell_type <- cluster_cell_types[seurat_obj@active.ident]