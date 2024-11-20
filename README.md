# single-cell-RNA-sequencing
# 1 download data
## 1.1 experiment data
GSE144240
```
mkdir -p ~/project/scRNA
cd ~/project/scRNA
mkdir genome sequence annotation
cd sequence
nohup prefetch SRR11832836 SRR11832837 -O . &

parallel -j 2 "
  fastq-dump --split-files --gzip {1}
" ::: $(ls *.sra)

rm *.sra
```
## 1.2 ref data
ensembl-human
```
cd ../genome
aria2c -d ./ -Z https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gzip -d Homo_sapiens.GRCh38.cdna.all.fa.gz
mv Homo_sapiens.GRCh38.cdna.all.fa refgenome.fa
```
## 1.3 annotation data
```
cd ../annotation
aria2c -d ./ -Z https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
gzip -d Homo_sapiens.GRCh38.113.gtf.gz
mv Homo_sapiens.GRCh38.113.gtf hg.gtf
```
# 2 quality control and trimming
* fastqc
```
cd ../sequence
fastqc --threads 3 *.fastq.gz
```
* trimmomaric
```
parallel -j 4 "
  java -jar ~/biosoft/Trimmomatic-0.39/trimmomatic-0.39.jar \
  PE -phred33 {1} ./trim/{1} \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:5:15 MINLEN:30 \
" ::: $(ls *.gz)
```
# 3 seq alignment
* hisat2
make index
```
cd ~/project/scRNA/genome
mkdir index
cd index

hisat2-build -p 6 ../*.fasta hg_index
```
seq alignment
```
mkdir ../sequence/align
cd ../sequence/trim
parallel -k -j 4 "
    hisat2 -p 8 -x ~/project/scRNA/genome/index/hg_index \
    -1 sra_data/SRR11832836_1_trimmed.fastq \
    -2 sra_data/SRR11832836_2_trimmed.fastq \
    -S ../align/{1}.sam 2>../align/{1}.log
" ::: $(ls *.gz | perl -p -e 's/.fastq.gz$//')

cd ../align
parallel -j 4 "
    samtools sort -@ 2 {1}.sam > {1}.sort.bam && samtools index {1}.sort.bam
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')

rm *.sam
ls
```
# 4 gene expression count
* HTseq
```
cd ~/project/scRNA/sequence
mkdir HTseq

cd align
parallel =j 4 "
    htseq-count -f bam -s yes {1}.sort.bam ../../annotation/hg.gtf \
        > ../HTseq/{1}.count 2>../HTseq/{1}.log
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')

cd ../HTseq
cat SRR11832836.count | head -n 10
```
# 5 merge & normalize
## 5.1 merge data
```
R
rm(list=ls())

files <- list.files(".", "*.count")
f_lists <- list()
for (i in files){
    perfix = gsub("(_\\w+)?\\.count", "", i, perl=TRUE)
    f_lists[[perfix]]= i
}

id_list <- names(f_lists)
data <- list()
count <- 0
for (i in id_list){
    count <- count + 1
    a <- read.table(f_lists[[i]], sep="\t", col.names = c("gene_id", i))
    data[[count]] <- a
}

data_merge <- data[[1]]
for (i in seq(2, length(id_list))){
    data_merge <- merge(data_merge, data[[i]], by="gene_id")
}

write.csv(data_merge, "merge.csv", quote = FALSE, row.names = FALSE)
```
## 5.2 normalize
```
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = data.frame(condition = rep("control", ncol(count_matrix))), design = ~ 1)

# 估计大小因子并进行标准化
dds <- estimateSizeFactors(dds)
normalized_counts_deseq <- counts(dds, normalized = TRUE)

# 查看标准化后的计数结果
head(normalized_counts_deseq)
```
# 6 sc-RNA analysis
## 6.1 create Seurat object
```
# 安装并加载Seurat包
if (!requireNamespace("Seurat", quietly = TRUE))
    install.packages("Seurat")
library(Seurat)

# 创建Seurat对象
seurat_obj <- CreateSeuratObject(count_matrix)

# 过滤细胞
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 10)  # 根据文章中过滤条件，可调整参数

# 数据标准化
seurat_obj <- NormalizeData(seurat_obj)

# 寻找高变基因
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)  # 可调整高变基因数量
```
## 6.2 cluster and visualize
```
# PCA降维
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# 确定PCA主成分数量（可根据肘部法则等方法选择）
ElbowPlot(seurat_obj)

# UMAP可视化
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)  # 根据确定的主成分数量调整

# 细胞聚类
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)  # 可调整聚类分辨率

# 可视化聚类结果
DimPlot(seurat_obj, reduction = "umap")
```