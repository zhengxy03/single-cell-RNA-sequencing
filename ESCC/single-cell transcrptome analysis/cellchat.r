library(Seurat)
library(CellChat)
library(patchwork)
library(dplyr)

data.input <- epi_fib[["RNA"]]$data

meta <- epi_fib@meta.data
cell.types <- meta$cell_type
names(cell.types) <- rownames(meta)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")

CellChatDB <- CellChatDB.human  # 人类用
# CellChatDB <- CellChatDB.mouse # 小鼠用

# 使用全部数据库或子集（可选）
cellchat@DB <- CellChatDB


cellchat <- subsetData(cellchat) # 必要时裁剪数据
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)

# 过滤低可信度的相互作用（可选）
cellchat <- filterCommunication(cellchat, min.cells = 10)

# 聚合细胞类型的通讯网络
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(meta$cell_type))
# 保存为PNG
png("cellchat_network.png", width=2000, height=1800, res=300)  # 高分辨率
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T)
dev.off()