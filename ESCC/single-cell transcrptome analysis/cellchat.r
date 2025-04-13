data.input <- GetAssayData(epi_fib, assay = "RNA", slot = "data")  # 使用 log-normalized 数据

meta <- epi_fib@meta.data
counts <- GetAssayData(epi_fib, assay = "RNA", slot = "counts")  # 原始计数
data <- GetAssayData(epi_fib, assay = "RNA", slot = "data")      # 标准化数据

# 创建对象（显式指定两个层）
cellchat <- createCellChat(
  object = list(raw = counts, normalized = data),  # 关键！
  meta = epi_fib@meta.data,
  group.by = "cell_type"
)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")
cellchat <- createCellChat(
  object = as.matrix(data.input),  # 确保输入是矩阵
  meta = meta,
  group.by = "cell_type"
)

# 设置细胞分组信息
cellchat <- setIdent(cellchat, ident.use = "cell_type")
levels(cellchat@idents)

CellChatDB <- CellChatDB.human  # 人类数据库
# CellChatDB <- CellChatDB.mouse  # 小鼠数据库

# 使用全部相互作用或子集
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")  # 可选：仅用分泌信号
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)  # 可选：子集数据加速计算
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = FALSE) 