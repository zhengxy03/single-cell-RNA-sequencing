Dedifferentiation_genes <- c("CDHR1", "COL17A1", "ITGA6", "NOTCH1", "ABCG2", "AFMID", "ALOX5", "ANPEP", "ANXA5", "APC", "ASCL2", "AXIN2", "AZGP1", "BOC", "CD38", "CDCA7", "CDK6", "CFTR", "CHD7", "CREM", "CRHBP", "CXCL5", "DNER", "DNMT3A", "DPP4", "EBF1", "EIF4B", "EMX2", "EPHB2", "ETS2", "ETV1", "FERMT1", "FGF9", "FOXA2", "FOXA3", "FOXG1", "GBP2", "GFAP", "GFI1", "GLIPR1", "GUCY1A1", "H3F3B", "HAPLN1", "HBB", "HES1", "HIST1H2AC", "HIST1H2BK", "HLX", "HOPX", "IDH1", "IKZF1", "IL18R1", "KLK10", "KRT14", "LATS2", "LRIG1", "MESP2", "METTL3", "MLLT10", "MLLT3", "MPZL2", "MYB", "NANOG", "NEK5", "NELL2", "NFAT5", "NFE2", "NFIA", "NKX2-5", "NODAL", "NR4A2", "NRIP1", "OPHN1", "OPTN", "PABPC1", "PAX6", "POU5F1", "PTMA", "PTPRG", "PTPRO", "RBM6", "RBPMS", "RGMB", "RGS1", "RHAG", "RNF43", "SET", "SLC12A2", "SMAD2", "SMOC2", "SOX1", "SOX17", "SOX3", "SOX9", "SPDYE1", "SPTBN1", "SVIL", "TCF12", "TCF4", "TCF7L2", "TFPI", "TFRC", "TGFB1I1", "TM4SF1", "TNFSF10", "TOX3", "TPBG", "TRA2B", "TSPAN31", "UGT8", "VNN1", "ZBTB10", "ZNF793")
fibroblasts_genes <- rownames(fibroblasts)
matched_genes <- Dedifferentiation_genes[Dedifferentiation_genes %in% fibroblasts_genes]
length(matched_genes)
Dedifferentiation_genes <- matched_genes

fibroblasts <- AddModuleScore(fibroblasts, features = list(Dedifferentiation_genes), name = "Dedifferentiation_Score")
head(fibroblasts@meta.data)
summary(fibroblasts@meta.data$Dedifferentiation_Score1)
fibroblasts@meta.data$Dedifferentiation_Score_Z <- scale(fibroblasts@meta.data$Dedifferentiation_Score1)
summary(fibroblasts@meta.data$Dedifferentiation_Score_Z)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(4)



p <- FeaturePlot(
  object = fibroblasts,
  features = "Dedifferentiation_Score1",
  label = TRUE,
  repel = TRUE,
  order = TRUE,
  pt.size = 1
) +
  scale_color_gradient2(
    low = "#4DBBD5FF",  # 低值颜色
    mid = "white",      # 中值颜色
    high = "#E64B35FF", # 高值颜色
    midpoint = 0,       # 中值点（对应Stemness_Score1的中值）
    limits = c(-0.6, 0.8)  # 固定图例范围
  )
ggsave("fibro_stem_umap.png", plot = p)

npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(15)
p <- VlnPlot(
  object = fibroblasts,
  features = "Dedifferentiation_Score_Z",
  group.by = "cell_type",
  pt.size = 0,
  cols = npg_extended
) +
  theme(legend.position = "none",
  axis.title.x = element_blank(), 
  axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("fibro_stem_vln.png", plot = p)