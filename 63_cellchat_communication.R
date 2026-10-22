# ============================================================
# 63_cellchat_communication.R
# CellChat 通讯分析：non-MPR（ICB 耐药）vs MPR/pCR（应答）
# 目标：锚定 P8 缺氧轴 / SPP1 程序到「免疫排斥/检查点」配受体轴
# CellChat v2 API（云端 2.2.0.9001）
# 运行：云端 RStudio（内存 ~50GB，两组各十几万细胞，computeCommunProb 每组约 5-15 分钟）
# 输出：./cellchat/
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
})

DIS     <- "."
out_dir <- file.path(DIS, "cellchat")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. 加载对象 + 归一化 + 分组 ----------
so <- readRDS(file.path(DIS, "GSE243013_LUAD_immune_cellchat.rds"))
so <- NormalizeData(so)
so$response <- ifelse(so$pathological_response == "non-MPR", "nonMPR", "MPR_pCR")
cat("分组:\n"); print(table(so$response))

# ---------- 2. CellChat 函数（v2 流程） ----------
run_cellchat <- function(so_sub, tag) {
  cat("\n===== 开始 CellChat:", tag, "=====\n")
  if (!requireNamespace("presto", quietly = TRUE)) {
    stop("缺少 presto 包：请先运行 remotes::install_github('immunogenomics/presto')，然后重跑本脚本")
  }
  data.input <- GetAssayData(so_sub, assay = "RNA", layer = "data")
  meta <- so_sub@meta.data
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)                     # 降内存
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = TRUE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cat("===== 完成 CellChat:", tag, "=====\n")
  return(cellchat)
}

# ---------- 3. 分别跑两组 ----------
so_nompr <- subset(so, response == "nonMPR")
so_mpr   <- subset(so, response == "MPR_pCR")
rm(so); gc()

cc_nompr <- run_cellchat(so_nompr, "nonMPR")
cc_mpr   <- run_cellchat(so_mpr, "MPR_pCR")

# ---------- 4. 合并对比 ----------
object.list <- list(nonMPR = cc_nompr, MPR_pCR = cc_mpr)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

saveRDS(list(nonMPR = cc_nompr, MPR_pCR = cc_mpr, merged = cellchat),
        file.path(out_dir, "cellchat_results.rds"))
cat("\n已保存 cellchat_results.rds\n")

# ---------- 5. 差异图 ----------
# (A) 差异互作数量（红=nonMPR 增强，蓝=MPR_pCR 增强）
ggA <- netVisual_diffInteraction(cellchat, weight.scale = TRUE)
ggsave(file.path(out_dir, "diff_interaction_number.png"), ggA, width = 10, height = 8, dpi = 300)

# (B) 差异互作强度
ggB <- netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
ggsave(file.path(out_dir, "diff_interaction_strength.png"), ggB, width = 10, height = 8, dpi = 300)

# (C) 差异通路排序（含统计）
pdf(file.path(out_dir, "rankNet_comparison.pdf"), width = 8, height = 12)
rankNet(cellchat, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE)
dev.off()

# ---------- 6. 提取目标配受体（SPP1/CXCL/LGALS9/C1Q 轴） ----------
# 合并两组的显著配受体对，筛目标基因
df_nompr <- subsetCommunication(cc_nompr)
df_mpr   <- subsetCommunication(cc_mpr)
df_nompr$group <- "nonMPR"; df_mpr$group <- "MPR_pCR"
df_all <- rbind(df_nompr, df_mpr)

target <- grepl("SPP1|CD44|ITGA4|ITGB1|C1Q|GPR17|CXCL9|CXCL10|CXCL11|CXCR3|LGALS9|HAVCR2",
                paste(df_all$ligand, df_all$receptor, df_all$pathway_name, sep = "|"))
df_target <- df_all[target, ]
if (nrow(df_target) > 0) {
  write.csv(df_target, file.path(out_dir, "target_ligand_receptor.csv"), row.names = FALSE)
  cat("\n目标配受体对数量:", nrow(df_target), "\n")
  print(unique(df_target[, c("ligand", "receptor", "pathway_name", "group")]))
} else {
  cat("\n未找到目标配受体对（SPP1/CD44/CXCL/LGALS9/C1Q），说明这些轴在本数据通讯弱\n")
}

# ---------- 7. 目标轴气泡图（巨噬 -> T 细胞的 SPP1/CXCL 等） ----------
# 巨噬类（Mφ_*）为 source，T 类为 target，看 SPP1/CXCL/LGALS9 等通路
mac_ids <- c("Mφ_CXCL10","Mφ_CXCL2","Mφ_DNAJB1","Mφ_FCGR3A","Mφ_FOLR2","Mφ_ISG15",
             "Mφ_MARCO","Mφ_MKI67","Mφ_MMP9","Mφ_S100A10","Mφ_VCAN")
t_ids   <- c("CD4T","CD8T","CD8T_Tex","Treg","NK","gdT")
sig_target <- c("SPP1", "C1Q", "CXCL", "LGALS9", "MHC-II")
for (sig in sig_target) {
  p <- tryCatch(
    netVisual_bubble(cellchat, sources.use = mac_ids, targets.use = t_ids,
                     signaling = sig, remove.isolate = FALSE),
    error = function(e) NULL)
  if (!is.null(p)) {
    ggsave(file.path(out_dir, paste0("bubble_", sig, ".png")), p,
           width = 10, height = 8, dpi = 300)
  } else {
    cat("\n", sig, " 通路在当前数据无显著通讯（跳过气泡图）\n")
  }
}

cat("\n===== 全部完成，输出在", out_dir, "=====\n")
