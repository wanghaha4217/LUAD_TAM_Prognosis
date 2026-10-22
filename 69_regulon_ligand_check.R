# ============================================================
# 69_regulon_ligand_check.R
# 方向一核心：查通讯配体/受体是不是 NF-κB / HIF / STAT 的靶基因
# 目的：把 master regulator（NF-κB）跟通讯轴（SPP1/TIM-3 排斥 + CXCL 趋化）串起来
# 运行：本地 RStudio（dorothea 已装），秒级
# ============================================================

suppressPackageStartupMessages(library(dorothea))

data(dorothea_hs, package = "dorothea")
reg <- dorothea_hs[dorothea_hs$confidence %in% c("A", "B"), ]

# 三轴的关键基因：免疫排斥(SPP1) / 检查点(LGALS9-HAVCR2) / 趋化(CXCL-CXCR3)
genes <- c("SPP1", "LGALS9", "HAVCR2", "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CXCR3")

# 候选 master TF（P8 的候选 + IFN 轴对照）
tfs <- c("RELA", "NFKB1", "HIF1A", "EPAS1", "STAT1", "STAT2", "IRF1")

cat("===== 通讯基因的 DoRothEA 调控关系（A+B 置信，mor=+1 激活 / -1 抑制）=====\n\n")
for (g in genes) {
  hits <- reg[reg$target == g & reg$tf %in% tfs, ]
  if (nrow(hits) > 0) {
    for (i in seq_len(nrow(hits))) {
      cat(sprintf("%-8s  <-  %-6s  mor=%+d  conf=%s\n",
                  g, hits$tf[i], hits$mor[i], hits$confidence[i]))
    }
  } else {
    cat(sprintf("%-8s  <-  无 %s 调控记录\n", g, paste(tfs, collapse = "/")))
  }
}

cat("\n===== 关键问题：SPP1 / LGALS9 是不是 NF-κB 的激活靶基因？ =====\n")
for (g in c("SPP1", "LGALS9")) {
  nfkb <- reg[reg$target == g & reg$tf %in% c("RELA", "NFKB1") & reg$mor == 1, ]
  if (nrow(nfkb) > 0) {
    cat(sprintf("✅ %s 是 NF-κB(RELA/NFKB1) 激活靶基因（%s）\n",
                g, paste(unique(nfkb$tf), collapse = ", ")))
  } else {
    cat(sprintf("❌ %s 不是 NF-κB 激活靶基因（A+B 置信里无记录）\n", g))
  }
}
