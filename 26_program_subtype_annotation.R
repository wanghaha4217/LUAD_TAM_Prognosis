# ============================================================
# 26_program_subtype_annotation.R
# 打印「程序 ↔ 亚型」注释表 + 每程序 top 基因，用于给 P1~P10 起生物学名字
# 依赖：Seurat（仅读 meta）；已验证 H 列顺序 = Cells 顺序（cor=1）
# 运行：本地 RStudio
# ============================================================
suppressPackageStartupMessages({ library(Seurat) })

base    <- "."
res_dir <- file.path(base, "results")

obj <- readRDS(file.path(base, "GSE243013_LUAD_TAM_clustered.rds"))
nmf <- readRDS(file.path(base, "nmf_k10.rds"))

W <- nmf$W; H <- nmf$H; genes <- nmf$genes; programs <- nmf$programs
subtype <- as.character(obj@meta.data$sub_cell_type)

usub <- sort(unique(subtype))
np   <- nrow(H)

# ---- 程序 ↔ 亚型 使用度 ----
mean_usage <- matrix(NA_real_, np, length(usub),
                     dimnames = list(paste0("P", seq_len(np)), usub))
frac_top   <- matrix(NA_real_, np, length(usub),
                     dimnames = list(paste0("P", seq_len(np)), usub))
for (s in usub) {
  idx <- which(subtype == s)
  Hs  <- H[, idx, drop = FALSE]
  mean_usage[, s] <- rowMeans(Hs)
  frac_top[, s]   <- tabulate(apply(Hs, 2, which.max), nbins = np) / length(idx)
}

cat("===== 各亚型细胞数 =====\n")
print(sort(table(subtype), decreasing = TRUE))

cat("\n===== 亚型 -> 主导程序（mean usage 最大）=====\n")
for (s in usub) {
  cat(sprintf("%-14s -> %s\n", s, rownames(mean_usage)[which.max(mean_usage[, s])]))
}

cat("\n===== 程序 -> 主导亚型（frac_top，每程序被哪些亚型当首选）=====\n")
for (j in seq_len(np)) {
  o <- order(frac_top[j, ], decreasing = TRUE)
  top3 <- usub[o][seq_len(min(3, length(o)))]
  cat(sprintf("P%d: %s\n", j,
      paste(sprintf("%s(%.1f%%)", top3, 100 * frac_top[j, top3]), collapse = ", ")))
}

cat("\n===== 每程序 top 基因（programs 列表，前 20）=====\n")
for (j in seq_len(np)) {
  cat(sprintf("P%d: %s\n", j, paste(head(programs[[j]], 20), collapse = ", ")))
}

# ---- 交叉核对：programs[[j]] 是否 = W[,j] 的 top 基因 ----
cat("\n===== 交叉核对（programs 列表 vs W 列 top 基因）=====\n")
for (j in seq_len(np)) {
  w_top <- genes[order(W[, j], decreasing = TRUE)[seq_len(length(programs[[j]]))]]
  cat(sprintf("P%d: programs 与 W 列重叠 %d/%d\n", j,
      length(intersect(programs[[j]], w_top)), length(programs[[j]])))
}

cat("\n===== mean usage 矩阵（程序 × 亚型）=====\n")
print(round(mean_usage, 3))

write.csv(mean_usage, file.path(res_dir, "tables", "program_mean_usage_by_subtype.csv"))
write.csv(frac_top,   file.path(res_dir, "tables", "program_frac_top_by_subtype.csv"))
cat("\n已保存到:", file.path(res_dir, "tables"), "\n")
