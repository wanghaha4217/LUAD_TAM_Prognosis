# ============================================================
# 41_nmf_rank_selection.R  (v2 — 本地可跑、内存安全版)
# NMF 秩(程序数 k)选择：cophenetic 系数 + dispersion + silhouette
# 目的：为「k=10」提供稳定性证据（Supplementary 证据，不是唯一依据）
#
# 依赖：NMF(CRAN，会自动带 Biobase/BiocGenerics) / Seurat / ggplot2
# 安装：install.packages("NMF")   # 若 Biobase 装不上，先 BiocManager::install("Biobase")
#
# ⚠️ 为什么抽 5000 细胞（2026-09-27 排雷结论，务必保留）：
#   全 38211 细胞跑不通 —— cophenetic/dispersion/silhouette 需要 n×n 共识矩阵
#   (38211² = 14.6 亿元素 ≈ 117 GB) + hclust 层次聚类 O(n³)，光 k=5 一个就跑 3 天，
#   连云端 256 GB 都 OOM。
#   修复：抽 5000 细胞 → 共识矩阵降到 5000² = 2500 万(≈200 MB)、hclust 秒级。
#   图注必须诚实写「在 5000 细胞随机子集上计算」。
#
# ⚠️ 内存安全要点（本地 31.6 GB 跑得动的关键）：
#   ① 先抽 5000 细胞「索引」，再取 data 槽 —— 从不物化 2000×38211 的全量稠密矩阵
#   ② 取完 X 后立刻 rm(so); gc() 删掉 Seurat 对象（约 3-6 GB），NMF 前不留大对象
#   ③ NMF 默认 2 核并行即最优，别去 force cores（20 核折腾不出，反而退化单核）
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(NMF)
  library(ggplot2)
})

# ---- 0) 参数 ----
INPUT_RDS  <- "./GSE243013_LUAD_TAM_clustered.rds"  # 本地完整对象(579 MB)
OUT_TAB    <- "./results/tables"
OUT_FIG    <- "./results/figures"

SUBSAMPLE_N_CELLS <- 5000      # 抽 5000 细胞子集（见头部说明，勿改回全量）
kvec  <- 5:13                  # 待评估的秩范围
nrun  <- 5                     # 每个 k 跑几次；正式投稿图可升到 20~30（时间线性增加）
seed  <- 123456
nfeatures <- 2000              # 必须与原 NMF 的 2000 HVG 一致

# ---- 1) 载入对象 → 先抽细胞索引 → 只取 2000 HVG × 5000 细胞 ----
so <- readRDS(INPUT_RDS)
cat("已载入对象；内存占用(对象) ≈", format(object.size(so), units = "auto"), "\n")

hv <- VariableFeatures(so)
stopifnot(length(hv) == nfeatures)   # 必须与原 NMF 的 2000 HVG 一致

set.seed(seed)
cells <- sample(colnames(so), SUBSAMPLE_N_CELLS)      # ① 先抽 5000 细胞索引

X <- GetAssayData(so, layer = "data")[hv, cells]      # ② 只取 2000×5000 = 80 MB（Seurat v5 用 layer，slot 已废弃）
rm(so); gc()                                          # ③ 立刻删大对象，NMF 前不留
cat("抽取子集后 X:", nrow(X), "基因 ×", ncol(X), "细胞\n")
gc()

X <- as.matrix(X)
storage.mode(X) <- "double"

# 抽样后可能有 HVG 在 5000 子集里全为 0（该基因只在少数细胞表达，随机抽到的 5000 细胞里没出现）。
# NMF 不接受全零行 → 剔除并如实报告（原 NMF 用全量 38211 细胞，此处是子集，属正常）。
zero_rows <- rowSums(X) == 0
if (any(zero_rows)) {
  cat("抽样后", sum(zero_rows), "个 HVG 在 5000 细胞子集里全为 0，已剔除：",
      paste(head(hv[zero_rows]), collapse = ", "),
      if (sum(zero_rows) > 6) " ..." else "", "\n")
  X <- X[!zero_rows, , drop = FALSE]
}
stopifnot(nrow(X) > 0)

cat("X 尺寸:", nrow(X), "×", ncol(X),
    "| 最小值 =", min(X), "| 大小 =", format(object.size(X), units = "auto"), "\n")
stopifnot(min(X) >= 0)   # NMF 硬前提：非负

# ---- 2) 秩调查(自动算 cophenetic/dispersion/silhouette/rss/sparseness) ----
cat("开始 nmfEstimateRank: k =", paste(kvec, collapse = ","),
    ", nrun =", nrun, "（默认 2 核并行，请勿 force cores）\n")
est <- nmfEstimateRank(X, range = kvec, method = "brunet",
                       nrun = nrun, seed = seed, .options = "v")

# ---- 3) 数值表 ----
s <- summary(est)
print(s)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
write.csv(s, file.path(OUT_TAB, "nmf_rank_survey.csv"), row.names = FALSE)

# ---- 4) 默认图(base，一次出全部指标) ----
png(file.path(OUT_FIG, "nmf_rank_survey.png"), width = 2000, height = 1200, res = 150)
plot(est)
dev.off()

# ---- 5) 论文图注（写进 txt 备用，投稿时直接贴） ----
caption_en <- paste0(
  "NMF stability metrics (cophenetic coefficient, dispersion, and silhouette) ",
  "computed on a random subset of 5,000 macrophages (n = ", SUBSAMPLE_N_CELLS,
  "; seed ", seed, "). The n x n consensus matrix is computationally infeasible ",
  "across all 38,211 cells. For each k, ", nrun,
  " independent NMF runs were performed (brunet algorithm)."
)
writeLines(caption_en, file.path(OUT_TAB, "nmf_rank_caption.txt"))
cat("\n完成。表:", file.path(OUT_TAB, "nmf_rank_survey.csv"),
    "\n图:", file.path(OUT_FIG, "nmf_rank_survey.png"),
    "\n图注:", file.path(OUT_TAB, "nmf_rank_caption.txt"), "\n")
cat("图注(英文，投稿用)：\n", caption_en, "\n")

# ============================================================
# 备选：手动跑每个 k，自己算 cophenetic + dispersion(可控、可自定义 ggplot)
# 若 nmfEstimateRank 输出/耗时不合意，用下面这段替代第 2~4 步。
# 注意：此路径会显式构建共识矩阵 C(n×n)，仅在 5000 子集下可行。
# ============================================================
# res <- data.frame(k = kvec, cophenetic = NA_real_, dispersion = NA_real_)
# for (k in kvec) {
#   cat("k =", k, " ...\n")
#   fit <- nmf(X, rank = k, method = "brunet", nrun = nrun, seed = seed, .options = "v")
#   C <- consensus(fit)                    # 共识(连接)矩阵，细胞 × 细胞
#   res$cophenetic[res$k == k] <- cophcor(C)          # Brunet cophenetic 相关系数
#   res$dispersion[res$k == k] <- dispersion(C)       # Kim & Park dispersion
# }
# print(res)
# write.csv(res, file.path(OUT_TAB, "nmf_coph_disp.csv"), row.names = FALSE)
#
# g <- ggplot(res, aes(k)) +
#   geom_line(aes(y = cophenetic),  color = "steelblue") +
#   geom_point(aes(y = cophenetic), color = "steelblue") +
#   labs(x = "Rank (k)", y = "Cophenetic coefficient",
#        title = "NMF stability: cophenetic coefficient vs rank") +
#   theme_bw()
# ggsave(file.path(OUT_FIG, "nmf_cophenetic.png"), g, width = 7, height = 5, dpi = 300)
