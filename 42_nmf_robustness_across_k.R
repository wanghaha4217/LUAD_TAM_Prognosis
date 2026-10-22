# ============================================================
# 42_nmf_robustness_across_k.R  (本地 Windows RStudio 运行)
# 目的：方案 A —— 回答「k 取其他值，签名在 3 独立队列是否更好/更差」
#   对 k ∈ {8, 9, 10, 11} 重跑 NMF 分解（围绕 10，直接回答「为什么是 10 不是 9/11」），四个 k 用完全相同的：
#     输入(全量 38211 细胞 × 2000 HVG, data 槽) / 方法(brunet) / seed / nrun
#   使差异只来自 k，不来自随机或预处理。
#   只做「分解」拿 W/H，不碰 cophenetic/silhouette（那才是 9 小时的来源）。
#
# 依赖：NMF / Seurat（均已装）
# 输入：./GSE243013_LUAD_TAM_clustered.rds
# 输出：CSV 直接写到 xinjihefenxi/ 根目录（脚本 54 会从这里读，免手动复制）；
#       中间 .rds 写到 xinjihefenxi/robustness/
#
# ⚠️ 重要说明：
#   - 原论文的 k=10 是交互式跑的、seed 未存档；本脚本的 k=10 是「全新重跑」，
#     可能与论文原 k=10 的程序略有差异。对「跨 k 比较」来说，本脚本这套才是公平基准。
#   - 关键：rm(so); gc() 必须在 nmf 前 —— 否则 NMF 并行会把 2.7G 对象
#     复制进每个 worker（OOM 的根因，见记忆）。
#   - 本地 Windows 默认「顺序」跑（没注册并行后端），内存安全但慢；
#     若你的 RStudio 里已注册过并行后端，脚本会尝试切回顺序，避免 PSOCK 复制 600MB 矩阵 OOM。
#   - 运行时长为数小时级别（4 个 k × 每次多 start）。可先把 nrun 设 3 试速，再决定是否加到 10。
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(NMF)
})

# ---- 0) 路径 + 参数（如需调整只改这里）----
data_dir <- "."               # 下载的 .rds 所在
so_file  <- "GSE243013_LUAD_TAM_clustered.rds"                  # 2000 HVG 的巨噬对象（非 _seurat）
proj_dir <- "."       # 输出根（脚本 54 从这里读 CSV）
out_dir  <- file.path(proj_dir, "robustness")                   # 中间 .rds 放这里
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

kvec      <- c(8, 9, 10, 11)   # 待评估的 k（围绕 10）
nrun      <- 5                 # 每个 k 的独立 NMF 次数（越大越稳，越慢；本地慢，先用 5，有耐心可改 10）
seed      <- 123456            # 四 k 统一，保证公平
nfeatures <- 2000              # 必须与原 NMF 一致

# ---- 1) 载入完整对象，取非负 data 矩阵（2000 HVG × 38211）----
so <- readRDS(file.path(data_dir, so_file))
cat("对象大小:", format(object.size(so), units = "auto"),
    "| 细胞:", ncol(so), "| 基因:", nrow(so), "\n")

hv <- VariableFeatures(so)
cat("VariableFeatures 数量:", length(hv), "（应 =", nfeatures, "）\n")
stopifnot(length(hv) == nfeatures)   # 与原 NMF 的 2000 HVG 一致；若这里报错，换 so_file = _seurat 再试

X <- GetAssayData(so, layer = "data")[hv, ]   # log1p 归一化数据，非负
rm(so); gc()                                  # ★ 关键：fork 前删大对象，防 OOM

# 确保顺序执行（Windows 下并行 PSOCK 会把 600MB 矩阵复制进每个 worker，易 OOM）
if (requireNamespace("foreach", quietly = TRUE)) {
  foreach::registerDoSEQ()
}
if (exists("setupParallel", where = asNamespace("NMF"))) {
  tryCatch(get("setupParallel", envir = asNamespace("NMF"))(FALSE), error = function(e) NULL)
}

X <- as.matrix(X)
storage.mode(X) <- "double"
stopifnot(min(X) >= 0)                        # NMF 硬前提：非负
cat("输入矩阵:", nrow(X), "基因 ×", ncol(X), "细胞 | min =", min(X), "\n")

# ---- 2) 逐 k 跑 NMF，导出 W/H + top-50 基因表 ----
for (k in kvec) {
  cat("\n===== k =", k, "（nrun =", nrun, "）=====\n")
  fit <- nmf(X, rank = k, method = "brunet", nrun = nrun,
             seed = seed, .options = "v")

  W <- basis(fit)                 # 2000 基因 × k（每列 = 一个程序的基因权重）
  H <- coef(fit)                  # k × 38211（每行 = 一个程序在各细胞的活性）
  rownames(W) <- hv

  # top-50 基因/程序（按 W 权重降序），宽表：列名 P1..Pk
  prog_genes <- apply(W, 2, function(w) names(sort(w, decreasing = TRUE))[1:50])
  df <- as.data.frame(prog_genes, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- paste0("P", seq_len(k))

  write.csv(df, file.path(proj_dir, paste0("program_genes_top50_k", k, ".csv")),
            row.names = FALSE)
  saveRDS(list(W = W, H = H, k = k, seed = seed, nrun = nrun),
          file.path(out_dir, paste0("nmf_k", k, "_fresh.rds")))

  cat("已存:", file.path(proj_dir, paste0("program_genes_top50_k", k, ".csv")),
      " / ", file.path(out_dir, paste0("nmf_k", k, "_fresh.rds")), "\n")
}

cat("\n全部完成。\n")
cat("CSV 已在:", proj_dir, "（脚本 54 直接读这里）\n")
cat("中间 rds 已在:", out_dir, "\n")
cat("下一步：本地跑 scripts/54_cross_k_validate_local.R，看「结论」输出。\n")
