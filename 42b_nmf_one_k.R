# ============================================================
# 42b_nmf_one_k.R  (本地/云端均可，一次只跑一个 k；手动 nrun=1×nrun，内存 ~2G)
# 目的：方案 A —— 把脚本 42 的「4 个 k 循环」拆成「一次一个 k」。
#   每个 k 是完全独立的 NMF 分解，彼此无依赖，拆开跑的好处：
#     - 单个 k 跑完即存盘，某一个 k 失败/崩溃不影响其他 k；
#     - 每次只占一份内存，不容易 OOM；
#     - 可以 k=8 跑完先看一眼，再改 K 跑 k=9 ……
#
# 用法（四步，只改第 32 行的 K）：
#   ① K <- 8   → 全选运行 → 得 program_genes_top50_k8.csv
#   ② K <- 9   → 全选运行 → 得 program_genes_top50_k9.csv
#   ③ K <- 10  → 全选运行 → 得 program_genes_top50_k10.csv
#   ④ K <- 11  → 全选运行 → 得 program_genes_top50_k11.csv
#   ⑤ 四个都跑完后，再跑一次 scripts/54_cross_k_validate_local.R（读全部 CSV 做打分+验证，无 NMF，快）
#
# ⚠️ 关键：本脚本「手动跑 nrun 次 nrun=1、取 objective 最小者」，等价于 nmf(nrun=nrun) 的
#   「取最优」，但完全不触发 n×n 共识矩阵(38211²≈11.7Gb×多份)——那才是本地 OOM、
#   云端飙 100G、且 consensus 阶段要几小时的元凶。手动版内存 ~2G，本地/云端都能跑。
#
# 依赖：NMF / Seurat（均已装）
# 输入：GSE243013_LUAD_TAM_clustered.rds（2000 HVG 的巨噬对象，非 _seurat）
# 输出：CSV 写到 xinjihefenxi/ 根目录（脚本 54 从这里读）；中间 .rds 写到 robustness/
#
# ⚠️ 重要说明：
#   - 原论文的 k=10 是交互式跑的、seed 未存档；本脚本的 k=10 是「全新重跑」，
#     可能与论文原 k=10 的程序略有差异。对「跨 k 比较」来说，本脚本这套才是公平基准。
#   - 关键：rm(so); gc() 必须在 nmf 前 —— 否则 NMF 并行会把大对象复制进每个 worker（OOM 根因）。
#   - 本地 Windows 默认「顺序」跑；若 RStudio 已注册过并行后端，脚本会切回顺序，避免 PSOCK 复制 600MB 矩阵 OOM。
# ============================================================

# 云上 openblas-pthread 多线程会与 NMF 的顺序计算互相抢，反而慢（记忆里的排雷经验）
Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1)

suppressPackageStartupMessages({
  library(Seurat)
  library(NMF)
})

# ---- 0) 路径 + 参数（只改这里）----
K         <- 8                    # ★ 每次只改这一个数字：8 / 9 / 10 / 11
nrun      <- 5                    # 手动跑 5 次 nrun=1 取最优（不触发共识矩阵，内存 ~2G）
seed      <- 123456               # 四 k 统一，保证公平（改 K 时不要改这个）
nfeatures <- 2000                 # 必须与原 NMF 一致

data_dir <- "."   # 下载的 .rds 所在（若你的在 xinjihefenxi 根目录，改这里）
so_file  <- "GSE243013_LUAD_TAM_clustered.rds"     # 2000 HVG 的巨噬对象（非 _seurat）
proj_dir <- "."   # 输出根（脚本 54 从这里读 CSV）
out_dir  <- file.path(proj_dir, "robustness")               # 中间 .rds 放这里
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1) 载入完整对象，取非负 data 矩阵（2000 HVG × 38211）----
so <- readRDS(file.path(data_dir, so_file))
cat("对象大小:", format(object.size(so), units = "auto"),
    "| 细胞:", ncol(so), "| 基因:", nrow(so), "\n")

hv <- VariableFeatures(so)
cat("VariableFeatures 数量:", length(hv), "（应 =", nfeatures, "）\n")
stopifnot(length(hv) == nfeatures)   # 与原 NMF 的 2000 HVG 一致；若这里报错，换 so_file = _seurat 再试

X <- GetAssayData(so, layer = "data")[hv, ]   # 全基因 sparse → 只取 2000 HVG 行，之后 as.matrix 才 611MB
rm(so); gc()                                  # ★ 关键：NMF 前删大对象，防 OOM

# ★ 强制顺序执行（Windows 下并行 PSOCK 会把 600MB 矩阵复制进每个 worker，易 OOM）
options(mc.cores = 1)                                  # 老 parallel 后端单核
if (requireNamespace("foreach", quietly = TRUE)) {
  foreach::registerDoSEQ()                             # foreach 退化为顺序
}
nmf.options(pbackend = NULL, cores = 1)                # 关掉 NMF 自身并行（清上一脚本残留）
# 跑之前先自查后端，若仍显示非顺序，立刻停手别撞 OOM
cat("NMF 后端检查 → pbackend:", nmf.getOption("pbackend"),
    "| cores:", nmf.getOption("cores"), "（应为 NULL / 1）\n")

tmat <- Sys.time()
X <- as.matrix(X)
storage.mode(X) <- "double"
cat("as.matrix(2000×38211) 耗时 =",
    round(as.numeric(difftime(Sys.time(), tmat, units = "mins")), 2), "min\n")
zero_rows <- rowSums(X) == 0                  # 剔除全零行（NMF 拒绝 null row）
if (any(zero_rows)) {
  cat("剔除", sum(zero_rows), "个全零 HVG\n")
  X <- X[!zero_rows, , drop = FALSE]
  hv <- rownames(X)                           # hv 同步到实际保留的基因
}
stopifnot(min(X) >= 0)                        # NMF 硬前提：非负
cat("输入矩阵:", nrow(X), "基因 ×", ncol(X), "细胞 | min =", min(X), "\n")

# ---- 2) 单个 k 跑 NMF（手动 nrun 次 × nrun=1，取 objective 最小者），导出 W/H + top-50 基因表 ----
# 为什么不用 nmf(nrun=5)：nrun>1 会算 n×n 共识矩阵(38211²≈11.7Gb × 多份)，
#   内存飙到 ~100G、consensus 阶段还要几小时。手动跑 nrun 次 nrun=1、各取 objective
#   最小的那一次，等价于「nrun 次取最优」，但完全不碰共识矩阵，内存 ~2G。
cat("\n===== k =", K, "（手动", nrun, "次 × nrun=1，取 objective 最小者）=====\n")
t0 <- Sys.time()
best_fit <- NULL; best_obj <- Inf
for (i in seq_len(nrun)) {
  cat("  第", i, "/", nrun, "次分解中...\n", sep = "")
  ti <- Sys.time()
  fit_i <- nmf(X, rank = K, method = "brunet", nrun = 1,
               seed = seed + i - 1, .options = "v")
  dt_i <- difftime(Sys.time(), ti, units = "mins")
  obj_i <- as.numeric(fit_i@residuals)   # objective（brunet 的 KL 散度），越小越好
  cat("    objective =", signif(obj_i, 6),
      "| 单次耗时 =", round(as.numeric(dt_i), 2), "min\n")
  if (obj_i < best_obj) { best_obj <- obj_i; best_fit <- fit_i }
}
fit <- best_fit
dt_tot <- difftime(Sys.time(), t0, units = "mins")
cat("  选定最优 objective =", signif(best_obj, 6), "\n")
cat("  k =", K, "总计 =", round(as.numeric(dt_tot), 2), "min（", nrun,
    " 次，平均 ", round(as.numeric(dt_tot) / nrun, 2), " min/次）\n")

W <- basis(fit)                 # 2000 基因 × K（每列 = 一个程序的基因权重）
H <- coef(fit)                  # K × 38211（每行 = 一个程序在各细胞的活性）
rownames(W) <- hv

# top-50 基因/程序（按 W 权重降序），宽表：列名 P1..PK
prog_genes <- apply(W, 2, function(w) names(sort(w, decreasing = TRUE))[1:50])
df <- as.data.frame(prog_genes, stringsAsFactors = FALSE, check.names = FALSE)
colnames(df) <- paste0("P", seq_len(K))

write.csv(df, file.path(proj_dir, paste0("program_genes_top50_k", K, ".csv")),
          row.names = FALSE)
saveRDS(list(W = W, H = H, k = K, seed = seed, nrun = nrun),
        file.path(out_dir, paste0("nmf_k", K, "_fresh.rds")))

cat("\nk =", K, "完成。\n")
cat("CSV 已存:", file.path(proj_dir, paste0("program_genes_top50_k", K, ".csv")), "\n")
cat("rds 已存:", file.path(out_dir, paste0("nmf_k", K, "_fresh.rds")), "\n")
cat("下一步：把开头参数区的 K 改成下一个值再跑；四个都跑完后，跑 scripts/54_cross_k_validate_local.R。\n")
