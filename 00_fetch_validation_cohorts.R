# ============================================================
# 00_fetch_validation_cohorts.R
# 目的：整理"验证层"生存队列 —— 训练 TCGA-LUAD + 独立验证 GSE68465/GSE72094/GSE31210
# 输出：每个队列一个干净 RDS（表达矩阵 + 生存数据对齐）
# 运行：本地 RStudio（需联网访问 GEO / 已下好的本地文件）
# ============================================================

# ---------- 依赖 ----------
# 若 GEOquery 未装：BiocManager::install("GEOquery")
suppressPackageStartupMessages(library(GEOquery))

out_dir <- "./validation_cohorts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 第一部分：三个 GEO 微阵列队列（表达 + 生存都在 series_matrix 里）
# ============================================================
geos <- c("GSE68465", "GSE72094", "GSE31210")

for (g in geos) {
  cat("\n==================== ", g, " ====================\n")

  # 下载 series_matrix（若已手动下载到 out_dir，会直接复用本地文件）
  gset <- getGEO(g, GSEMatrix = TRUE, destdir = out_dir, getGPL = FALSE)
  eset <- gset[[1]]

  # 表达矩阵（探针 × 样本）
  expr <- exprs(eset)
  cat("表达矩阵维度（探针 × 样本）:", dim(expr), "\n")

  # 样本元数据（含生存）
  p <- pData(eset)
  cat("样本数:", nrow(p), "\n")
  cat("--- 所有临床字段名 ---\n")
  print(colnames(p))

  # 自动探测生存相关字段（不依赖确切命名）
  surv_cols <- grep("survival|os\\.|os_|vital|status|death|months|days|time|fustat|futime",
                    colnames(p), ignore.case = TRUE, value = TRUE)
  cat("--- 疑似生存字段:", paste(surv_cols, collapse = " | "), " ---\n")

  saveRDS(list(expr = expr, pdata = p), file.path(out_dir, paste0(g, "_eset.rds")))
  cat("已保存:", file.path(out_dir, paste0(g, "_eset.rds")), "\n")
}

# ============================================================
# 第二部分：TCGA-LUAD（Xena 直接下载的两个文件）
# 直接 URL（粘贴到浏览器/迅雷即可，无需在 xenabrowser 网站点选）：
#   1) 表达矩阵: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.htseq_counts.tsv.gz
#   2) 生存数据: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.survival.tsv
# 下好后把本地路径填到下面两个变量再跑。
# ============================================================
xena_counts_file <- ""   # 例如 "D:/.../TCGA-LUAD.htseq_counts.tsv.gz"
xena_surv_file   <- ""   # 例如 "D:/.../TCGA-LUAD.survival.tsv"

if (nzchar(xena_counts_file) && nzchar(xena_surv_file)) {
  cat("\n==================== TCGA-LUAD ====================\n")

  # 读表达矩阵：第一列 Ensembl_ID，列名是 TCGA barcode
  cnt <- read.delim(gzfile(xena_counts_file), header = TRUE, check.names = FALSE,
                    row.names = 1, stringsAsFactors = FALSE)
  # Xena 的 htseq_counts 末尾 5 行是汇总/QC 行，去掉
  cnt <- cnt[1:(nrow(cnt) - 5), , drop = FALSE]
  cnt_mat <- as.matrix(cnt)
  storage.mode(cnt_mat) <- "numeric"
  cat("表达矩阵维度（基因 × 样本，去汇总行后）:", dim(cnt_mat), "\n")

  # Xena 的 htseq_counts 是 log2(x+1) 变换的（值带小数）；若你下到的是原始整数则自动跳过
  if (!all(cnt_mat == round(cnt_mat))) {
    cat("检测到 log2(x+1) 变换，反变换回原始计数 ...\n")
    cnt_mat <- round(2^cnt_mat - 1)
  }

  # 只保留原发肿瘤（barcode 第 14-15 位 == "01"）
  is_tumor <- substr(colnames(cnt_mat), 14, 15) == "01"
  cnt_mat <- cnt_mat[, is_tumor, drop = FALSE]
  cat("原发肿瘤样本数:", ncol(cnt_mat), "\n")

  # 读生存文件（列含 sample / OS / OS.time）
  surv <- read.delim(xena_surv_file, header = TRUE, check.names = FALSE,
                     stringsAsFactors = FALSE)
  cat("生存文件列名:", paste(colnames(surv), collapse = ", "), "\n")

  # 对齐：计数列名取前 15 位（TCGA-XX-XXXX-01）与生存 sample 对齐
  key_cnt  <- substr(colnames(cnt_mat), 1, 15)
  key_surv <- substr(surv$sample, 1, 15)
  common <- intersect(key_cnt, key_surv)
  cat("对齐成功样本数:", length(common), "\n")

  m_cnt  <- match(common, key_cnt)
  m_surv <- match(common, key_surv)
  cnt_mat <- cnt_mat[, m_cnt, drop = FALSE]
  surv    <- surv[m_surv, , drop = FALSE]
  rownames(surv) <- common

  saveRDS(list(counts = cnt_mat, survival = surv), file.path(out_dir, "TCGA-LUAD.rds"))
  cat("已保存:", file.path(out_dir, "TCGA-LUAD.rds"), "\n")
} else {
  cat("\n[TCGA-LUAD] 尚未指定 Xena 文件路径，跳过（下好后填 xena_counts_file / xena_surv_file 再跑）。\n")
}

cat("\n全部完成。输出目录:", out_dir, "\n")
