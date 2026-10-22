# ============================================================
# 01_build_validation_cohorts.R
# 目的：把验证层生存队列整理成干净 RDS
#   - 训练队列：TCGA-LUAD（RNA-seq）
#   - 验证队列：GSE68465 / GSE72094 / GSE31210（微阵列）
# 只读本地文件，不需要 GEOquery、不需要联网
# 运行：本地 RStudio，路径已按你的实际下载位置填好
# ============================================================

base    <- "."
out_dir <- file.path(base, "validation_cohorts")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 工具函数 1：解析 GEO series_matrix（不依赖 GEOquery）
# 返回 list(expr = 探针×样本矩阵, pdata = 样本元数据 data.frame)
# ============================================================
read_series_matrix <- function(path) {
  cat("读取:", path, "\n")
  con <- gzfile(path, "rt")
  lines <- readLines(con, warn = FALSE)
  close(con)

  tb <- which(lines == "!series_matrix_table_begin")[1]
  te <- which(lines == "!series_matrix_table_end")[1]
  if (is.na(tb) || is.na(te)) stop("未找到表达表格标记：", path)

  # 表达矩阵（表体：表头 + 探针行）
  expr <- read.delim(text = paste(lines[(tb + 1):(te - 1)], collapse = "\n"),
                     header = TRUE, check.names = FALSE, row.names = 1,
                     stringsAsFactors = FALSE)
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"

  # 元数据（! 开头的行，在表之前）
  meta_lines <- lines[seq_len(tb - 1)]
  meta_lines <- meta_lines[grepl("^!", meta_lines)]
  p <- parse_geo_meta(meta_lines, ncol(expr))

  list(expr = expr, pdata = p)
}

# 工具函数 2：把 !Sample_* 元数据行转成 data.frame（行=样本，列=字段）
# characteristics 逐样本按 "key: value" 前缀拆列 —— 位置无关，能抗 GEO 系列矩阵的串列/空行
parse_geo_meta <- function(meta_lines, nsample) {
  meta_lines <- sub("^!", "", meta_lines)
  parts <- strsplit(meta_lines, "\t", fixed = TRUE)

  # 只保留"样本级"字段（值数量 == 样本数）
  keep <- vapply(parts, function(x) length(x) == nsample + 1, logical(1))
  parts <- parts[keep]
  keys  <- vapply(parts, function(x) x[1], character(1))
  vals  <- lapply(parts, function(x) sub('^"', '', sub('"$', '', x[-1])))

  # --- 非 characteristics 字段：整列直接收 ---
  is_char <- grepl("^Sample_characteristics", keys)
  out <- list()
  for (k in unique(keys[!is_char])) {
    if (!nzchar(k)) next
    out[[k]] <- vals[[match(k, keys)]]
  }

  # --- characteristics：逐样本、按每个值自带的 "key:" 前缀拆到正确列（抗串列/空行）---
  char_vals <- vals[is_char]
  sample_list <- replicate(nsample, list(), simplify = FALSE)
  for (vi in char_vals) {
    for (j in seq_along(vi)) {
      v <- vi[j]
      if (is.na(v) || !nzchar(v) || !grepl(":", v, fixed = TRUE)) next
      k <- trimws(strsplit(v, ":", fixed = TRUE)[[1]][1])
      if (!nzchar(k)) next
      val <- trimws(sub("^[^:]+:[[:space:]]*", "", v))
      sample_list[[j]][[k]] <- val
    }
  }
  allkeys <- unique(unlist(lapply(sample_list, names)))
  if (length(allkeys)) {
    char_df <- as.data.frame(
      do.call(rbind, lapply(sample_list, function(x)
        vapply(allkeys, function(k) if (!is.null(x[[k]])) x[[k]] else NA_character_,
               character(1)))),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    out <- c(out, as.list(char_df))
  }

  out <- as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)

  # 样本 ID（GSM）作为行名
  if ("Sample_geo_accession" %in% names(out)) {
    rownames(out) <- out[["Sample_geo_accession"]]
  }
  out
}

# ============================================================
# 1. TCGA-LUAD（训练队列）
# ============================================================
cat("\n==================== TCGA-LUAD（训练） ====================\n")
cnt_file <- file.path(base, "TCGA-LUAD.star_counts (1).tsv.gz")  # 若你改过名，改这里
surv_file <- file.path(base, "TCGA-LUAD.survival.tsv.gz")
clin_file <- file.path(base, "TCGA-LUAD.clinical.tsv.gz")

if (file.exists(cnt_file)) {
  # 表达矩阵：行 Ensembl_ID，列 TCGA barcode；值是 log2(count+1)
  cnt <- read.delim(gzfile(cnt_file), header = TRUE, check.names = FALSE,
                    row.names = 1, stringsAsFactors = FALSE)
  cat("计数维度（基因 × 样本）:", dim(cnt), "\n")

  # 只留原发肿瘤（barcode 第 14-15 位 == "01"；"11" 是正常组织）
  is_tumor <- substr(colnames(cnt), 14, 15) == "01"
  cnt <- cnt[, is_tumor, drop = FALSE]
  cat("保留原发肿瘤样本:", ncol(cnt), "\n")

  # 生存数据（列：sample / OS.time / OS / _PATIENT）
  surv <- read.delim(gzfile(surv_file), header = TRUE, check.names = FALSE,
                     stringsAsFactors = FALSE)
  cat("生存数据样本数:", nrow(surv), " | 列:", paste(colnames(surv), collapse = ", "), "\n")

  # 对齐：barcode 前 15 位（去掉末位 vial 字母）
  key_cnt  <- substr(colnames(cnt), 1, 15)
  key_surv <- substr(surv$sample, 1, 15)
  common <- intersect(key_cnt, key_surv)
  cat("对齐成功样本数:", length(common), "\n")

  m_cnt  <- match(common, key_cnt)
  m_surv <- match(common, key_surv)
  cnt  <- cnt[, m_cnt, drop = FALSE]
  surv <- surv[m_surv, , drop = FALSE]
  rownames(surv) <- common

  # 临床表型（可选，供后面 nomogram / 多因素校正用）
  clin <- NULL
  if (file.exists(clin_file)) {
    clin <- read.delim(gzfile(clin_file), header = TRUE, check.names = FALSE,
                       stringsAsFactors = FALSE)
    cat("临床表型样本数:", nrow(clin), "\n")
  }

  saveRDS(list(counts = cnt, survival = surv, clinical = clin),
          file.path(out_dir, "TCGA-LUAD.rds"))
  cat("已保存: TCGA-LUAD.rds\n")
} else {
  cat("[跳过] 未找到计数文件:", cnt_file, "\n")
}

# ============================================================
# 2. GEO 微阵列队列（表达 + 生存都在 series_matrix 里）
# ============================================================
geos <- list(
  GSE68465 = list(file = "GSE68465/GSE68465_series_matrix.txt.gz",
                  url  = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68465/matrix/GSE68465_series_matrix.txt.gz"),
  GSE72094 = list(file = "GSE72094/GSE72094_series_matrix.txt.gz",
                  url  = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72094/matrix/GSE72094_series_matrix.txt.gz"),
  GSE31210 = list(file = "GSE31210/GSE31210_series_matrix.txt.gz",
                  url  = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31210/matrix/GSE31210_series_matrix.txt.gz")
)

for (g in names(geos)) {
  cat("\n==================== ", g, " ====================\n")
  f <- file.path(base, geos[[g]]$file)
  if (!file.exists(f)) {
    cat("[缺失] 未找到文件，请下载：\n  ", geos[[g]]$url, "\n", sep = "")
    next
  }
  res <- read_series_matrix(f)
  cat("表达矩阵维度（探针 × 样本）:", dim(res$expr), "\n")
  cat("样本元数据字段:", paste(colnames(res$pdata), collapse = ", "), "\n")

  # 自动探测生存字段
  surv_cols <- grep("survival|os\\.|os_|vital|status|death|month|day|time|fustat|futime",
                    colnames(res$pdata), ignore.case = TRUE, value = TRUE)
  if (length(surv_cols)) {
    cat(">>> 疑似生存字段:", paste(surv_cols, collapse = ", "), "\n")
  } else {
    cat(">>> 注意：series_matrix 里【没有】生存字段，需另找临床文件补充 OS。\n")
  }

  saveRDS(res, file.path(out_dir, paste0(g, "_eset.rds")))
  cat("已保存:", paste0(g, "_eset.rds"), "\n")
}

cat("\n全部完成。输出目录:", out_dir, "\n")
