# ============================================================
# 20_bulk_program_scores.R
# 目的：把 scRNA NMF 抽出的程序基因，在 4 个 bulk 队列里算「程序打分」
#   - 跨平台基因映射：TCGA(Ensembl→symbol) / 微阵列(探针→symbol via GPL)
#   - 打分：每个基因跨样本 z-score，再按程序取均值
#   - 生存：时间统一成「天」，事件统一成 0/1
# 依赖：org.Hs.eg.db / AnnotationDbi / GEOquery（已装）
# 运行：本地 RStudio 全选运行（文件都在 D 盘）
# ============================================================

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GEOquery)
})

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
dir.create(val_dir, showWarnings = FALSE)

# ---------- 0. 读程序基因（云上复制来的 CSV，在 xinjihefenxi/ 根目录） ----------
# CSV 是「宽表」：列名 P1~P10，每列 = 该程序的 top-50 基因（按 NMF 权重降序）
prog_df <- read.csv(file.path(base, "program_genes_top50.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
cat("程序基因表：", nrow(prog_df), "行 ×", ncol(prog_df), "列\n")
cat("列名（程序）:", paste(colnames(prog_df), collapse = ", "), "\n")

# 每列各自取该程序的基因列表
prog_genes <- lapply(colnames(prog_df), function(p) {
  g <- as.character(prog_df[[p]])
  unique(g[!is.na(g) & g != ""])
})
names(prog_genes) <- colnames(prog_df)

# 排除 P1（混合细胞周期程序：SPP1+TUBA1B+STMN1 等混杂，非稳定 TAM 功能；且与 P2 高度重叠会共线性）
# 若想保留 P1，注释掉下面这一行
prog_genes <- prog_genes[setdiff(names(prog_genes), "P1")]

cat("程序:", paste(names(prog_genes), collapse = ", "), "\n")
cat("每程序基因数:", paste(sprintf("%s=%d", names(prog_genes), lengths(prog_genes)), collapse = " "), "\n")

# ---------- 工具函数 ----------
# 按 symbol 折叠矩阵（同一 symbol 取均值）；自动丢 NA 分组、强制矩阵
collapse_mean <- function(mat, group) {
  mat  <- as.matrix(mat)
  g    <- as.character(group)
  keep <- !is.na(g)
  mat  <- mat[keep, , drop = FALSE]
  g    <- g[keep]
  gf   <- factor(g)
  s <- rowsum(mat, gf, reorder = TRUE)
  n <- as.vector(table(gf)[rownames(s)])
  s / n
}

# 对单个程序打分：基因跨样本 z-score → 均值
score_program <- function(mat_sym, genes) {
  g <- intersect(genes, rownames(mat_sym))
  if (length(g) < 3) return(list(score = rep(NA_real_, ncol(mat_sym)), n = length(g)))
  z <- t(scale(t(mat_sym[g, , drop = FALSE])))
  list(score = colMeans(z, na.rm = TRUE), n = length(g))
}

# 对全部程序打分，返回 samples×programs 的 data.frame + 覆盖度
score_all <- function(mat_sym, prog_genes) {
  scores <- matrix(NA_real_, nrow = ncol(mat_sym), ncol = length(prog_genes))
  colnames(scores) <- names(prog_genes)
  rownames(scores) <- colnames(mat_sym)
  cov <- numeric(length(prog_genes)); names(cov) <- names(prog_genes)
  for (p in names(prog_genes)) {
    s <- score_program(mat_sym, prog_genes[[p]])
    scores[, p] <- s$score
    cov[p] <- s$n
  }
  list(scores = as.data.frame(scores), coverage = cov)
}

# 事件文本 → 0/1
to_event <- function(x) {
  x <- trimws(as.character(x))
  dead  <- c("dead", "Dead", "1", "yes", "Yes", "TRUE", "True")
  alive <- c("alive", "Alive", "0", "no", "No", "FALSE", "False")
  ev <- rep(NA_real_, length(x))
  ev[x %in% dead]  <- 1
  ev[x %in% alive] <- 0
  ev
}
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# 读平台注释，返回 探针ID -> symbol 的命名向量（1:many 取第一个）
# 带重试：NCBI GEO 经常 SSL 连接重置，多试几次
fetch_gpl <- function(gpl_id, max_tries = 6) {
  soft <- file.path(val_dir, paste0(gpl_id, ".soft.gz"))
  if (file.exists(soft) && file.info(soft)$size < 1000) file.remove(soft)  # 清掉半截/空缓存

  gpl <- NULL
  for (i in 1:max_tries) {
    cat("  [", gpl_id, "] 尝试 ", i, "/", max_tries, " ...\n", sep = "")
    ok <- tryCatch({
      gpl <- getGEO(gpl_id, destdir = val_dir, AnnotGPL = FALSE, getGPL = TRUE)
      TRUE
    }, error = function(e) {
      cat("    失败：", conditionMessage(e), "\n", sep = "")
      FALSE
    })
    if (ok) break
    Sys.sleep(12)
  }
  if (is.null(gpl)) stop("  ", gpl_id, " 重试多次仍下载失败，稍后重跑本脚本或换网络")

  tab <- Table(gpl)
  id_col  <- grep("^ID$", colnames(tab), ignore.case = TRUE)[1]
  sym_col <- grep("symbol", colnames(tab), ignore.case = TRUE)[1]
  if (is.na(sym_col)) sym_col <- grep("gene", colnames(tab), ignore.case = TRUE)[1]
  if (is.na(sym_col)) stop("  找不到 symbol 列，实际列名：",
                           paste(colnames(tab), collapse = ", "))
  id  <- as.character(tab[[id_col]])
  sym <- as.character(tab[[sym_col]])
  sym <- trimws(sub("\\s*///.*$", "", sym))     # "A /// B" 取 A
  out <- sym; names(out) <- id
  out <- out[nzchar(out)]
  cat("  ", length(out), " 个探针有 symbol\n", sep = "")
  out
}

# ---------- 结果容器 ----------
scores_list   <- list()
surv_list     <- list()
coverage_list <- list()

# ============================================================
# 1. TCGA-LUAD（Ensembl → symbol；训练队列）
# ============================================================
cat("\n===== TCGA-LUAD（训练） =====\n")
tcga <- readRDS(file.path(val_dir, "TCGA-LUAD.rds"))
cnt  <- tcga$counts
rn   <- rownames(cnt)
if (grepl("^ENSG", rn[1])) {
  ens <- sub("\\..*$", "", rn)
  sym <- mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first")
  names(sym) <- rn
  sym <- sym[!is.na(sym) & sym != ""]
} else {
  sym <- setNames(rn, rn)
}
cat("Ensembl → symbol 映射成功:", length(sym), " 个基因\n")

keep <- !is.na(sym[rownames(cnt)])          # 丢映射失败的基因
cnt  <- cnt[keep, , drop = FALSE]
grp  <- sym[rownames(cnt)]
tcga_sym <- collapse_mean(cnt, grp)
colnames(tcga_sym) <- substr(colnames(tcga_sym), 1, 15)   # 完整 barcode → 15 位 key，与 survival 行名对齐
cat("TCGA symbol 矩阵（基因 × 样本）:", dim(tcga_sym), "\n")

sc <- score_all(tcga_sym, prog_genes)

# 生存：OS.time 天 + OS 0/1（rownames 已对齐到 15 位 barcode）
surv <- tcga$survival
tcga_surv <- data.frame(time  = to_num(surv$OS.time),
                        event = to_num(surv$OS),
                        row.names = rownames(surv))
common <- intersect(rownames(sc$scores), rownames(tcga_surv))
scores_list[["TCGA-LUAD"]]   <- sc$scores[common, , drop = FALSE]
surv_list[["TCGA-LUAD"]]     <- tcga_surv[common, , drop = FALSE]
coverage_list[["TCGA-LUAD"]] <- sc$coverage
cat("TCGA 样本:", length(common),
    " | 死亡事件:", sum(surv_list[["TCGA-LUAD"]]$event == 1, na.rm = TRUE), "\n")

# ============================================================
# 2-4. 三个 GEO 微阵列队列
# ============================================================
array_cohorts <- list(
  GSE68465 = list(gpl = "GPL96",    time = "months_to_last_contact_or_death", time_unit = "months", event = "vital_status"),
  GSE72094 = list(gpl = "GPL15048", time = "survival_time_in_days",           time_unit = "days",   event = "vital_status"),
  GSE31210 = list(gpl = "GPL570",   time = "days before death/censor",        time_unit = "days",   event = "death")
)

for (g in names(array_cohorts)) {
  cat("\n===== ", g, " =====\n", sep = "")
  f <- file.path(val_dir, paste0(g, "_eset.rds"))
  if (!file.exists(f)) { cat("  [缺失] 没有 ", f, "，跳过\n", sep = ""); next }
  res <- readRDS(f)
  expr <- res$expr
  pdata <- res$pdata
  cat("表达矩阵（探针 × 样本）:", dim(expr), "\n")

  # 探针 → symbol
  probe2sym <- fetch_gpl(array_cohorts[[g]]$gpl)
  keep <- !is.na(probe2sym[rownames(expr)])
  expr <- expr[keep, , drop = FALSE]
  sym  <- probe2sym[rownames(expr)]
  mat_sym <- collapse_mean(expr, sym)
  cat(g, " symbol 矩阵（基因 × 样本）:", dim(mat_sym), "\n", sep = "")

  sc <- score_all(mat_sym, prog_genes)

  # ---- 生存字段（精确列名） ----
  cfg <- array_cohorts[[g]]
  if (!(cfg$time %in% colnames(pdata))) {
    cat("  [警告] 找不到时间列 '", cfg$time, "'，实际列名: ",
        paste(colnames(pdata), collapse = ", "), "\n", sep = "")
    next
  }
  if (!(cfg$event %in% colnames(pdata))) {
    cat("  [警告] 找不到事件列 '", cfg$event, "'，实际列名: ",
        paste(colnames(pdata), collapse = ", "), "\n", sep = "")
    next
  }
  time <- to_num(pdata[[cfg$time]])
  if (cfg$time_unit == "months") time <- time * 30.44   # 统一成「天」
  event <- to_event(pdata[[cfg$event]])

  # GSE31210 有「排除预后分析」标记，剔除
  excl_col <- grep("exclude", colnames(pdata), ignore.case = TRUE, value = TRUE)
  excl_flag <- rep(FALSE, nrow(pdata))
  if (length(excl_col)) {
    excl_val <- trimws(tolower(as.character(pdata[[excl_col[1]]])))
    excl_flag <- grepl("exclude|^1$|^yes$|true", excl_val)
    cat("  排除列 '", excl_col[1], "' 取值分布: ",
        paste(names(table(excl_val)), collapse = ", "), "\n", sep = "")
  }

  surv_df <- data.frame(time = time, event = event, row.names = rownames(pdata))
  surv_df <- surv_df[!excl_flag & !is.na(time) & !is.na(event), , drop = FALSE]

  common <- intersect(rownames(sc$scores), rownames(surv_df))
  scores_list[[g]]   <- sc$scores[common, , drop = FALSE]
  surv_list[[g]]     <- surv_df[common, , drop = FALSE]
  coverage_list[[g]] <- sc$coverage
  cat(g, " 样本:", length(common),
      " | 死亡事件:", sum(surv_list[[g]]$event == 1, na.rm = TRUE), "\n", sep = "")
}

# ============================================================
# 汇总 + 保存
# ============================================================
cat("\n===== 汇总 =====\n")
for (g in names(scores_list)) {
  cat(sprintf("%-12s 样本 %4d | 事件 %4d | 覆盖度: %s\n",
              g, nrow(scores_list[[g]]),
              sum(surv_list[[g]]$event == 1, na.rm = TRUE),
              paste(sprintf("%s=%d", names(coverage_list[[g]]), coverage_list[[g]]),
                    collapse = " ")))
}

out <- list(scores = scores_list, surv = surv_list, coverage = coverage_list)
saveRDS(out, file.path(val_dir, "bulk_program_scores.rds"))
cat("\n已保存:", file.path(val_dir, "bulk_program_scores.rds"), "\n")
