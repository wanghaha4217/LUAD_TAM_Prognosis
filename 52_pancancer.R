# ============================================================
# 52_pancancer.R
# ⑨ 泛癌验证：P8（缺氧）程序打分 在 TCGA 33 癌种的预后一致性
#   口径与 20_bulk_program_scores.R 完全一致：program_genes_top50.csv 的
#   top-50 基因 → 每基因跨样本 z-score → 取均值；P1 排除。
# 数据：UCSC Xena PanCanAtlas（S3 直链，首次下载 ~1.5-2GB，本地缓存）
#   - 表达：EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena（log2(norm+1)，HUGO symbol）
#   - 生存：Survival_SupplementalTable_S1_20171025_xena_sp（OS/OS.time 天）
#   - 表型：TCGA_phenotype_denseDataOnlyDownload.tsv（sample → acronym 癌种缩写）
# 输出：tables/pancancer_program_hr.csv（程序 × 癌种 HR）
#       tables/pancancer_p8_summary.csv（P8 各癌种，排序）
#       figures/pancancer_p8_forest.png + pancancer_program_heatmap.png
# 依赖：survival、data.table、ggplot2
# 运行：本地 RStudio（首次运行需联网下载）
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
})

base    <- "."
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures")
xena_dir <- file.path(base, "xena_pancan")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(xena_dir, recursive = TRUE, showWarnings = FALSE)

has_dt <- requireNamespace("data.table", quietly = TRUE)

# ---------- 0) 下载（本地缓存；gzip 校验 + curl 带浏览器头兜底）----------
# 注：expression 文件名含 "++"，URL 里须写成 %2B%2B。
#     S3 直链已失效（403），只走 pancanatlas.xenahubs.net；若 R 内置下载被拦，
#     用 curl 带 Referer + User-Agent（常见 403 解法）。
# 判断下载文件是否有效：只按大小判断即可（403 错误页只有 ~110 字节，会被挡掉）。
# 不再强制 gzip magic —— 浏览器下下来的 survival 是未压缩纯文本。
is_valid_file <- function(path, min_size = 1000) {
  file.exists(path) && file.info(path)$size >= min_size
}

try_download <- function(urls, dest, min_size = 1000) {
  if (is_valid_file(dest, min_size)) { cat("[已有] 跳过:", basename(dest), "\n"); return(invisible(TRUE)) }
  if (file.exists(dest)) { file.remove(dest); cat("[清除损坏缓存]", basename(dest), "\n") }
  for (u in urls) {
    for (m in c("auto", "libcurl")) {
      ok <- tryCatch({
        download.file(u, dest, mode = "wb", method = m, quiet = TRUE)
        is_valid_file(dest, min_size)
      }, error = function(e) FALSE)
      if (ok) { cat("[下载成功]", basename(dest), "\n"); return(invisible(TRUE)) }
      if (file.exists(dest)) file.remove(dest)
    }
  }
  if (nzchar(Sys.which("curl"))) {
    for (u in urls) {
      cat("[curl] 尝试:", u, "\n")
      r <- system2("curl", c("-f", "-L", "-A", "Mozilla/5.0",
                             "-e", "https://xenabrowser.net/",
                             "--retry", "2", "-o", dest, u))
      if (isTRUE(r == 0) && is_valid_file(dest, min_size)) {
        cat("[下载成功 via curl]", basename(dest), "\n"); return(invisible(TRUE))
      }
      if (file.exists(dest)) file.remove(dest)
    }
  } else {
    cat("（系统没有 curl，已跳过 curl 兜底）\n")
  }
  stop("下载失败：", basename(dest),
       "。把上面 [curl] 那几行的输出/报错贴给我。若确认被墙，改走 GDC/TCGAbiolinks。")
}

hub <- "https://pancanatlas.xenahubs.net/download/"
paths <- c(expr  = file.path(xena_dir, "expr"),
           surv  = file.path(xena_dir, "surv"),
           pheno = file.path(xena_dir, "pheno"))

try_download(paste0(hub, "EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"),
             paths["expr"], min_size = 100e6)
try_download(paste0(hub, "Survival_SupplementalTable_S1_20171025_xena_sp.gz"),
             paths["surv"], min_size = 1000)
try_download(paste0(hub, "TCGA_phenotype_denseDataOnlyDownload.tsv.gz"),
             paths["pheno"], min_size = 1000)

# ---------- 1) 程序基因表（P1 排除）----------
prog_df <- read.csv(file.path(base, "program_genes_top50.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
prog_genes <- lapply(colnames(prog_df), function(p) {
  g <- as.character(prog_df[[p]]); unique(g[!is.na(g) & g != ""])
})
names(prog_genes) <- colnames(prog_df)
prog_genes <- prog_genes[setdiff(names(prog_genes), "P1")]
cat("程序（P1 已排除）:", paste(names(prog_genes), collapse = ", "), "\n")

# ---------- 2) 表达矩阵 ----------
cat("\n读取表达矩阵（大文件，稍候）...\n")
if (has_dt) {
  expr <- data.table::fread(paths["expr"], header = TRUE, sep = "\t", showProgress = TRUE)
  gene_names <- as.character(expr[[1]])
  samp_names <- colnames(expr)[-1]
  expr_mat <- as.matrix(expr[, -1, with = FALSE])
  rownames(expr_mat) <- gene_names; colnames(expr_mat) <- samp_names
  rm(expr); gc()
} else {
  expr <- read.delim(gzfile(paths["expr"]), header = TRUE, check.names = FALSE,
                     stringsAsFactors = FALSE, sep = "\t")
  gene_names <- expr[[1]]
  expr_mat <- as.matrix(expr[, -1]); rownames(expr_mat) <- gene_names
  rm(expr); gc()
}
cat("表达矩阵:", nrow(expr_mat), "基因 ×", ncol(expr_mat), "样本\n")

# ---------- 3) 生存 + 表型 ----------
# 浏览器下载的 survival 是未压缩纯文本，phenotype 是 .gz；自动识别再读
read_auto <- function(path, ...) {
  con <- file(path, "rb"); h <- readBin(con, "raw", n = 2); close(con)
  if (identical(as.integer(h), c(31L, 139L))) read.delim(gzfile(path), ...) else read.delim(path, ...)
}
surv <- read_auto(paths["surv"], header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE, sep = "\t")
pheno <- read_auto(paths["pheno"], header = TRUE, check.names = FALSE,
                   stringsAsFactors = FALSE, sep = "\t", quote = "")
cat("生存文件列:", paste(head(colnames(surv)), collapse = ", "), " ...\n")
cat("表型文件列:", paste(head(colnames(pheno)), collapse = ", "), " ...\n")

# 表型：找 sample ID 列 / 癌种缩写列 / 样本类型列
sid_col <- grep("sampleID|^sample$|barcode", colnames(pheno), ignore.case = TRUE, value = TRUE)[1]
acr_col <- grep("acronym", colnames(pheno), ignore.case = TRUE, value = TRUE)[1]
if (is.na(acr_col)) acr_col <- grep("primary_disease|^disease$|cancer_type|project_id|^project$",
                                     colnames(pheno), ignore.case = TRUE, value = TRUE)[1]
typ_col <- grep("^sample_type$", colnames(pheno), ignore.case = TRUE, value = TRUE)[1]
if (is.na(typ_col)) typ_col <- grep("sample_type_id", colnames(pheno), ignore.case = TRUE, value = TRUE)[1]
if (is.na(sid_col) || is.na(acr_col)) {
  stop("表型文件找不到 sample/acronym 列，实际列名：", paste(colnames(pheno), collapse = ", "))
}
cat("表型映射列: sid=", sid_col, " acronym=", acr_col, " type=", typ_col, "\n", sep = "")

sid2acr <- setNames(as.character(pheno[[acr_col]]), as.character(pheno[[sid_col]]))
sid2typ <- setNames(as.character(pheno[[typ_col]]), as.character(pheno[[sid_col]]))

# 生存：OS + OS.time
if (!all(c("sample", "OS", "OS.time") %in% colnames(surv))) {
  stop("生存文件缺 sample/OS/OS.time 列，实际：", paste(colnames(surv), collapse = ", "))
}
surv$OS.time <- suppressWarnings(as.numeric(surv$OS.time))
surv$OS      <- suppressWarnings(as.numeric(surv$OS))
surv <- surv[!is.na(surv$OS.time) & !is.na(surv$OS), ]

# ---------- 4) 对齐：样本 → 癌种，只留原发肿瘤 ----------
samples <- colnames(expr_mat)
acr <- sid2acr[samples]
typ <- sid2typ[samples]
cat("诊断: 表达样本", length(samples),
    "| 癌种匹配", sum(!is.na(acr)),
    "| 类型匹配", sum(!is.na(typ)), "\n")
cat("sample_type 取值示例:", paste(head(unique(typ), 12), collapse = ", "), "\n")
is_primary <- is.na(typ) | grepl("Primary", typ, ignore.case = TRUE)
keep <- !is.na(acr) & acr != "" & is_primary
expr_mat <- expr_mat[, keep, drop = FALSE]
acr <- acr[keep]
cat("原发肿瘤样本:", sum(keep), "| 癌种数:", length(unique(acr)), "\n")

# ---------- 5) 打分（同脚本 20：基因跨样本 z-score → 均值）----------
score_program <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < 3) return(rep(NA_real_, ncol(mat)))
  z <- t(scale(t(mat[g, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}
prog_scores <- sapply(prog_genes, function(g) score_program(expr_mat, g))
rownames(prog_scores) <- colnames(expr_mat)
cov <- sapply(prog_genes, function(g) length(intersect(g, rownames(expr_mat))))
cat("各程序基因覆盖度:", paste(sprintf("%s=%d", names(cov), cov), collapse = " "), "\n")
rm(expr_mat); gc()

# ---------- 6) Cox：程序 × 癌种 ----------
min_samples <- 30
min_events  <- 10
surv_map <- setNames(surv$OS.time, surv$sample)
evt_map  <- setNames(surv$OS,      surv$sample)

rows <- list()
for (ct in sort(unique(acr))) {
  ct_samp <- rownames(prog_scores)[acr == ct]
  t <- surv_map[ct_samp]; e <- evt_map[ct_samp]
  ok <- !is.na(t) & !is.na(e) & t >= 0
  if (sum(ok) < min_samples || sum(e[ok] == 1) < min_events) next
  for (p in colnames(prog_scores)) {
    s <- prog_scores[ct_samp, p]
    z <- as.numeric(scale(s[ok]))
    if (any(is.na(z)) || sd(z) < 1e-6) next
    fit <- tryCatch(coxph(Surv(t[ok], e[ok]) ~ z), error = function(e) NULL)
    if (is.null(fit)) next
    hr <- exp(coef(fit)); ci <- exp(confint(fit)); pv <- summary(fit)$coef[1, 5]
    rows[[paste(ct, p)]] <- data.frame(cancer = ct, program = p,
      n = sum(ok), events = sum(e[ok] == 1),
      HR = hr, lower = ci[1], upper = ci[2], P = pv, stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, rows)
res$padj <- p.adjust(res$P, "BH")
write.csv(res, file.path(tab_dir, "pancancer_program_hr.csv"), row.names = FALSE)

# ---------- 7) P8 汇总 + 森林图 ----------
p8 <- res[res$program == "P8", ]
p8 <- p8[order(p8$HR), ]   # 按 HR 升序：保护在底、不良在顶
p8$sig <- ifelse(p8$padj < 0.05, ifelse(p8$HR > 1, "adverse", "protective"), "ns")
write.csv(p8, file.path(tab_dir, "pancancer_p8_summary.csv"), row.names = FALSE)

cat("\n===== P8 泛癌汇总 =====\n")
cat("HR>1（不良）:", sum(p8$HR > 1), "/", nrow(p8),
    "| HR<1（保护）:", sum(p8$HR < 1), "/", nrow(p8),
    "| 显著不良:", sum(p8$padj < 0.05 & p8$HR > 1),
    "| 显著保护:", sum(p8$padj < 0.05 & p8$HR < 1), "\n\n")
print(p8, row.names = FALSE)

p8$cancer <- factor(p8$cancer, levels = p8$cancer)
cols <- c(adverse = "#C0392B", protective = "#2E86AB", ns = "grey70")
p <- ggplot(p8, aes(x = HR, y = cancer, color = sig)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey40") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.25, linewidth = 0.7) +
  geom_point(size = 2.2) +
  scale_color_manual(values = cols, name = NULL,
                     labels = c(adverse = "P8 不良", protective = "P8 保护", ns = "ns")) +
  labs(x = "Hazard ratio (per SD of P8 score)", y = NULL,
       title = "P8 hypoxia-macrophage program across TCGA cancers (OS)",
       subtitle = "HR>1 = higher P8 → worse OS; colored = BH-adjusted P<0.05") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "pancancer_p8_forest.png"), p, width = 8, height = 9, dpi = 300)

# ---------- 8) 程序 × 癌种 HR 热图（log2 HR，RdBu + 白色分隔线）----------
res$logHR <- log2(res$HR)
# 癌种按 P8 HR 排序（不良→保护），与森林图方向对齐
cancer_lv <- rev(as.character(p8$cancer))
cancer_lv <- cancer_lv[cancer_lv %in% unique(res$cancer)]
res$cancer <- factor(res$cancer, levels = cancer_lv)
m_hr <- max(abs(res$logHR), na.rm = TRUE)
p2 <- ggplot(res, aes(x = cancer, y = program, fill = logHR)) +
  geom_tile(color = "white", linewidth = 0.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-m_hr, m_hr), name = "log2 HR") +
  labs(x = NULL, y = NULL,
       title = "Program prognostic direction across TCGA cancers (OS)",
       subtitle = "red = adverse (HR>1), blue = protective (HR<1)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(fig_dir, "pancancer_program_heatmap.png"), p2, width = 12, height = 5, dpi = 300)

cat("\n已保存 tables/pancancer_program_hr.csv + pancancer_p8_summary.csv\n")
cat("      figures/pancancer_p8_forest.png + pancancer_program_heatmap.png\n")
cat("说明：若 P8 在多数癌种 HR>1（尤其显著），即证明「缺氧巨噬轴」是泛癌不良预后机制。方向如实读，勿 cherry-pick。\n")
