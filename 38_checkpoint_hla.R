# ============================================================
# 38_checkpoint_hla.R
# 目的：签名高风险 vs 低风险组的免疫检查点 + HLA 表达谱（4 bulk 队列）
#   连接「预后签名」→「免疫治疗」：看高风险组是否免疫冷（检查点/HLA 均下调）
#   注：实测方向为「高风险→检查点下调（非上调）」，即高风险=免疫冷表型；
#       与 P8(缺氧轴)主导、IFN/检查点臂在 P6(无预后价值)一致；ICB 阴性见脚本 39。
#   基因映射：TCGA Ensembl→symbol；GEO 探针→symbol（本地 GPL，不联网）
# 依赖：org.Hs.eg.db / AnnotationDbi / survival / ggplot2
# 运行：本地 RStudio（约 2-5 分钟）
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
})

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 基因集 ----------
checkpoint_genes <- c("CD274","PDCD1","PDCD1LG2","CTLA4","LAG3","HAVCR2",
                      "TIGIT","IDO1","CD276","VSIR","LGALS9","CD80","CD86")
hla_genes <- c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","TAPBP","CIITA",
               "HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1")

# ---------- 工具（与 34/35 一致）----------
collapse_mean <- function(mat, group) {
  mat  <- as.matrix(mat); g <- as.character(group)
  keep <- !is.na(g); mat <- mat[keep, , drop = FALSE]; g <- g[keep]
  gf <- factor(g)
  s <- rowsum(mat, gf, reorder = TRUE)
  n <- as.vector(table(gf)[rownames(s)])
  s / n
}
read_gpl_soft <- function(path) {
  con <- gzfile(path, "rt"); lines <- readLines(con, warn = FALSE); close(con)
  tb <- which(lines == "!platform_table_begin")[1]
  te <- which(lines == "!platform_table_end")[1]
  if (is.na(tb) || is.na(te)) stop("找不到平台表标记：", path)
  tab <- read.delim(text = paste(lines[(tb + 1):(te - 1)], collapse = "\n"),
                    header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  id_col  <- grep("^ID$", colnames(tab), ignore.case = TRUE)[1]
  sym_col <- grep("symbol", colnames(tab), ignore.case = TRUE)[1]
  if (is.na(sym_col)) sym_col <- grep("gene", colnames(tab), ignore.case = TRUE)[1]
  id  <- as.character(tab[[id_col]]); sym <- as.character(tab[[sym_col]])
  sym <- trimws(sub("\\s*///.*$", "", sym))
  out <- sym; names(out) <- id; out[nzchar(out)]
}

# ---------- 1) 各队列 symbol 矩阵 ----------
cat("===== 构建 symbol 矩阵 =====\n")

tcga <- readRDS(file.path(val_dir, "TCGA-LUAD.rds"))
cnt <- tcga$counts; rn <- rownames(cnt)
if (grepl("^ENSG", rn[1])) {
  ens <- sub("\\..*$", "", rn)
  sym <- mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first")
  names(sym) <- rn; sym <- sym[!is.na(sym) & sym != ""]
} else sym <- setNames(rn, rn)
keep <- !is.na(sym[rownames(cnt)])
cnt <- cnt[keep, , drop = FALSE]; grp <- sym[rownames(cnt)]
tcga_sym <- collapse_mean(cnt, grp)
colnames(tcga_sym) <- substr(colnames(tcga_sym), 1, 15)
cat("  TCGA symbol 矩阵:", dim(tcga_sym), "\n")

array_cohorts <- c(GSE68465 = "GPL96", GSE72094 = "GPL15048", GSE31210 = "GPL570")
geo_sym <- list()
for (g in names(array_cohorts)) {
  es <- readRDS(file.path(val_dir, paste0(g, "_eset.rds")))
  expr <- es$expr
  probe2sym <- read_gpl_soft(file.path(val_dir, paste0(array_cohorts[[g]], ".soft.gz")))
  keep <- !is.na(probe2sym[rownames(expr)])
  expr <- expr[keep, , drop = FALSE]; sym <- probe2sym[rownames(expr)]
  geo_sym[[g]] <- collapse_mean(expr, sym)
  cat("  ", g, " symbol 矩阵:", dim(geo_sym[[g]]), "\n", sep = "")
}
mat_list <- c(list("TCGA-LUAD" = tcga_sym), geo_sym)

# ---------- 2) 签名风险 + 检查点/HLA 分 ----------
d <- readRDS(file.path(val_dir, "bulk_program_scores.rds"))
sig <- readRDS(file.path(res_dir, "rds", "signature_results.rds"))
risk_list <- sig$risk

gene_score <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < 3) return(NULL)
  z <- t(scale(t(mat[g, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

rows <- list(); per_gene_rows <- list(); long <- list()

for (g in names(mat_list)) {
  mat <- mat_list[[g]]
  risk <- risk_list[[g]]
  common <- intersect(colnames(mat), names(risk))
  mat <- mat[, common, drop = FALSE]
  risk <- risk[common]
  grp <- factor(ifelse(risk > median(risk), "High", "Low"), levels = c("Low", "High"))

  cp <- gene_score(mat, checkpoint_genes)
  hl <- gene_score(mat, hla_genes)
  if (is.null(cp) || is.null(hl)) { cat("  [", g, "] 可匹配基因太少，跳过\n", sep = ""); next }

  cp_p <- wilcox.test(cp ~ grp)$p.value
  hl_p <- wilcox.test(hl ~ grp)$p.value
  rows[[g]] <- data.frame(
    cohort = g,
    cp_med_Low = round(median(cp[grp == "Low"]), 3),
    cp_med_High = round(median(cp[grp == "High"]), 3),
    cp_P = signif(cp_p, 3),
    hla_med_Low = round(median(hl[grp == "Low"]), 3),
    hla_med_High = round(median(hl[grp == "High"]), 3),
    hla_P = signif(hl_p, 3),
    stringsAsFactors = FALSE)

  # 单个检查点基因（供参考）
  for (gg in intersect(checkpoint_genes, rownames(mat))) {
    z <- as.numeric(mat[gg, ])
    p <- wilcox.test(z ~ grp)$p.value
    per_gene_rows[[paste(g, gg)]] <- data.frame(
      cohort = g, gene = gg,
      med_Low = round(median(z[grp == "Low"]), 3),
      med_High = round(median(z[grp == "High"]), 3),
      P = signif(p, 3), stringsAsFactors = FALSE)
  }

  long[[g]] <- rbind(
    data.frame(cohort = g, score_type = "Checkpoint", value = cp, group = grp),
    data.frame(cohort = g, score_type = "HLA",        value = hl, group = grp))
}

res <- do.call(rbind, rows)
per_gene <- do.call(rbind, per_gene_rows)
long_df <- do.call(rbind, long)

cat("\n===== 检查点/HLA 分：High vs Low 风险 =====\n")
print(res, row.names = FALSE)
write.csv(res, file.path(tab_dir, "checkpoint_hla_by_risk.csv"), row.names = FALSE)
per_gene$padj <- p.adjust(per_gene$P, "BH")
write.csv(per_gene, file.path(tab_dir, "checkpoint_genes_by_risk.csv"), row.names = FALSE)

# ---------- 3) 箱线图 ----------
long_df$cohort <- factor(long_df$cohort, levels = names(mat_list))
long_df$score_type <- factor(long_df$score_type, levels = c("Checkpoint", "HLA"))

p <- ggplot(long_df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.4, color = "grey30") +
  facet_grid(score_type ~ cohort, scales = "free_y") +
  scale_fill_manual(values = c("Low" = "#2E86AB", "High" = "#C0392B"), name = "Risk") +
  labs(x = NULL, y = "Score (gene z-score mean)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "grey92", color = NA),
        panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "checkpoint_hla_boxplot.png"), p, width = 10, height = 5.5, dpi = 300)

cat("\n已保存 checkpoint_hla_by_risk.csv / checkpoint_genes_by_risk.csv / checkpoint_hla_boxplot.png\n")
cat("说明：风险分中位二分（High/Low）；P 为 Wilcoxon；检查点/HLA 分为各基因跨样本 z-score 的均值。\n")
cat("     如实读：高风险组 HLA 显著下调（4/4 队列），检查点也下调（方向一致，2/4 显著）→ 高风险=免疫冷。\n")
cat("     即「高风险→检查点上调」的初始假设被数据否掉，勿在正文写反方向。\n")
