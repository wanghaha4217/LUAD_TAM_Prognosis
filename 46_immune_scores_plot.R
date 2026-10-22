# ============================================================
# 46_immune_scores_plot.R
# 目的：补「免疫评分」箱线图（脚本 45 只存了中位数 CSV、没画图）
#   CYT(Rooney)=GZMA+PRF1；IFN-γ(Ayers)=IDO1/CXCL10/CXCL9/HLA-DRA/STAT1/IFNG；
#   IPS(Charoentong) 用 hacksig（缺包则跳过 IPS，CYT/IFN-γ 照跑）
#   口径与 45 一致（基因跨样本 z-score 后取均值；风险中位二分）；
#   TCGA 计数本脚本先 log2(CPM+1) 归一化，使 CYT/IFN-γ/IPS 跨队列量纲可比
# 输出：figures/immune_scores_boxplot.png + tables/immune_scores_long.csv
# 运行：本地 RStudio（需先跑过 45，依赖 org.Hs.eg.db / AnnotationDbi / ggplot2）
# ============================================================

suppressPackageStartupMessages({
  library(org.Hs.eg.db); library(AnnotationDbi); library(ggplot2)
})

base <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
fig_dir <- file.path(res_dir, "figures"); tab_dir <- file.path(res_dir, "tables")

has_ips <- requireNamespace("hacksig", quietly = TRUE)

collapse_mean <- function(mat, g) {
  mat <- as.matrix(mat); gf <- factor(g); s <- rowsum(mat, gf); s / as.vector(table(gf)[rownames(s)])
}
read_gpl_soft <- function(path) {                      # 与 45 一致：GPL soft 探针→symbol
  con <- gzfile(path, "rt"); L <- readLines(con, warn = FALSE); close(con)
  tb <- which(L == "!platform_table_begin")[1]; te <- which(L == "!platform_table_end")[1]
  tab <- read.delim(text = paste(L[(tb+1):(te-1)], collapse = "\n"),
                    header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  id  <- as.character(tab[[grep("^ID$", colnames(tab), ignore.case = TRUE)[1]]])
  sc  <- grep("symbol", colnames(tab), ignore.case = TRUE)[1]
  sym <- trimws(sub("\\s*///.*$", "", as.character(tab[[sc]])))
  setNames(sym, id)[nzchar(sym)]
}
score <- function(mat, genes) {
  g <- intersect(genes, rownames(mat)); if (length(g) < 2) return(NULL)
  colMeans(t(scale(t(mat[g, , drop = FALSE]))), na.rm = TRUE)
}
CYT_genes  <- c("GZMA", "PRF1")
IFNG_genes <- c("IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG")
ips_of <- function(mat) {
  if (!has_ips) return(NULL)
  tryCatch({
    zmat <- t(scale(t(mat)))                                  # 逐基因 z-score（跨样本）去量纲，跨平台可比
    res  <- hacksig::hack_immunophenoscore(expr_data = zmat)  # 返回 tibble(sample_id, raw_score, ips_score)
    setNames(as.numeric(res$raw_score), as.character(res$sample_id))[colnames(mat)]
  }, error = function(e) { cat("  IPS 失败:", conditionMessage(e), "\n"); NULL })
}

# ---------- symbol 矩阵（4 bulk，与 45 同口径）----------
tcga <- readRDS(file.path(val_dir, "TCGA-LUAD.rds"))
cnt <- tcga$counts; rn <- rownames(cnt)
sym <- mapIds(org.Hs.eg.db, sub("\\..*$", "", rn), "SYMBOL", "ENSEMBL", multiVals = "first")
names(sym) <- rn; sym <- sym[!is.na(sym) & sym != ""]
keep <- !is.na(sym[rn])
cnt_sym <- collapse_mean(cnt[keep, , drop = FALSE], sym[rn[keep]])   # 先按 symbol 合并原始 count
tcga_sym <- log2(sweep(cnt_sym, 2, colSums(cnt_sym), "/") * 1e6 + 1) # 再 log2(CPM+1)，与 GEO log2 表达对齐
colnames(tcga_sym) <- substr(colnames(tcga_sym), 1, 15)

geo_sym <- list()
array_cohorts <- c(GSE68465 = "GPL96", GSE72094 = "GPL15048", GSE31210 = "GPL570")
for (g in names(array_cohorts)) {
  es <- readRDS(file.path(val_dir, paste0(g, "_eset.rds")))
  p2s <- read_gpl_soft(file.path(val_dir, paste0(array_cohorts[[g]], ".soft.gz")))
  keep <- !is.na(p2s[rownames(es$expr)])
  geo_sym[[g]] <- collapse_mean(es$expr[keep, , drop = FALSE], p2s[rownames(es$expr)[keep]])
}
bulk <- c(list("TCGA-LUAD" = tcga_sym), geo_sym)

# ---------- 逐样本评分 → 长表 ----------
sig <- readRDS(file.path(res_dir, "rds", "signature_results.rds"))
long <- list()
for (g in names(bulk)) {
  risk <- sig$risk[[g]]; common <- intersect(colnames(bulk[[g]]), names(risk))
  mat <- bulk[[g]][, common, drop = FALSE]; rk <- risk[common]
  grp <- factor(ifelse(rk > median(rk), "High", "Low"), levels = c("Low", "High"))
  for (tp in c("CYT", "IFNg", "IPS")) {
    v <- switch(tp, CYT = score(mat, CYT_genes), IFNg = score(mat, IFNG_genes), IPS = ips_of(mat))
    if (!is.null(v)) long[[paste(g, tp)]] <- data.frame(cohort = g, score_type = tp, value = as.numeric(v), group = grp)
  }
}
df <- do.call(rbind, long)
df$cohort     <- factor(df$cohort, levels = names(bulk))
df$score_type <- factor(df$score_type, levels = c("CYT", "IFNg", "IPS"))
write.csv(df, file.path(tab_dir, "immune_scores_long.csv"), row.names = FALSE)

# ---------- 箱线图（与 checkpoint_hla_boxplot.png 同风格）----------
p <- ggplot(df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.4, color = "grey30") +
  facet_grid(score_type ~ cohort, scales = "free_y") +
  scale_fill_manual(values = c("Low" = "#2E86AB", "High" = "#C0392B"), name = "Risk") +
  labs(x = NULL, y = "Score") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", strip.background = element_rect(fill = "grey92", color = NA),
        panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "immune_scores_boxplot.png"), p, width = 10, height = 6, dpi = 300)

cat("已保存 figures/immune_scores_boxplot.png + tables/immune_scores_long.csv\n")
cat(sprintf("说明：高风险组 CYT/IFN-γ/IPS 应整体偏低（免疫冷）。本次 IPS：%s。\n",
            if ("IPS" %in% levels(df$score_type)) "已纳入(raw_score)" else "缺包/失败，图中仅 CYT/IFN-γ"))
