# ============================================================
# 48_metabolism.R
# 目的：代谢重编程（糖酵解/OXPHOS/缺氧/脂肪酸代谢）——验证 P8 缺氧轴的下游代谢表型
#   实测（写后校正）：糖酵解↑、缺氧↑，但 OXPHOS 亦↑（非经典 Warburg，而是「双升」高代谢表型）
# 方法：msigdbr 取 HALLMARK 基因集；打分 = 各基因跨样本 z-score 后取均值
#       （与脚本 38/45 的 CYT/IFN-γ/HLA 同口径，方法一致、可解释）
# 队列：4 bulk（TCGA-LUAD + GSE68465/GSE72094/GSE31210）
# 输出：tables/metabolism_by_risk.csv + figures/metabolism_boxplot.png
# 依赖：msigdbr（CRAN，可选，缺则用内置紧凑基因表）；org.Hs.eg.db / AnnotationDbi / ggplot2
# 运行：本地 RStudio
# ============================================================

suppressPackageStartupMessages({
  library(org.Hs.eg.db); library(AnnotationDbi); library(ggplot2)
})

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

has_msig <- requireNamespace("msigdbr", quietly = TRUE)
if (!has_msig) cat("【缺包】msigdbr -> install.packages('msigdbr')（缺则用内置紧凑基因表）\n")

# ---------- 基因集 ----------
pathways <- c(GLYCOLYSIS = "HALLMARK_GLYCOLYSIS",
              OXPHOS     = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
              HYPOXIA    = "HALLMARK_HYPOXIA",
              FAO        = "HALLMARK_FATTY_ACID_METABOLISM")
if (has_msig) {
  mdb <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  gsets <- lapply(pathways, function(p) unique(mdb$gene_symbol[mdb$gs_name == p]))
} else {
  gsets <- list(
    GLYCOLYSIS = c("HK2","PFKL","PFKP","PKM","LDHA","SLC2A1","GPI","ALDOA","GAPDH","ENO1","PGK1","PGAM1","TPI1","LDHB"),
    OXPHOS = c("NDUFA1","NDUFB1","NDUFS1","SDHA","SDHB","UQCRC1","UQCRC2","CYC1","COX5A","COX6A1","ATP5A1","ATP5B","ATP5C1"),
    HYPOXIA = c("VEGFA","CA9","SLC2A1","LDHA","PGK1","ENO1","BNIP3","BNIP3L","ADM","LOX","PDK1","NDRG1","P4HA1"),
    FAO = c("CPT1A","CPT2","ACADM","ACADVL","HADHA","HADHB","ACAA2","ACSL1","ACADL","ETFDH"))
}
cat("基因集大小:", paste(names(gsets), lengths(gsets), sep = "=", collapse = ", "), "\n")

# ---------- 工具 ----------
collapse_mean <- function(mat, group) {
  mat <- as.matrix(mat); g <- as.character(group)
  k <- !is.na(g); mat <- mat[k, , drop = FALSE]; g <- g[k]
  gf <- factor(g); s <- rowsum(mat, gf, reorder = TRUE)
  n <- as.vector(table(gf)[rownames(s)]); s / n
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
  sym <- trimws(sub("\\s*///.*$", "", sym)); out <- sym; names(out) <- id; out[nzchar(out)]
}
zscore_mean <- function(mat, genes) {
  g <- intersect(genes, rownames(mat)); if (length(g) < 5) return(NULL)
  z <- t(scale(t(mat[g, , drop = FALSE]))); colMeans(z, na.rm = TRUE)
}

# ---------- 4 bulk symbol 矩阵 ----------
tcga <- readRDS(file.path(val_dir, "TCGA-LUAD.rds"))
cnt <- tcga$counts; rn <- rownames(cnt)
ens <- sub("\\..*$", "", rn)
sym <- mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first")
names(sym) <- rn; sym <- sym[!is.na(sym) & sym != ""]
keep <- !is.na(sym[rownames(cnt)])
tcga_sym <- collapse_mean(cnt[keep, , drop = FALSE], sym[rownames(cnt)[keep]])
colnames(tcga_sym) <- substr(colnames(tcga_sym), 1, 15)

array_cohorts <- c(GSE68465 = "GPL96", GSE72094 = "GPL15048", GSE31210 = "GPL570")
geo_sym <- list()
for (g in names(array_cohorts)) {
  es <- readRDS(file.path(val_dir, paste0(g, "_eset.rds")))
  expr <- es$expr
  probe2sym <- read_gpl_soft(file.path(val_dir, paste0(array_cohorts[[g]], ".soft.gz")))
  keep <- !is.na(probe2sym[rownames(expr)])
  geo_sym[[g]] <- collapse_mean(expr[keep, , drop = FALSE], probe2sym[rownames(expr)[keep]])
}
bulk <- c(list("TCGA-LUAD" = tcga_sym), geo_sym)

sig <- readRDS(file.path(res_dir, "rds", "signature_results.rds"))
risk_list <- sig$risk

# ---------- 打分 + 比较 ----------
rows <- list(); long <- list()
for (g in names(bulk)) {
  mat <- bulk[[g]]
  risk <- risk_list[[g]]
  common <- intersect(colnames(mat), names(risk))
  mat <- mat[, common, drop = FALSE]; risk <- risk[common]
  grp <- factor(ifelse(risk > median(risk), "High", "Low"), levels = c("Low", "High"))
  for (p in names(gsets)) {
    v <- zscore_mean(mat, gsets[[p]])
    if (is.null(v)) { cat("  [", g, "]", p, "可匹配基因太少，跳过\n", sep = ""); next }
    w <- wilcox.test(v ~ grp)$p.value
    rows[[paste(g, p)]] <- data.frame(cohort = g, pathway = p,
      Low = round(median(v[grp == "Low"]), 3), High = round(median(v[grp == "High"]), 3),
      P = signif(w, 3), stringsAsFactors = FALSE)
    long[[paste(g, p)]] <- data.frame(cohort = g, pathway = p, value = v, group = grp)
  }
}
res <- do.call(rbind, rows)
res$padj <- p.adjust(res$P, "BH")
cat("\n===== 代谢通路分：High vs Low 风险（4 bulk）=====\n")
print(res, row.names = FALSE)
write.csv(res, file.path(tab_dir, "metabolism_by_risk.csv"), row.names = FALSE)

# ---------- 箱线图 ----------
long_df <- do.call(rbind, long)
long_df$pathway <- factor(long_df$pathway, levels = names(pathways))
long_df$cohort  <- factor(long_df$cohort,  levels = names(bulk))
p <- ggplot(long_df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 0.4, alpha = 0.3, color = "grey30") +
  facet_grid(pathway ~ cohort, scales = "free_y") +
  scale_fill_manual(values = c("Low" = "#2E86AB", "High" = "#C0392B"), name = "Risk") +
  labs(x = NULL, y = "Pathway score (gene z-score mean)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", strip.background = element_rect(fill = "grey92", color = NA),
        panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "metabolism_boxplot.png"), p, width = 11, height = 8, dpi = 300)

cat("\n已保存 tables/metabolism_by_risk.csv + figures/metabolism_boxplot.png\n")
cat("说明：High 组糖酵解↑、缺氧↑（4/4 一致）；但 OXPHOS 非↓而是↑（4/4 方向一致、3/4 显著）——非经典 Warburg，\n")
cat("     而是糖酵解+氧化磷酸化「双升」的高代谢/高增殖表型，与药敏「高风险=增殖=更敏感」的增殖悖论一致。FAO 方向不一致，勿硬编。\n")
