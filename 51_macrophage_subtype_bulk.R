# ============================================================
# 51_macrophage_subtype_bulk.R
# 目的：把 Scissor 找到的巨噬亚型→预后方向（MKI67/VCAN/S100A10/MMP9=坏、
#       FOLR2=好）搬到 bulk 层独立验证。不用反卷积（参考只有巨噬会误用），
#       改用「亚型 marker 打分」：从 scRNA 提各亚型 marker 基因 → 4 bulk 队列打分
#       → High vs Low 风险比较。
# 方法：marker = 各亚型 log2(mean+1) 相对其他亚型的 logFC 前 50；
#       bulk 打分 = marker 基因跨样本 z-score 后取均值（与脚本 38/45/48 同口径）。
# 输入：GSE243013_LUAD_TAM_seurat.rds（巨噬 scRNA，meta$sub_cell_type）
#       4 bulk 符号矩阵（TCGA-LUAD + GSE68465/GSE72094/GSE31210，脚本 48 同款）
# 输出：tables/macrophage_subtype_markers.csv（marker 清单，透明可查）
#       tables/macrophage_subtype_by_risk.csv（亚型打分 High vs Low）
#       figures/macrophage_subtype_heatmap.png + macrophage_subtype_key_boxplot.png
# 依赖：Seurat、org.Hs.eg.db、AnnotationDbi、ggplot2、Matrix
# 运行：本地 RStudio
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
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

# ============================================================
# Part 1：从 scRNA 提各亚型 marker 基因（伪 bulk logFC，非 Seurat FindAllMarkers）
# ============================================================
seu <- readRDS(file.path(base, "GSE243013_LUAD_TAM_seurat.rds"))
cat("scRNA 细胞:", ncol(seu), "| 基因:", nrow(seu), "\n")

cnt     <- GetAssayData(seu, assay = "RNA", layer = "counts")
libsize <- if ("nCount_RNA" %in% colnames(seu@meta.data)) {
  seu@meta.data$nCount_RNA
} else {
  Matrix::colSums(cnt)
}
cpm <- sweep(cnt, 2, libsize + 1, "/") * 1e6   # 稀疏，保持稀疏

sub  <- as.character(seu@meta.data$sub_cell_type)
sf   <- factor(sub)
design <- model.matrix(~ 0 + sf)
colnames(design) <- levels(sf)

counts_per <- cpm %*% design          # genes × subtypes：各亚型 cpm 之和
ns         <- colSums(design)         # 各亚型细胞数
tot_sum    <- rowSums(counts_per)     # 各基因总 cpm 之和
tot_n      <- ncol(seu)

min_cells <- 50
n_top     <- 50
logfc_min <- 0.25
subs <- names(ns)[ns >= min_cells]
cat("参与打分的亚型（细胞>=50）:", length(subs), "\n")

marker_rows <- list(); markers <- list()
for (s in subs) {
  mean_s <- counts_per[, s] / ns[s]
  mean_o <- (tot_sum - counts_per[, s]) / (tot_n - ns[s])
  logfc  <- log2(mean_s + 1) - log2(mean_o + 1)
  ord    <- names(sort(logfc, decreasing = TRUE))
  top    <- head(ord, n_top)
  top    <- top[logfc[top] > logfc_min]
  if (length(top) < 10) top <- head(ord, n_top)   # 兜底：至少 10 个
  markers[[s]] <- top
  marker_rows[[s]] <- data.frame(subtype = s, rank = seq_along(top),
                                 gene = top, logFC = round(logfc[top], 3),
                                 stringsAsFactors = FALSE)
}
marker_df <- do.call(rbind, marker_rows)
write.csv(marker_df, file.path(tab_dir, "macrophage_subtype_markers.csv"), row.names = FALSE)
cat("各亚型 marker 数:\n"); print(vapply(markers, length, integer(1)))

# ============================================================
# Part 2：4 bulk 符号矩阵（与脚本 48 同款）
# ============================================================
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

# ============================================================
# Part 3：打分 + High vs Low 比较
# ============================================================
rows <- list(); long <- list()
for (g in names(bulk)) {
  mat <- bulk[[g]]
  risk <- risk_list[[g]]
  common <- intersect(colnames(mat), names(risk))
  mat <- mat[, common, drop = FALSE]; risk <- risk[common]
  grp <- factor(ifelse(risk > median(risk), "High", "Low"), levels = c("Low", "High"))
  for (s in names(markers)) {
    v <- zscore_mean(mat, markers[[s]])
    if (is.null(v)) next
    w <- wilcox.test(v ~ grp)$p.value
    rows[[paste(g, s)]] <- data.frame(cohort = g, subtype = s,
      Low = round(median(v[grp == "Low"]), 3), High = round(median(v[grp == "High"]), 3),
      diff = round(median(v[grp == "High"]) - median(v[grp == "Low"]), 3),
      P = signif(w, 3), stringsAsFactors = FALSE)
    long[[paste(g, s)]] <- data.frame(cohort = g, subtype = s, value = v, group = grp)
  }
}
res <- do.call(rbind, rows)
res$padj <- p.adjust(res$P, "BH")

# 标注 Scissor 方向（脚本 49：OR>1=坏预后，OR<1=好预后），按后缀匹配
poor_keys <- c("MKI67", "VCAN", "S100A10", "MMP9", "ISG15", "FCGR3A", "CXCL2")
good_keys <- c("FOLR2", "MARCO", "DNAJB1", "CXCL10")
res$scissor <- ifelse(grepl(paste(poor_keys, collapse = "|"), res$subtype), "poor",
               ifelse(grepl(paste(good_keys, collapse = "|"), res$subtype), "good", ""))

cat("\n===== 巨噬亚型打分：High vs Low（diff>0 = High 风险更高）=====\n")
print(res, row.names = FALSE)
write.csv(res, file.path(tab_dir, "macrophage_subtype_by_risk.csv"), row.names = FALSE)

# ============================================================
# Part 4：图
# ============================================================
# 4a 热图：diff（High-Low），星号 = 显著性
res$sig <- ifelse(res$padj < 0.001, "***", ifelse(res$padj < 0.01, "**",
            ifelse(res$padj < 0.05, "*", "")))
p1 <- ggplot(res, aes(x = cohort, y = subtype, fill = diff)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), color = "black", size = 4) +
  scale_fill_gradient2(low = "#2E86AB", mid = "white", high = "#C0392B",
                       midpoint = 0, name = "High-Low\nmedian diff") +
  labs(x = NULL, y = NULL,
       title = "Macrophage subtype signature in high- vs low-risk (4 bulk)",
       subtitle = "diff>0 = subtype signature higher in high-risk; * P<0.05") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())
ggsave(file.path(fig_dir, "macrophage_subtype_heatmap.png"), p1, width = 7, height = 7, dpi = 300)

# 4b 关键亚型箱线图（4 坏 + 1 好）
key_subs <- subs[grepl(paste(c(poor_keys, "FOLR2"), collapse = "|"), subs)]
long_df <- do.call(rbind, long)
long_df <- long_df[long_df$subtype %in% key_subs, ]
long_df$subtype <- factor(long_df$subtype, levels = key_subs)
long_df$cohort  <- factor(long_df$cohort, levels = names(bulk))
p2 <- ggplot(long_df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 0.3, alpha = 0.25, color = "grey30") +
  facet_grid(subtype ~ cohort, scales = "free_y") +
  scale_fill_manual(values = c("Low" = "#2E86AB", "High" = "#C0392B"), name = "Risk") +
  labs(x = NULL, y = "Subtype signature score (gene z-score mean)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "top", strip.background = element_rect(fill = "grey92", color = NA),
        panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "macrophage_subtype_key_boxplot.png"), p2, width = 11, height = 10, dpi = 300)

cat("\n已保存 tables/macrophage_subtype_markers.csv + macrophage_subtype_by_risk.csv\n")
cat("      figures/macrophage_subtype_heatmap.png + macrophage_subtype_key_boxplot.png\n")
cat("说明：若 MKI67/VCAN/S100A10/MMP9 的 diff>0（High 更高）且 FOLR2 的 diff<0（Low 更高），\n")
cat("     即 bulk 层独立验证 Scissor 的「坏亚型=不良、好亚型=良好」结论。方向如实读，勿预设。\n")
