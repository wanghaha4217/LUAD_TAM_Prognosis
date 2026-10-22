# ============================================================
# 49_scissor.R (v4：补 SNN 图 + NMF 位置对齐)
# ⑤Scissor：找「与不良预后关联」的巨噬（scRNA + TCGA bulk 生存）
#   再验证 Scissor+ 是否富集在 P8 程序 / 特定 Mφ 亚型
# v4 改动（针对前两次报错）：
#   - Scissor 需要 @graphs$RNA_snn（细胞邻居图），子集上跑 FindNeighbors 生成
#   - NMF H 无细胞名 → 位置对齐（H 第 i 列 = Seurat 第 i 个细胞）
#   - 细胞子集 ~8000 + 基因取与 bulk 交集（避免 9.1 GiB 稠密化）
# 输出：tables/scissor_enrichment.csv + figures/scissor_enrichment.png
# 依赖：Scissor、Seurat、org.Hs.eg.db、AnnotationDbi、ggplot2
# ============================================================

suppressPackageStartupMessages({
  library(Scissor)
  library(Seurat)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
})

base <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 0) bulk + 生存 ----------
collapse_mean <- function(mat, group) {
  mat <- as.matrix(mat); g <- as.character(group)
  k <- !is.na(g); mat <- mat[k, , drop = FALSE]; g <- g[k]
  gf <- factor(g); s <- rowsum(mat, gf, reorder = TRUE)
  n <- as.vector(table(gf)[rownames(s)]); s / n
}
tcga <- readRDS(file.path(val_dir, "TCGA-LUAD.rds"))
cnt <- tcga$counts; rn <- rownames(cnt)
ens <- sub("\\..*$", "", rn)
sym <- mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first")
names(sym) <- rn; sym <- sym[!is.na(sym) & sym != ""]
keep <- !is.na(sym[rownames(cnt)])
cpm <- sweep(cnt[keep, , drop = FALSE], 2, colSums(cnt) + 1, "/") * 1e6
bulk <- log2(collapse_mean(cpm, sym[rownames(cnt)[keep]]) + 1)
colnames(bulk) <- substr(colnames(bulk), 1, 15)

d <- readRDS(file.path(val_dir, "bulk_program_scores.rds"))
surv <- d$surv[["TCGA-LUAD"]]
common <- intersect(colnames(bulk), rownames(surv))
bulk <- bulk[, common, drop = FALSE]
phenotype <- data.frame(time = as.numeric(surv[common, "time"]),
                        status = as.integer(surv[common, "event"]),
                        row.names = common)
cat("bulk 样本:", length(common), "| 事件:", sum(phenotype$status), "\n")

# ---------- 1) scRNA + 基因交集 ----------
seu <- readRDS(file.path(base, "GSE243013_LUAD_TAM_seurat.rds"))
cat("scRNA 细胞:", ncol(seu), "| 基因:", nrow(seu), "\n")
genes_common <- intersect(rownames(seu), rownames(bulk))
cat("scRNA ∩ bulk 基因数:", length(genes_common), "\n")

# ---------- 2) 细胞子集（分层抽样 ~8000）----------
set.seed(1)
n_target <- 8000
sub <- as.character(seu@meta.data$sub_cell_type)
cells_keep <- unlist(lapply(unique(sub), function(s) {
  idx <- which(sub == s)
  n <- max(1, ceiling(n_target * length(idx) / ncol(seu)))
  sample(idx, min(n, length(idx)))
}))
cat("子集细胞数:", length(cells_keep), "\n")

# ---------- 3) 建子集对象 + 完整预处理（生成 SNN 图）----------
cnt_sub <- GetAssayData(seu, assay = "RNA", layer = "counts")[genes_common, cells_keep]
meta_sub <- seu@meta.data[cells_keep, , drop = FALSE]

seu_sub <- CreateSeuratObject(counts = cnt_sub, meta.data = meta_sub,
                              min.cells = 0, min.features = 0)
seu_sub <- NormalizeData(seu_sub, verbose = FALSE)
seu_sub <- FindVariableFeatures(seu_sub, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu_sub <- ScaleData(seu_sub, verbose = FALSE)
seu_sub <- RunPCA(seu_sub, npcs = 30, verbose = FALSE)
seu_sub <- FindNeighbors(seu_sub, dims = 1:30, verbose = FALSE)
cat("graphs:", paste(names(seu_sub@graphs), collapse = ", "), "\n")

# 把 RNA assay 转 v4（Scissor 读 @assays$RNA@data），保留 @graphs$RNA_snn
assay_v4 <- CreateAssayObject(counts = GetAssayData(seu_sub, assay = "RNA", layer = "counts"))
assay_v4@data <- GetAssayData(seu_sub, assay = "RNA", layer = "data")
seu_sub@assays[["RNA"]] <- assay_v4
cat("v4 对象: 细胞", ncol(seu_sub), "基因", nrow(seu_sub),
    "| assay 类:", class(seu_sub@assays$RNA),
    "| graphs:", paste(names(seu_sub@graphs), collapse = ", "), "\n")

# ---------- 4) NMF 程序（位置对齐）----------
nmf <- readRDS(file.path(base, "nmf_k10.rds"))
H <- nmf$H
prog_names <- nmf$programs
if (is.list(prog_names) && !is.null(names(prog_names))) prog_names <- names(prog_names)
if (is.null(prog_names) || length(prog_names) != nrow(H)) prog_names <- paste0("P", seq_len(nrow(H)))
prog_of_cell <- prog_names[apply(H, 2, which.max)]
# H 无细胞名 → 位置对齐：H 第 i 列 = Seurat 第 i 个细胞（ncol 相同）
if (ncol(H) == ncol(seu)) {
  names(prog_of_cell) <- colnames(seu)
} else {
  stop("ncol(H) != ncol(seu)，无法位置对齐，先检查 nmf 与 Seurat 的细胞数。")
}
cat("程序名:", paste(prog_names, collapse = ", "), "\n")
cat("H 维度:", paste(dim(H), collapse = " x "), "\n")

# ---------- 5) 跑 Scissor（Cox）----------
set.seed(1)
infos <- Scissor(bulk_dataset = bulk, sc_dataset = seu_sub, phenotype = phenotype,
                 tag = "scissor_tcga", alpha = NULL, cutoff = 0.2, family = "cox")

extract_cells <- function(x) if (is.data.frame(x)) x[[1]] else x
pos_cells <- extract_cells(infos$Scissor_pos)
neg_cells <- extract_cells(infos$Scissor_neg)
cat("Scissor+（不良预后关联）:", length(pos_cells),
    "| Scissor-（良好预后关联）:", length(neg_cells), "\n")

# ---------- 6) 富集 ----------
enrich_grp <- function(labels, pos, neg, min_cells = 5) {
  cells <- intersect(names(labels), c(pos, neg))
  pos <- pos[pos %in% cells]; neg <- neg[neg %in% cells]
  lab <- labels[cells]
  grp <- setNames(ifelse(cells %in% pos, "pos", "neg"), cells)
  do.call(rbind, lapply(sort(unique(lab)), function(p) {
    in_p <- lab == p
    n_tot <- sum(in_p); n_pos <- sum(in_p & grp == "pos"); n_neg <- sum(in_p & grp == "neg")
    if (n_tot < min_cells || n_pos == 0 || n_neg == 0) return(NULL)
    ft <- fisher.test(table(in_p, grp == "pos"))
    data.frame(group = p, n_total = n_tot, n_pos = n_pos, n_neg = n_neg,
               pos_frac = round(n_pos / n_tot, 3),
               OR = round(unname(ft$estimate), 3),
               P = signif(ft$p.value, 3), stringsAsFactors = FALSE)
  }))
}

prog_labels <- prog_of_cell
res_prog <- enrich_grp(prog_labels, pos_cells, neg_cells)
if (is.null(res_prog)) stop("程序富集为空，检查 NMF 细胞与 Scissor 细胞对齐。")
res_prog$var <- "NMF program"

res_sub <- NULL
if ("sub_cell_type" %in% colnames(seu@meta.data)) {
  sub_labels <- setNames(as.character(seu@meta.data$sub_cell_type), colnames(seu))
  res_sub <- enrich_grp(sub_labels, pos_cells, neg_cells)
  if (!is.null(res_sub)) res_sub$var <- "Mφ subtype"
}
res <- res_prog
if (!is.null(res_sub)) res <- rbind(res, res_sub)
overall_pos <- round(length(pos_cells) / (length(pos_cells) + length(neg_cells)), 3)
cat("\n===== Scissor+ 富集（OR>1 = 富集于不良预后巨噬）=====\n")
cat("总体 Scissor+ 比例:", overall_pos, "\n")
print(res, row.names = FALSE)
write.csv(res, file.path(tab_dir, "scissor_enrichment.csv"), row.names = FALSE)

# ---------- 7) 图 ----------
plot_df <- res_prog
plot_df$group <- factor(plot_df$group, levels = intersect(prog_names, unique(plot_df$group)))
plot_df$sig <- ifelse(plot_df$P < 0.05, "P<0.05", "ns")

# Wilson 95% CI（二项比例 n_pos/n_total 的置信区间）
wilson_ci <- function(k, n) {
  z <- qnorm(0.975); p <- k / n; d <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / d
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(lo = center - half, hi = center + half)
}
ciw <- t(mapply(wilson_ci, plot_df$n_pos, plot_df$n_total))
plot_df$lo <- ciw[, "lo"]; plot_df$hi <- ciw[, "hi"]

p <- ggplot(plot_df, aes(x = group, y = pos_frac, fill = sig)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, linewidth = 0.4, color = "grey25") +
  geom_hline(yintercept = overall_pos, linetype = 2, color = "grey40") +
  geom_text(aes(y = hi, label = paste0("n=", n_total, "\nP=", P)),
            vjust = 0, size = 2.8, lineheight = 0.9) +
  scale_fill_manual(values = c("P<0.05" = "#C0392B", "ns" = "grey70"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = "NMF program", y = "Fraction Scissor+ (poor-prognosis)",
       title = "Program enrichment in poor-prognosis macrophages (Scissor)",
       subtitle = paste0("dashed line = overall Scissor+ fraction ", overall_pos)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 9, color = "grey30"))
ggsave(file.path(fig_dir, "scissor_enrichment.png"), p, width = 8, height = 5, dpi = 300)

cat("\n已保存 tables/scissor_enrichment.csv + figures/scissor_enrichment.png\n")
cat("说明：若 P8(Program_8) 的 OR 最大/pos_frac 最高，即独立验证『P8 缺氧程序巨噬 = 不良预后巨噬』。\n")
