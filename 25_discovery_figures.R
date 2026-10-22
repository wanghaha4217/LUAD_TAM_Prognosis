# ============================================================
# 25_discovery_figures.R
# 发现层三图（不依赖表达矩阵，只靠 UMAP + NMF 的 W/H）：
#   Fig1 UMAP：11 个巨噬亚型着色
#   Fig2 NMF 程序热图：每程序 top 基因的 W 权重（块对角结构）
#   Fig3 程序↔亚型 DotPlot：每亚型主导使用哪个程序
# 依赖：Seurat / ggplot2 / pheatmap / RColorBrewer
# 运行：本地 RStudio（先跑 24 探针确认过对象结构）
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

base    <- "."
fig_dir <- file.path(base, "results", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(base, "GSE243013_LUAD_TAM_clustered.rds"))
nmf <- readRDS(file.path(base, "nmf_k10.rds"))

W       <- nmf$W          # 2000 基因 × 10 程序
H       <- nmf$H          # 10 程序 × 38211 细胞
genes   <- nmf$genes      # 2000 个基因名（W 的行）
programs <- nmf$programs  # list：每程序的 top 基因（可选用）

meta    <- obj@meta.data
umap    <- Embeddings(obj, "umap")
subtype <- as.character(meta$sub_cell_type)

# ============================================================
# 诊断（先打印，据此核对 DotPlot 的 H↔细胞对齐）
# ============================================================
cat("===== 诊断 =====\n")
cat("W dim:", paste(dim(W), collapse = "x"),
    " | H dim:", paste(dim(H), collapse = "x"), "\n")
cat("genes 长度:", length(genes), " | head:", paste(head(genes), collapse = ", "), "\n")
cat("W rownames 有?", !is.null(rownames(W)), "| head:",
    paste(head(rownames(W)), collapse = ", "), "\n")
cat("H colnames 有?", !is.null(colnames(H)), "| head:",
    paste(head(colnames(H)), collapse = ", "), "\n")
cat("meta 行名 head:", paste(head(rownames(meta)), collapse = ", "), "\n")
cat("meta cellID head:", paste(head(meta$cellID), collapse = ", "), "\n")
cat("programs 元素数:", length(programs),
    "| 各长度:", paste(lengths(programs), collapse = ","), "\n")
cat("subtype 11 类:", paste(sort(unique(subtype)), collapse = ", "), "\n")

# H 列 ↔ 细胞对齐
if (!is.null(colnames(H))) {
  m1 <- match(colnames(H), rownames(meta))
  m2 <- match(colnames(H), meta$cellID)
  cat("H 列名 == meta 行名 匹配:", sum(!is.na(m1)), "/", ncol(H), "\n")
  cat("H 列名 == cellID   匹配:", sum(!is.na(m2)), "/", ncol(H), "\n")
} else {
  cat("H 无列名 → 假设列顺序 = Cells(obj) 顺序（需你确认）\n")
}
cat("===== 诊断结束 =====\n\n")

# ---- H 对齐验证：用 program_score_per_patient.rds（61×10）反证 ----
# 若 H 列顺序 = Cells(obj) 顺序，按 sampleID 聚合 H 应能复现该表
cat("===== H 对齐验证 =====\n")
psp_file <- file.path(base, "program_score_per_patient.rds")
if (file.exists(psp_file)) {
  psp <- readRDS(psp_file)                      # 61 × 10（患者 × 程序）
  patients <- sort(unique(meta$sampleID))
  pat_H <- matrix(NA_real_, nrow(H), length(patients),
                  dimnames = list(paste0("P", seq_len(nrow(H))), patients))
  for (p in patients) {
    idx <- which(meta$sampleID == p)
    pat_H[, p] <- rowMeans(H[, idx, drop = FALSE])
  }
  common_p <- intersect(rownames(psp), colnames(pat_H))
  cat("患者交集:", length(common_p), "/", nrow(psp), "\n")
  cors <- sapply(seq_len(nrow(H)), function(j) {
    suppressWarnings(cor(as.numeric(psp[common_p, j]),
                         as.numeric(pat_H[j, common_p]), use = "complete.obs"))
  })
  names(cors) <- paste0("P", seq_len(nrow(H)))
  cat("各程序「H聚合」↔「program_score」相关系数:\n")
  print(round(cors, 4))
  cat("（全 ≈1 → H 列顺序 = Cells 顺序，Fig3 对齐正确；≈0/乱 → 顺序错需重排）\n")
} else {
  cat("没有 program_score_per_patient.rds，跳过验证\n")
}
cat("===== 验证结束 =====\n\n")

# ============================================================
# Fig1 UMAP（11 亚型）
# ============================================================
udf <- data.frame(UMAP1 = umap[, 1], UMAP2 = umap[, 2], subtype = subtype)
pal <- setNames(brewer.pal(11, "Set3"), sort(unique(subtype)))
p1 <- ggplot(udf, aes(UMAP1, UMAP2, color = subtype)) +
  geom_point(size = 0.2, alpha = 0.55) +
  scale_color_manual(values = pal) +
  coord_fixed() +
  labs(title = "LUAD TAM subtypes (n = 38,211 cells)") +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
ggsave(file.path(fig_dir, "discovery_1_umap_subtypes.png"), p1,
       width = 9, height = 7, dpi = 300)
cat("Fig1 UMAP 已保存\n")

# ============================================================
# Fig2 NMF 程序热图（W：每程序 top 基因权重，ggplot2 版，绕开 pheatmap 的坑）
# ============================================================
gname <- if (!is.null(rownames(W)) && length(rownames(W)) == nrow(W)) {
  rownames(W)
} else {
  genes
}
topN <- 20
row_idx <- unlist(lapply(seq_len(ncol(W)),
                         function(j) order(W[, j], decreasing = TRUE)[seq_len(topN)]))
subW <- W[row_idx, , drop = FALSE]
lab   <- make.unique(gname[row_idx])            # 基因可重复进多程序，标签加 .1
rownames(subW) <- lab
colnames(subW) <- paste0("P", seq_len(ncol(W)))

df2 <- expand.grid(gene = lab, program = colnames(subW), stringsAsFactors = FALSE)
df2$weight  <- as.vector(subW)
df2$gene    <- factor(df2$gene, levels = lab)
df2$program <- factor(df2$program, levels = colnames(subW))

p2 <- ggplot(df2, aes(program, gene, fill = weight)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("white", "steelblue", "firebrick")) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "NMF k=10 program genes (top 20, weight-scaled)",
       x = NULL, y = NULL, fill = "NMF weight") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 9),
        panel.grid = element_blank())
ggsave(file.path(fig_dir, "discovery_2_nmf_heatmap.png"), p2,
       width = 7, height = 9, dpi = 300)
cat("Fig2 NMF 热图已保存\n")

# ============================================================
# Fig3 程序↔亚型 DotPlot（H：细胞程序使用度）
# ============================================================
# 每个 meta 行（细胞）对应 H 的第几列
if (!is.null(colnames(H))) {
  if (sum(!is.na(match(colnames(H), meta$cellID))) == ncol(H)) {
    h_col <- match(meta$cellID, colnames(H))       # 用 cellID 对齐
  } else {
    h_col <- match(rownames(meta), colnames(H))    # 用行名对齐
  }
} else {
  h_col <- seq_len(ncol(H))                         # 假设顺序一致
}
keep <- !is.na(h_col)
H_use <- H[, h_col[keep], drop = FALSE]
sub_use <- subtype[keep]
cat("Fig3 参与细胞:", sum(keep), "（NA 丢弃", sum(!keep), "）\n")

# 列归一化：每细胞程序使用度总和=1，否则 H 尺度极小、mean usage 被压成 ~0 颜色失效
H_norm <- sweep(H_use, 2, colSums(H_use), "/")

usub <- sort(unique(sub_use))
np   <- nrow(H_norm)
mean_usage <- matrix(NA_real_, np, length(usub),
                     dimnames = list(paste0("P", seq_len(np)), usub))
frac_top   <- matrix(NA_real_, np, length(usub),
                     dimnames = list(paste0("P", seq_len(np)), usub))
for (s in usub) {
  idx  <- which(sub_use == s)
  Hs   <- H_norm[, idx, drop = FALSE]
  mean_usage[, s] <- rowMeans(Hs)
  frac_top[, s]   <- tabulate(apply(Hs, 2, which.max), nbins = np) / length(idx)
}

long <- expand.grid(program = rownames(mean_usage),
                    subtype = colnames(mean_usage),
                    stringsAsFactors = FALSE)
long$mean <- as.vector(mean_usage)
long$frac <- as.vector(frac_top)
long$program <- factor(long$program, levels = rev(paste0("P", seq_len(np))))  # P1 在上、P10 在下

p3 <- ggplot(long, aes(subtype, program)) +
  geom_point(aes(size = frac, color = mean)) +
  scale_size_continuous(range = c(0.5, 9), name = "fraction of cells\n(top program)") +
  scale_color_gradient(low = "white", high = "firebrick", name = "mean\nusage") +
  labs(x = "Macrophage subtype", y = "NMF program") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
ggsave(file.path(fig_dir, "discovery_3_program_subtype_dotplot.png"), p3,
       width = 9, height = 6, dpi = 300)
cat("Fig3 DotPlot 已保存\n")

cat("\n全部完成。输出目录:", fig_dir, "\n")
