# ============================================================
# 61_crossk_figures.R
# 目的：补 Fig S9 —— 跨 k 稳健性 + 独立算法复现 两张 overlap 热图
#   背景：Discussion 第 143 行已引用两组数字（跨 k mean 0.58–0.87；
#         独立算法 union 0.80–1.00），但此前只有 crossk_overlap.py /
#         align_k10.py 的控制台打印，没有可视化。
#   Panel A（跨 k）：以 RcppML k=10 为基准，算每个程序在
#         k=5,6,7,8,9,11 的 best-match overlap（= |A∩B| / min(|A|,|B|)）。
#   Panel B（算法复现）：原始 NMF 包 k=10 vs RcppML k=10 的 overlap 矩阵。
#   数据：./ 下 8 个 top50 基因 CSV。
#   输出：results/figures/crossk_overlap_heatmap.png
#         results/figures/algorithm_overlap_heatmap.png
#         results/tables/crossk_overlap_matrix.csv
#         results/tables/algorithm_overlap_matrix.csv
#         results/tables/algorithm_union_overlap.csv（顺带核对 Discussion 0.80–1.00）
#   依赖：ggplot2（仅此一个；长表手搓，不用 tidyr）
#   运行：本地 RStudio，秒级。
# ============================================================

suppressPackageStartupMessages(library(ggplot2))

input_dir <- "."
out_dir   <- "./results"
fig_dir   <- file.path(out_dir, "figures")
tab_dir   <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 读 CSV：程序名 -> top50 基因向量 ----------
read_prog <- function(path) {
  m   <- read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  out <- lapply(m, function(col) col[!is.na(col) & col != ""])
  # 按 P1, P2, ... 排序
  out[order(as.integer(sub("P", "", names(out))))]
}

# overlap = 交集大小 / 较小集合大小（与两个 py 脚本一致）
ov <- function(a, b) length(intersect(a, b)) / min(length(a), length(b))

# ---------- 读取所有程序集合 ----------
orig <- read_prog(file.path(input_dir, "program_genes_top50.csv"))      # 原始 NMF 包 k=10（论文用）
rcpp <- read_prog(file.path(input_dir, "program_genes_top50_k10.csv"))  # RcppML k=10（独立算法）
other_ks <- c(5, 6, 7, 8, 9, 11)
others <- lapply(other_ks, function(k)
  read_prog(file.path(input_dir, sprintf("program_genes_top50_k%d.csv", k))))
names(others) <- paste0("k=", other_ks)

# ---------- Panel A：跨 k best-match overlap（基准 = RcppML k=10） ----------
mat_A <- matrix(NA_real_, nrow = length(rcpp), ncol = length(others),
                dimnames = list(names(rcpp), names(others)))
for (i in seq_along(rcpp)) {
  for (j in seq_along(others)) {
    mat_A[i, j] <- max(vapply(others[[j]], function(o) ov(rcpp[[i]], o), numeric(1)))
  }
}

# ---------- Panel B：原始 NMF k=10 vs RcppML k=10 overlap 矩阵 ----------
mat_B <- matrix(NA_real_, nrow = length(orig), ncol = length(rcpp),
                dimnames = list(names(orig), names(rcpp)))
for (i in seq_along(orig)) {
  for (j in seq_along(rcpp)) {
    mat_B[i, j] <- ov(orig[[i]], rcpp[[j]])
  }
}

# ---------- union overlap（原始程序 = 最佳 + 次佳 RcppML 的并集）核对 Discussion 0.80–1.00 ----------
union_B <- sapply(seq_along(orig), function(i) {
  o    <- sort(mat_B[i, ], decreasing = TRUE)
  top2 <- names(o)[1:2]
  ov(orig[[i]], unique(unlist(rcpp[top2])))
})
names(union_B) <- names(orig)

# ---------- 生物学标注（与两个 py 脚本一致，用于行标签） ----------
lab_rcpp <- c(
  P1 = "LAM/APOE-TREM2",  P2 = "C1q-complement", P3 = "NR4A-early-response",
  P4 = "Metallothionein/mono", P5 = "SPP1-hypoxic-prolif", P6 = "MHC-II/IFN",
  P7 = "PLIN2-lipid-inflam", P8 = "Heat-shock", P9 = "S100A8/9-mono",
  P10 = "NFkB-chemokine")
lab_orig <- c(
  P1 = "SPP1-hypoxic(+C1q)", P2 = "LAM/APOE-TREM2", P3 = "NR4A-early",
  P4 = "S100A8/9-mono", P5 = "JUN/AP1-inflam", P6 = "MHC-II/IFN",
  P7 = "glycolytic-lipid", P8 = "hypoxia-NFkB-glyco", P9 = "Heat-shock",
  P10 = "NFkB-chemokine")

row_label_A <- paste0(names(rcpp), "  ", lab_rcpp[names(rcpp)])
row_label_B <- paste0(names(orig), "  ", lab_orig[names(orig)])
col_label_A <- names(others)                     # k=5 ... k=11
col_label_B <- paste0("RcppML ", names(rcpp))    # RcppML P1 ... P10

# ---------- 画 overlap 热图 ----------
heat <- function(mat, row_labels, col_labels, xlab, ylab, text_size = 3.0) {
  d <- expand.grid(r = seq_len(nrow(mat)), c = seq_len(ncol(mat)))
  d$value   <- mat[cbind(d$r, d$c)]
  d$row_lab <- factor(row_labels[d$r], levels = row_labels)
  d$col_lab <- factor(col_labels[d$c], levels = col_labels)
  d$label   <- sprintf("%.2f", d$value)
  d$txtcol  <- ifelse(d$value > 0.7, "white", "black")

  ggplot(d, aes(x = col_lab, y = row_lab, fill = value)) +
    geom_tile(color = "grey92", linewidth = 0.3) +
    geom_text(aes(label = label, color = txtcol), size = text_size, show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_gradient(low = "#F7F7F7", high = "#2166AC", limits = c(0, 1),
                        name = "Overlap") +
    scale_y_discrete(limits = rev(levels(d$row_lab))) +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text  = element_text(color = "black", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 11, face = "bold"),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 9)
    )
}

pA <- heat(mat_A, row_label_A, col_label_A, "Best match in", "RcppML k=10 program")
pB <- heat(mat_B, row_label_B, col_label_B, "RcppML program", "Original NMF program",
           text_size = 2.8)

ggsave(file.path(fig_dir, "crossk_overlap_heatmap.png"), pA,
       width = 6.5, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "algorithm_overlap_heatmap.png"), pB,
       width = 7.5, height = 6, dpi = 300)

# ---------- 存档 overlap 矩阵（便于核对 / 投稿附表） ----------
write.csv(mat_A, file.path(tab_dir, "crossk_overlap_matrix.csv"))
write.csv(mat_B, file.path(tab_dir, "algorithm_overlap_matrix.csv"))
write.csv(data.frame(program = names(union_B), union_overlap = unname(union_B),
                     stringsAsFactors = FALSE),
          file.path(tab_dir, "algorithm_union_overlap.csv"), row.names = FALSE)

# ---------- 控制台核对（应与 crossk_overlap.py / align_k10.py 输出一致） ----------
cat("Panel A（跨 k）每程序 mean overlap 范围：",
    paste(round(range(rowMeans(mat_A)), 3), collapse = " ~ "), "\n")
cat("Panel A 各程序 mean overlap：\n")
print(round(sort(rowMeans(mat_A)), 3))
cat("\nPanel B（算法复现）union overlap（应 0.80–1.00）：\n")
print(round(sort(union_B), 3))
cat("\n关键对应（写图注用）：\n")
cat("  原始 P8(hypoxia-NFkB-glyco) -> RcppML",
    names(sort(mat_B["P8", ], decreasing = TRUE))[1], "\n")
cat("  原始 P9(Heat-shock)          -> RcppML",
    names(sort(mat_B["P9", ], decreasing = TRUE))[1], "\n")
cat("\nDone.\n")
