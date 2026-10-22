# ============================================================
# 50_viper.R
# ⑦ SCENIC→VIPER（SCENIC-lite）：TF 活性分析，找 P8（缺氧）程序的主控转录因子
#   口径说明：真 SCENIC 需 Python + cisTarget 基序库（从头推断 regulon，重）。
#            本脚本用 DoRothEA 预注释 regulon（TF→靶基因，A-E 置信度）+ VIPER
#            算「每个巨噬的 TF 活性」，等价于 SCENIC 的活性输出、但省去从头建模。
#   问题：哪些 TF 驱动 P8 程序？（已由脚本 49 Scissor 证明 P8=最强不良预后程序，
#         故 P8 主控 TF ≈ 不良预后巨噬的主控 TF）
#   预期（若缺氧+NF-κB 故事成立）：HIF1A/ARNT/NFKB1/RELA 应排前列。
# 输入：GSE243013_LUAD_TAM_seurat.rds（巨噬 scRNA）+ nmf_k10.rds（H 矩阵）
# 输出：tables/viper_p8_tf_activity.csv + figures/viper_p8_tf_barplot.png
# 依赖：viper + dorothea（BiocManager::install(c("viper","dorothea"))）；Seurat
# 运行：本地 RStudio
# ============================================================

suppressPackageStartupMessages({
  library(viper)
  library(dorothea)
  library(Seurat)
  library(ggplot2)
})

base    <- "."
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 0) 取 DoRothEA regulon（跨版本兼容，多方法兜底）----------
get_regulon_df <- function() {
  norm <- function(x) {
    if (!is.data.frame(x)) return(NULL)
    tf <- if ("tf" %in% colnames(x)) as.character(x$tf)
          else if ("source" %in% colnames(x)) as.character(x$source) else NULL
    tg <- if ("target" %in% colnames(x)) as.character(x$target) else NULL
    if (is.null(tf) || is.null(tg)) return(NULL)
    mor <- if ("mor" %in% colnames(x)) as.numeric(x$mor) else rep(1, nrow(x))
    if (!"mor" %in% colnames(x) && all(c("is_stimulation", "is_inhibition") %in% colnames(x))) {
      mor <- as.numeric(x$is_stimulation) - as.numeric(x$is_inhibition)
    }
    conf <- if ("confidence" %in% colnames(x)) as.character(x$confidence)
            else if ("dorothea_level" %in% colnames(x)) as.character(x$dorothea_level)
            else rep("B", nrow(x))
    lik <- if ("likelihood" %in% colnames(x)) as.numeric(x$likelihood)
           else c(A = 0.9, B = 0.7, C = 0.5, D = 0.3, E = 0.1)[toupper(conf)]
    d <- data.frame(tf = tf, target = tg, mor = mor, confidence = conf,
                    likelihood = as.numeric(lik), stringsAsFactors = FALSE)
    d[!is.na(d$tf) & !is.na(d$target) & d$tf != "" & d$target != "", ]
  }
  saved <- file.path(tab_dir, "dorothea_regulons.rds")
  if (file.exists(saved)) { d <- readRDS(saved); if (is.data.frame(d)) return(d) }

  # 1) 老版 dorothea_hs（lazy data，不在 namespace，必须用 data() 加载）
  z <- tryCatch({
    data("dorothea_hs", package = "dorothea", envir = environment())
    if (exists("dorothea_hs", inherits = FALSE)) dorothea_hs else NULL
  }, error = function(e) NULL)
  d <- norm(z); if (!is.null(d)) return(d)

  # 2) get_dorothea 函数（2.x）
  z <- tryCatch({
    f <- get("get_dorothea", envir = asNamespace("dorothea"))
    f(organism = "human", levels = c("A", "B", "C"))
  }, error = function(e) NULL)
  d <- norm(z); if (!is.null(d)) return(d)

  # 3) Collection 对象（2.x）
  z <- tryCatch({
    o <- get("dorothea", envir = asNamespace("dorothea"))
    if (is.function(o)) NULL else as.data.frame(o)
  }, error = function(e) NULL)
  d <- norm(z); if (!is.null(d)) return(d)

  # 4) Omnipath 网络下载（版本无关，需联网）
  z <- tryCatch(
    read.delim("https://omnipathdb.org/interactions?datasets=dorothea&dorothea_levels=A,B",
               stringsAsFactors = FALSE), error = function(e) NULL)
  d <- norm(z); if (!is.null(d)) return(d)

  stop("所有方法都无法获取 DoRothEA regulon。")
}
df <- get_regulon_df()
cat("原始 regulon 行数:", nrow(df), "| 列:", paste(colnames(df), collapse = ", "), "\n")

# 只留高置信 A/B（A=实验验证最强，B=次强）；需要时可用 A+B+C
if ("confidence" %in% colnames(df)) {
  df <- df[df$confidence %in% c("A", "B"), ]
}
# mor（调控方向：+1 激活 / -1 抑制），缺失则默认为激活
if (!"mor" %in% colnames(df)) df$mor <- 1
# likelihood（VIPER 权重），缺失则按 confidence 映射
if (!"likelihood" %in% colnames(df)) {
  conf_val <- c(A = 0.9, B = 0.7, C = 0.5, D = 0.3, E = 0.1)
  df$likelihood <- conf_val[as.character(df$confidence)]
}
df <- df[!is.na(df$mor) & !is.na(df$likelihood), ]
cat("A+B regulon：TF 数", length(unique(df$tf)), "| 靶基因", length(unique(df$target)), "\n")

# 转成 viper 需要的 regulon 列表
regulon <- lapply(split(df, df$tf), function(x) {
  list(tfmode     = setNames(as.numeric(x$mor), x$target),
       likelihood = setNames(as.numeric(x$likelihood), x$target))
})
regulon <- regulon[sapply(regulon, function(r) length(r$tfmode) >= 10)]  # 至少 10 个靶
cat("过滤后 TF（靶>=10）:", length(regulon), "\n")

# ---------- 1) scRNA + NMF P8 ----------
seu <- readRDS(file.path(base, "GSE243013_LUAD_TAM_seurat.rds"))
cat("scRNA 细胞:", ncol(seu), "| 基因:", nrow(seu), "\n")

nmf <- readRDS(file.path(base, "nmf_k10.rds"))
H <- nmf$H
prog_names <- nmf$programs
if (is.list(prog_names) && !is.null(names(prog_names))) prog_names <- names(prog_names)
if (is.null(prog_names) || length(prog_names) != nrow(H)) prog_names <- paste0("P", seq_len(nrow(H)))
cat("程序名:", paste(prog_names, collapse = ", "), "\n")

p8_idx <- which(prog_names == "Program_8")
if (length(p8_idx) != 1) p8_idx <- grep("8", prog_names)[1]
if (length(p8_idx) != 1) stop("找不到 Program_8，程序名:", paste(prog_names, collapse = ", "))
p8_name <- prog_names[p8_idx]

# H 无细胞名 → 位置对齐（同脚本 49）
if (ncol(H) == ncol(seu)) {
  prog_of_cell <- prog_names[apply(H, 2, which.max)]
  names(prog_of_cell) <- colnames(seu)
  p8_score <- H[p8_idx, ]; names(p8_score) <- colnames(seu)
} else {
  stop("ncol(H) != ncol(seu)，无法位置对齐。")
}
cat("P8 程序名:", p8_name, "| P8 主导细胞数:", sum(prog_of_cell == p8_name), "\n")

# ---------- 2) 细胞子集：保留全部 P8 + 抽样非 P8（控内存，8000）----------
set.seed(1)
n_total <- 8000
p8_all    <- which(prog_of_cell == p8_name)
nonp8_all <- which(prog_of_cell != p8_name)
nonp8_n   <- max(1, n_total - length(p8_all))
nonp8_keep <- sample(nonp8_all, min(nonp8_n, length(nonp8_all)))
cells_keep <- c(p8_all, nonp8_keep)
is_p8_cells <- setNames(prog_of_cell[cells_keep] == p8_name, names(prog_of_cell[cells_keep]))
cat("子集细胞:", length(cells_keep), "| P8:", sum(is_p8_cells), "| 非P8:", sum(!is_p8_cells), "\n")

# ---------- 3) 表达矩阵：regulon 基因（TF+靶）→ log2(CPM+1) ----------
target_genes <- unique(unlist(lapply(regulon, function(r) names(r$tfmode))))
tf_genes     <- names(regulon)
use_genes    <- intersect(unique(c(target_genes, tf_genes)), rownames(seu))
cat("regulon 基因 ∩ scRNA:", length(use_genes), "\n")

cnt_sub <- as.matrix(GetAssayData(seu, assay = "RNA", layer = "counts")[use_genes, cells_keep])
libsize <- if ("nCount_RNA" %in% colnames(seu@meta.data)) {
  seu@meta.data$nCount_RNA[cells_keep]
} else {
  colSums(GetAssayData(seu, assay = "RNA", layer = "counts")[, cells_keep])
}
expr <- log2(sweep(cnt_sub, 2, libsize + 1, "/") * 1e6 + 1)

# ---------- 4) VIPER 算 TF 活性 ----------
set.seed(1)
vpres <- viper::viper(expr, regulon, method = "scale", minsize = 10, verbose = FALSE)
cat("VIPER 输出 TF 数:", nrow(vpres), "| 细胞:", ncol(vpres), "\n")

# ---------- 5) 差异：P8 vs 非P8 ----------
tfs <- rownames(vpres)
res <- do.call(rbind, lapply(tfs, function(tf) {
  a <- vpres[tf, is_p8_cells]; b <- vpres[tf, !is_p8_cells]
  w <- wilcox.test(a, b)$p.value
  sp <- suppressWarnings(cor.test(as.numeric(vpres[tf, ]), p8_score[cells_keep],
                                  method = "spearman"))
  data.frame(TF = tf, P8_mean = mean(a), nonP8_mean = mean(b),
             diff = mean(a) - mean(b), Wilcoxon_P = w,
             Spearman_rho = unname(sp$estimate), Spearman_P = unname(sp$p.value),
             stringsAsFactors = FALSE)
}))
res$padj <- p.adjust(res$Wilcoxon_P, "BH")
res <- res[order(res$diff, decreasing = TRUE), ]
cat("\n===== P8 主控 TF（diff = P8 活性 - 非P8 活性）=====\n")
cat("前 10（P8 上调 TF）:\n");  print(head(res, 10), row.names = FALSE)
cat("\n后 10（P8 下调 TF）:\n"); print(tail(res, 10), row.names = FALSE)
write.csv(res, file.path(tab_dir, "viper_p8_tf_activity.csv"), row.names = FALSE)

# 关键 TF 检查（缺氧 + NF-κB 轴）
key_tfs <- c("HIF1A","ARNT","EPAS1","HIF3A","NFKB1","RELA","RELB","NFKB2","STAT3","JUN","FOS")
cat("\n===== 缺氧/NF-κB 相关 TF 在 P8 中的位置（按 diff 排名）=====\n")
kk <- res[res$TF %in% key_tfs, ]
if (nrow(kk)) { print(kk[, c("TF","diff","P8_mean","nonP8_mean","padj")], row.names = FALSE) } else cat("（这些 TF 未在 regulon 中或全被过滤）\n")

# ---------- 6) 图：前 20 上调 + 前 20 下调 ----------
n_show <- 20
top <- head(res, n_show); bot <- tail(res, n_show)
plot_df <- rbind(
  data.frame(top, direction = "up"),
  data.frame(bot, direction = "down")
)
plot_df <- plot_df[order(plot_df$diff), ]
plot_df$TF <- factor(plot_df$TF, levels = plot_df$TF)
plot_df$direction <- factor(plot_df$direction, levels = c("up", "down"))
plot_df$nlog <- -log10(plot_df$padj)

p <- ggplot(plot_df, aes(x = TF, y = diff, fill = direction)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("up" = "#C0392B", "down" = "#2E86AB"),
                    labels = c("up" = "P8 上调", "down" = "P8 下调"), name = NULL) +
  labs(x = NULL, y = "TF activity difference (P8 - non-P8)",
       title = paste0("Master regulators of ", p8_name, " (hypoxia) program"),
       subtitle = "DoRothEA A+B regulons + VIPER; top/bottom 20 by activity difference") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "viper_p8_tf_barplot.png"), p, width = 7, height = 8, dpi = 300)

cat("\n已保存 tables/viper_p8_tf_activity.csv + figures/viper_p8_tf_barplot.png\n")
cat("说明：若 HIF1A/ARNT/NFKB1/RELA 排前列，即独立佐证 P8=缺氧+NF-κB 程序、且其活性由这些 TF 驱动。\n")
cat("     这是 VIPER 活性（regulon 富集），非表达量；方向如实读，勿预设。\n")
