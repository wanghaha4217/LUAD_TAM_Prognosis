# ============================================================
# 76_program_correlation.R
# 目的：程序间相关（Pearson）热图 —— 解释 LASSO 为何「空模型/选子集」
#       9 程序（P2-P10）在 TCGA 训练集的打分高度共线 → LASSO 只保留条件贡献
#       最大的程序，共线伙伴被吸收（对应正文 4.2 节「空模型」的解释）
# 输入：validation_cohorts/bulk_program_scores.rds（TCGA-LUAD 的 9 程序打分）
# 输出：results/figures/program_correlation_heatmap.png
# 依赖：ggplot2
# 运行：本地 RStudio（秒级）
# ============================================================

suppressPackageStartupMessages({ library(ggplot2) })

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
fig_dir <- file.path(base, "results", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

d <- readRDS(file.path(val_dir, "bulk_program_scores.rds"))
S <- as.matrix(d$scores[["TCGA-LUAD"]])

progs <- sort(colnames(S))           # P2..P10
S <- S[, progs, drop = FALSE]

C <- cor(S, method = "pearson")
cat("程序 × 程序 Pearson 相关矩阵（TCGA-LUAD）:\n")
print(round(C, 2))

long <- expand.grid(row = rownames(C), col = colnames(C), stringsAsFactors = FALSE)
long$r   <- as.vector(C)
long$row <- factor(long$row, levels = progs)
long$col <- factor(long$col, levels = progs)

gg <- ggplot(long, aes(col, row, fill = r)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3.0, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  coord_fixed() +
  labs(x = NULL, y = NULL,
       title = "Program-program Pearson correlation (TCGA-LUAD, 9 programs)",
       subtitle = "Collinearity explains why LASSO selects a sparse subset") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(fig_dir, "program_correlation_heatmap.png"), gg,
       width = 8, height = 7, dpi = 300, bg = "white")
cat("已保存 program_correlation_heatmap.png\n")
