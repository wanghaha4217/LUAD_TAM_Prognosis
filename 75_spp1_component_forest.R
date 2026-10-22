# ============================================================
# 75_spp1_component_forest.R  (v3 — 改用 ggplot2，彻底绕开 forestplot 的空白 PDF 坑)
# 目的：回答「为什么 P1（SPP1+C1q+细胞周期）整体无预后信号」——
#       把 P1 拆到组分，画 SPP1 / C1QA/B/C / CXCL9 在 4 队列的单变量森林图。
#       核心结论可视化：SPP1 4/4 方向不良（GSE31210 显著），
#                      C1q 方向随队列漂移（TCGA/GSE72094 保护、GSE31210 不良）。
# 输入：results/tables/spp1_prognosis.csv（由 74_spp1_prognosis.R 生成）
# 输出：results/figures/spp1_component_forest.png / .pdf
# 依赖：ggplot2（与全项目其余图一致，保证能正常打开）
# 运行：本地 RStudio
# ============================================================

suppressPackageStartupMessages({ library(ggplot2) })

base    <- "."
res_dir <- file.path(base, "results")
fig_dir <- file.path(res_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path(res_dir, "tables", "spp1_prognosis.csv"), stringsAsFactors = FALSE)

FEAT_ORDER <- c("SPP1", "SPP1_minus_CXCL9", "C1QA", "C1QB", "C1QC", "CXCL9")
FEAT_LABEL <- c(
  SPP1              = "SPP1",
  SPP1_minus_CXCL9  = "SPP1 \u2212 CXCL9",
  C1QA              = "C1QA",
  C1QB              = "C1QB",
  C1QC              = "C1QC",
  CXCL9             = "CXCL9"
)
COHORT_ORDER <- c("TCGA-LUAD", "GSE68465", "GSE72094", "GSE31210")

d$feature <- factor(d$feature, levels = FEAT_ORDER, labels = FEAT_LABEL)
d$cohort  <- factor(d$cohort, levels = rev(COHORT_ORDER))   # rev → TCGA 在每格最上
d$sig     <- factor(ifelse(d$P < 0.05, "P < 0.05", "n.s."),
                    levels = c("P < 0.05", "n.s."))
d$star <- ifelse(d$P < 0.001, "***",
                 ifelse(d$P < 0.01, "**",
                        ifelse(d$P < 0.05, "*", "")))
d$lab <- sprintf("%.2f (%.2f\u2013%.2f)%s", d$HR, d$lower, d$upper, d$star)

p <- ggplot(d, aes(x = HR, y = cohort)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3,
                 color = "grey35", linewidth = 0.7) +
  geom_point(aes(color = sig), size = 2.4) +
  geom_text(aes(x = 5.5, label = lab), hjust = 1, size = 3.0) +
  scale_color_manual(values = c("P < 0.05" = "#2E86AB", "n.s." = "grey55"),
                     drop = FALSE) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4, 8), limits = c(0.5, 6.2),
                expand = expansion(mult = c(0, 0.02))) +
  facet_grid(feature ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Hazard ratio (95% CI, log scale)", color = NULL) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.placement  = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 11),
    panel.spacing.y   = unit(0.7, "lines"),
    panel.grid.major.y = element_blank(),
    legend.position    = "bottom",
    axis.title.y       = element_blank(),
    axis.text.y        = element_text(size = 10),
    axis.title.x       = element_text(size = 11),
    plot.margin        = margin(5, 5, 5, 12)
  )

ggsave(file.path(fig_dir, "spp1_component_forest.png"), p,
       width = 8, height = 7, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "spp1_component_forest.pdf"), p,
       width = 8, height = 7)

cat("森林图已保存到:", fig_dir, "\n")
cat("特征组:", length(FEAT_ORDER), "，数据行:", nrow(d), "\n")
