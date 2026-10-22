# ============================================================
# 70_spp1_program_response.R  (v2 — 修 P=0 硬伤 + 零膨胀可视化 + 英文标签)
# SPP1（单基因）+ P1 程序打分 在 non-MPR vs MPR/pCR 的巨噬表达差异
# 目的：方向一收尾 —— SPP1/P1 程序是否在 ICB 耐药组巨噬里升高（呼应 CellChat 通讯结果）
# v2 修改：
#   (1) P 值不再显示 P=0：p < 2.2e-16 一律标 "P < 2.2e-16"（R 机器精度极限）
#   (2) 零膨胀：箱线图 + 半透明抖动散点（jitter），MPR 组"大量 0 + 少数高值"一眼可见
#   (3) 标签英文化：non-MPR / MPR/pCR（去掉下划线）
#   (4) 加图例、theme_classic、去掉 top/right 边框
#   (5) y 轴上限封顶到 99 分位，削掉极端离群点造成的留白
# 运行：本地 RStudio（加载 552MB 巨噬对象，~1 分钟）
# 输出：results/tables/spp1_p1_response.csv + results/figures/spp1_p1_response.png
# ============================================================

suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })

so <- readRDS("./GSE243013_LUAD_TAM_clustered.rds")

# P1 程序 top 50 基因（原始 NMF 版 P1 = SPP1-hypoxic-C1q）
prog <- read.csv("./program_genes_top50.csv", check.names = FALSE)
p1 <- prog$P1; p1 <- p1[!is.na(p1) & p1 != ""]
p1 <- p1[p1 %in% rownames(so)]
cat("P1 程序基因（匹配到对象）:", length(p1), "\n")

# 分组（英文标签，去下划线）
so$response <- factor(ifelse(so$pathological_response == "non-MPR", "non-MPR", "MPR/pCR"),
                      levels = c("non-MPR", "MPR/pCR"))
cat("分组:\n"); print(table(so$response))

# SPP1 单基因 + P1 程序打分
spp1 <- as.numeric(GetAssayData(so, assay = "RNA", layer = "data")["SPP1", ])
so <- AddModuleScore(so, features = list(P1 = p1), name = "P1score")

df <- data.frame(SPP1 = spp1, P1score = so$P1score1, response = so$response)

res <- data.frame(
  gene = c("SPP1", "P1_program"),
  mean_nonMPR = c(mean(df$SPP1[df$response == "non-MPR"]), mean(df$P1score[df$response == "non-MPR"])),
  mean_MPR = c(mean(df$SPP1[df$response == "MPR/pCR"]), mean(df$P1score[df$response == "MPR/pCR"])),
  p_wilcoxon = c(wilcox.test(SPP1 ~ response, df)$p.value,
                 wilcox.test(P1score ~ response, df)$p.value)
)
cat("\n===== SPP1 / P1 程序：non-MPR vs MPR/pCR =====\n")
print(res, digits = 4)
write.csv(res, "./results/tables/spp1_p1_response.csv", row.names = FALSE)

# ---- P 值格式化：杜绝 "P = 0" ----
fmt_p <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < 2.2e-16) return("P < 2.2e-16")
  if (p < 0.001)    return(sprintf("P = %.2e", p))
  sprintf("P = %.3f", p)
}

COL <- c("non-MPR" = "#E64B35", "MPR/pCR" = "#4DBBD5")

# ---- 箱线图 + 半透明抖动散点 ----
box <- function(y, ylab, pval) {
  cap <- quantile(y, 0.995, na.rm = TRUE)         # 99.5 分位封顶，只削最极端 0.5% 离群点，留出呼吸空间
  ggplot(df, aes(x = response, y = y)) +
    geom_boxplot(aes(fill = response), outlier.shape = NA,
                 linewidth = 0.35, width = 0.62, alpha = 0.8) +
    geom_jitter(aes(fill = response), shape = 21, color = "white", stroke = 0.15,
                width = 0.20, size = 1.1, alpha = 0.45) +
    scale_fill_manual(values = COL) +
    coord_cartesian(ylim = c(0, cap)) +
    labs(x = NULL, y = ylab,
         subtitle = paste0("Wilcoxon rank-sum test, ", fmt_p(pval))) +
    theme_classic(base_size = 13) +
    theme(plot.subtitle = element_text(size = 11, color = "grey25"),
          axis.title.y = element_text(size = 12))
}

g1 <- box(df$SPP1, "SPP1 expression", res$p_wilcoxon[1])
g2 <- box(df$P1score, "P1 program score", res$p_wilcoxon[2])

# 合并图例到下方
g <- (g1 + g2) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave("./results/figures/spp1_p1_response.png",
       g, width = 8, height = 5, dpi = 300, bg = "white")
cat("\n已保存 spp1_p1_response.csv / spp1_p1_response.png\n")
