# ============================================================
# 64_regulon_enrichment.R
# 目的：Regulon 富集 —— 谁的靶基因显著富集在 P8？证明 HIF1A/EPAS1/RELA 是 P8 候选 master regulator
# 方法：DoRothEA regulon（A+B 置信）超几何检验（Fisher 精确）
# 内存：极低（<100MB），纯 R，本地 RStudio 跑
# 输入：program_genes_top50.csv + GSE243013 genes.csv + dorothea_hs
# 输出：results/tables/regulon_enrichment_p8.csv + results/figures/regulon_enrichment_p8.png
# ============================================================

suppressPackageStartupMessages({
  library(dorothea)
  library(ggplot2)
})

base    <- "."
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables"); fig_dir <- file.path(res_dir, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 1. 读 P8 top 50 基因 ----------
prog <- read.csv(file.path(base, "program_genes_top50.csv"), check.names = FALSE)
p8 <- prog$P8
p8 <- p8[!is.na(p8) & p8 != ""]
cat("P8 top 基因数:", length(p8), "\n")
cat("P8 基因:", paste(p8, collapse = ", "), "\n\n")

# ---------- 2. 读全基因（背景） ----------
genes <- readLines(gzfile(file.path(base, "GSE243013/GSE243013_genes.csv(1).gz")))
genes <- genes[-1]                                   # 去表头 geneSymbol
N <- length(unique(genes))
cat("背景基因数 N =", N, "\n")

# ---------- 3. DoRothEA regulon（A+B 置信） ----------
data(dorothea_hs, package = "dorothea")
reg <- dorothea_hs[dorothea_hs$confidence %in% c("A", "B"), ]
reg <- reg[reg$target %in% genes, ]                  # 只留背景基因里的靶
cat("DoRothEA A+B regulon 条目:", nrow(reg), "\n\n")

# ---------- 4. 超几何检验：每个 TF 的靶基因在 P8 的富集 ----------
# 候选 TF（P8 上调） + 对照（P8 下调的 IFN 轴）
candidate_tf <- c("HIF1A", "EPAS1", "RELA", "NFKB1", "ARNT", "STAT1", "STAT2", "IRF1", "IRF9")

enrich <- function(tf, mor_val) {
  targets <- unique(reg$target[reg$tf == tf & reg$mor == mor_val])
  K <- length(targets)                       # 背景里该 TF 靶基因数
  k <- sum(p8 %in% targets)                  # P8 里该 TF 靶基因数
  n <- length(p8)                            # P8 基因数
  if (K == 0) return(c(K = 0, k = 0, fold = NA, p = NA))
  p_val <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  fold  <- (k / n) / (K / N)
  c(K = K, k = k, fold = fold, p = p_val)
}

res <- do.call(rbind, lapply(candidate_tf, function(tf) {
  act <- enrich(tf, 1)   # mor=+1 激活靶基因
  rep <- enrich(tf, -1)  # mor=-1 抑制靶基因
  data.frame(
    TF = tf,
    K_act = act["K"], k_act = act["k"], fold_act = act["fold"], p_act = act["p"],
    K_rep = rep["K"], k_rep = rep["k"], fold_rep = rep["fold"], p_rep = rep["p"],
    stringsAsFactors = FALSE
  )
}))

res$padj_act <- p.adjust(res$p_act, method = "BH")
res$padj_rep <- p.adjust(res$p_rep, method = "BH")
res <- res[order(res$p_act), ]
write.csv(res, file.path(tab_dir, "regulon_enrichment_p8.csv"), row.names = FALSE)

cat("===== Regulon 富集结果（mor=+1 激活靶基因）=====\n")
print(res[, c("TF", "K_act", "k_act", "fold_act", "p_act", "padj_act")])
cat("\n===== mor=-1 抑制靶基因 =====")
print(res[, c("TF", "K_rep", "k_rep", "fold_rep", "p_rep", "padj_rep")])

# ---------- 5. 图：激活靶基因富集（-log10 P + 富集倍数） ----------
p_act <- res[!is.na(res$p_act), ]
p_act$TF <- factor(p_act$TF, levels = p_act$TF[order(p_act$p_act)])
p_act$sig <- ifelse(p_act$padj_act < 0.05, "BH<0.05", "ns")
p_act$star <- ifelse(p_act$padj_act < 0.001, "***",
                     ifelse(p_act$padj_act < 0.01, "**",
                            ifelse(p_act$padj_act < 0.05, "*", "")))

gg <- ggplot(p_act, aes(x = TF, y = -log10(p_act), fill = sig)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = sprintf("%s\nk=%d/%d\nfold=%.1f", star, k_act, K_act, fold_act)),
            vjust = -0.2, size = 2.8) +
  scale_fill_manual(values = c("BH<0.05" = "#C0392B", "ns" = "#2E86AB")) +
  labs(x = "Transcription factor", y = "-log10(P)  (activated targets in P8)",
       fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.4),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "top")

ggsave(file.path(fig_dir, "regulon_enrichment_p8.png"), gg, width = 8, height = 5.5, dpi = 300)
cat("\n已保存 regulon_enrichment_p8.csv / .png\n")
