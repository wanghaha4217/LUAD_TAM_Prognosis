# ============================================================
# 31_icb_response.R
# 任务 E：TAM 程序 ↔ ICB 应答（GSE243013 新辅助 病理应答）
#   发现层 61 例 LUAD 的 NMF 程序（program_score_per_patient.rds = 61×10，行名 sampleID）
#   应答分组：Responder = MPR + pCR；Non-responder = non-MPR
#   输出：
#     results/tables/icb_program_response.csv   每程序 Wilcoxon + 序数相关
#     results/figures/icb_response_boxplot.png  Signature + P8/P2/P3 分组箱线图
# 依赖：ggplot2 + Seurat（仅取 meta.data）；AUC 用 Mann-Whitney U 直接算，无需 pROC
# 运行：本地 RStudio
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
})

base    <- "."
res_dir <- file.path(base, "results")
fig_dir <- file.path(res_dir, "figures")
tab_dir <- file.path(res_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

p_lab <- function(p) ifelse(p < 0.001, "P<0.001", sprintf("P=%.3f", p))

# ---- 1) 程序分（61×10）----
psp <- as.data.frame(readRDS(file.path(base, "program_score_per_patient.rds")))
prog_names <- colnames(psp)

# ---- 2) 患者级应答注释（slim meta.data 去重取一行/患者）----
obj  <- readRDS(file.path(base, "GSE243013_LUAD_TAM_clustered.rds"))
meta <- obj@meta.data
keep <- intersect(c("sampleID","pathological_response","anti-PD1_therapy",
                    "chemotherapy","pre_treatment_staging","gender","age","smoking_history"),
                  colnames(meta))
pat  <- meta[!duplicated(meta$sampleID), keep, drop = FALSE]
rownames(pat) <- pat$sampleID

pat$response <- factor(
  ifelse(pat$pathological_response %in% c("MPR","pCR"), "Responder", "Non-responder"),
  levels = c("Non-responder","Responder"))
pat$response_ord <- ifelse(pat$pathological_response == "non-MPR", 0L,
                    ifelse(pat$pathological_response == "MPR", 1L, 2L))

common <- intersect(rownames(psp), rownames(pat))
cat("对齐患者数:", length(common), "/", nrow(psp), "\n")
cat("应答分布:\n"); print(table(pat[common, "response"]))
cat("应答(序数 non-MPR=0/MPR=1/pCR=2)分布:\n"); print(table(pat[common, "response_ord"]))

# ---- 3) 每程序 Wilcoxon（应答 vs 不应答）+ 序数相关（Spearman）----
res <- do.call(rbind, lapply(prog_names, function(pp) {
  x <- psp[common, pp]
  r <- pat[common, "response"]
  wr <- wilcox.test(x ~ r)
  n_r   <- sum(r == "Responder"); n_nr <- sum(r == "Non-responder")
  U     <- as.numeric(wr$statistic)
  auc_nr <- U / (n_r * n_nr)   # = P(不应答者分数 > 应答者)；>0.5 表示"高分→不应答"
  sp <- cor.test(x, pat[common, "response_ord"], method = "spearman")
  data.frame(program = pp,
             med_NonResponder = median(x[r == "Non-responder"]),
             med_Responder    = median(x[r == "Responder"]),
             W = U, P = wr$p.value,
             AUC_high_to_NonResponder = auc_nr,
             rho_ordinal = as.numeric(sp$estimate), P_ordinal = sp$p.value,
             stringsAsFactors = FALSE)
}))
res$padj     <- p.adjust(res$P, "BH")
res$padj_ord <- p.adjust(res$P_ordinal, "BH")
res <- res[order(res$P), ]
rownames(res) <- NULL
cat("\n===== 每程序 ↔ 应答（按 P 升序）=====\n"); print(res)

# ---- 4) 签名分（7 程序 LASSO 系数）----
# 注：NMF H 均值 ≠ bulk 基因表达程序分，属跨尺度映射；先在 61 例内 z-score 再加权，
#     作为队列内相对"程序失衡"分，方向与生存一致（分越高越差）。
lc <- read.csv(file.path(tab_dir, "lasso_coefficients.csv"), stringsAsFactors = FALSE)
lc <- lc[lc$selected == TRUE, c("program","beta")]
cat("\n===== 用于签名的 LASSO 系数 =====\n"); print(lc)

sig <- rep(0, length(common)); names(sig) <- common
for (i in seq_len(nrow(lc))) {
  pp <- lc$program[i]
  if (!pp %in% prog_names) next
  z <- scale(psp[common, pp])[, 1]
  if (anyNA(z)) z[is.na(z)] <- 0
  sig <- sig + lc$beta[i] * z
}
pat2 <- pat[common, ]
sig_wr <- wilcox.test(sig ~ pat2$response)
sig_sp <- cor.test(sig, pat2$response_ord, method = "spearman")
sig_auc <- as.numeric(sig_wr$statistic /
                        (sum(pat2$response == "Responder") * sum(pat2$response == "Non-responder")))
cat("\n===== 签名分 ↔ 应答 =====\n")
cat("  中位数：Non-responder =", round(median(sig[pat2$response=="Non-responder"]),3),
    "| Responder =", round(median(sig[pat2$response=="Responder"]),3), "\n")
cat("  Wilcoxon P =", signif(sig_wr$p.value, 3),
    "| AUC(高分→不应答) =", round(sig_auc, 3),
    "| rho_ord =", round(as.numeric(sig_sp$estimate), 3),
    "| P_ord =", signif(sig_sp$p.value, 3), "\n")

# ---- 5) 敏感性：剔除 anti-PD1 为 No / unknowm 的 3 例 ----
bad <- pat2[["anti-PD1_therapy"]] %in% c("No","unknowm")
cat("\n===== 敏感性（剔除 anti-PD1=No/unknowm 的", sum(bad), "例）=====\n")
if (sum(!bad) >= 6) {
  for (nm in c("Signature", "P8")) {
    v  <- if (nm == "Signature") sig else psp[common, "P8"]
    r2 <- pat2$response[!bad]
    w2 <- wilcox.test(v[!bad] ~ r2)
    cat(" ", nm, ": Wilcoxon P =", signif(w2$p.value, 3),
        "| med(NoRsp/Rsp) =", round(median(v[!bad][r2=="Non-responder"]),3),
        "/", round(median(v[!bad][r2=="Responder"]),3), "\n")
  }
} else {
  cat("  剔除后样本过少，跳过\n")
}

# ---- 6) 箱线图：Signature + P8/P2/P3 ----
measures <- list(
  "Signature"                = sig,
  "P8 (glycolysis-angiogenesis)" = psp[common, "P8"],
  "P2 (lipid/TREM2)"         = psp[common, "P2"],
  "P3 (NR4A)"                = psp[common, "P3"]
)
long <- do.call(rbind, lapply(names(measures), function(nm) {
  data.frame(measure = nm, value = measures[[nm]], response = pat2$response)
}))
long$measure <- factor(long$measure, levels = names(measures))

ann <- data.frame(
  measure = names(measures),
  lab = sapply(measures, function(x) p_lab(wilcox.test(x ~ pat2$response)$p.value)),
  y    = sapply(measures, function(x) max(x) * 1.08)
)
ann$measure <- factor(ann$measure, levels = names(measures))

p <- ggplot(long, aes(x = response, y = value, fill = response)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.3, alpha = 0.55, color = "grey25") +
  geom_text(data = ann, aes(x = 1.5, y = y, label = lab), inherit.aes = FALSE,
            size = 3.2, vjust = 0) +
  facet_wrap(~ measure, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Non-responder" = "#C0392B", "Responder" = "#2E86AB"),
                    name = NULL) +
  labs(x = NULL, y = "Program / signature score (within-cohort)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "grey92", color = NA),
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "icb_response_boxplot.png"), p,
       width = 7.5, height = 6.5, dpi = 300)

# ---- 7) 写表 ----
write.csv(res, file.path(tab_dir, "icb_program_response.csv"), row.names = FALSE)

cat("\n===== 完成 =====\n")
cat("表已保存:", file.path(tab_dir, "icb_program_response.csv"), "\n")
cat("图已保存:", file.path(fig_dir, "icb_response_boxplot.png"), "\n")
