# ============================================================
# 39_icb_external.R
# 目的：TAM 程序预后签名 → 免疫治疗(ICB) 外部验证（两个独立队列）
#   GSE135222 = 27 例晚期 NSCLC anti-PD-1/PD-L1（RNA-seq，Ensembl 基因ID，线性表达值），
#              结局 = PFS（"progression-free survival (pfs)": 1=进展事件，pfs.time=天）→ Cox + KM
#   GSE91061  = 黑色素瘤 anti-PD-1（Riaz 2017，FPKM，Entrez 基因ID），
#              只用 pre-treatment 样本，结局 = response(PRCR vs PD/SD，剔除 UNK)→ Wilcoxon
#   签名分 = Σ β_prog ×（队列内 z-score 后的程序分），与脚本 31 同口径（队列内相对失衡）
#   程序分 = 各程序 top50 基因表达 log2(x+1) 后跨样本 z-score 取均值
# 依赖：org.Hs.eg.db / AnnotationDbi / survival / ggplot2
# 运行：本地 RStudio（约 1-2 分钟）
# ============================================================

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(survival)
  library(ggplot2)
})

base    <- "."
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 工具 ----------
collapse_mean <- function(mat, group) {
  mat  <- as.matrix(mat); g <- as.character(group)
  keep <- !is.na(g); mat <- mat[keep, , drop = FALSE]; g <- g[keep]
  gf <- factor(g)
  s <- rowsum(mat, gf, reorder = TRUE)
  n <- as.vector(table(gf)[rownames(s)])
  s / n
}
# 读本地 series_matrix.txt.gz → 每样本一行（title/geo_accession/各 characteristics 拆列）
read_series_matrix <- function(path) {
  con <- gzfile(path, "rt"); L <- readLines(con, warn = FALSE); close(con)
  L <- L[grepl("^!Sample_", L)]
  parts <- strsplit(L, "\t", fixed = TRUE)
  tag <- vapply(parts, function(p) p[1], character(1))
  ga <- parts[[which(tag == "!Sample_geo_accession")[1]]][-1]
  out <- data.frame(geo_accession = gsub('"', '', ga), stringsAsFactors = FALSE)
  ti <- which(tag == "!Sample_title")[1]
  if (length(ti)) out$title <- gsub('"', '', parts[[ti]][-1])
  for (i in which(tag == "!Sample_characteristics_ch1")) {
    v <- parts[[i]][-1]
    v <- gsub('^"|"$', '', v)                 # 去首尾引号
    nm <- sub(":.*$", "", v[1])               # 特性名（第一个冒号前）
    out[[nm]] <- sub("^[^:]*:\\s*", "", v)    # 值（去 "名: " 前缀）
  }
  out
}
# 程序打分：基因先 log2(x+1)，再跨样本 z-score，取均值
score_program <- function(mat_sym, genes) {
  g <- intersect(genes, rownames(mat_sym))
  if (length(g) < 3) return(rep(NA_real_, ncol(mat_sym)))
  z <- t(scale(t(mat_sym[g, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}
# 签名分 = Σ β × z(程序分)
signature_of <- function(score_mat, coefs) {
  s <- rep(0, nrow(score_mat)); names(s) <- rownames(score_mat)
  for (i in seq_len(nrow(coefs))) {
    pp <- coefs$program[i]
    if (!pp %in% colnames(score_mat)) next
    z <- scale(score_mat[, pp])[, 1]
    z[is.na(z)] <- 0
    s <- s + coefs$beta[i] * z
  }
  s
}

# ---------- 程序基因 + LASSO 系数 ----------
prog_df <- read.csv(file.path(base, "program_genes_top50.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
lc <- read.csv(file.path(tab_dir, "lasso_coefficients.csv"), stringsAsFactors = FALSE)
lc <- lc[lc$selected == TRUE, c("program", "beta")]
cat("签名系数（7 程序）:\n"); print(lc)

# ============================================================
# 1. GSE135222（NSCLC，PFS）
# ============================================================
cat("\n===== GSE135222（NSCLC anti-PD-1/PD-L1，PFS）=====\n")
dir1 <- file.path(base, "GSE135222")
exp1 <- read.delim(gzfile(file.path(dir1, "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz")),
                   row.names = 1, check.names = FALSE)
ens1 <- sub("\\..*$", "", rownames(exp1))
sym1 <- mapIds(org.Hs.eg.db, ens1, "SYMBOL", "ENSEMBL", multiVals = "first")
keep1 <- !is.na(sym1)
mat1 <- log2(as.matrix(exp1[keep1, , drop = FALSE]) + 1)
mat1 <- collapse_mean(mat1, sym1[keep1])
cat("  symbol 矩阵:", dim(mat1), "\n")

ps1 <- sapply(colnames(prog_df), function(p) score_program(mat1, prog_df[[p]]))
rownames(ps1) <- colnames(mat1)
sig1 <- signature_of(ps1, lc)

clin1 <- read_series_matrix(file.path(dir1, "GSE135222_series_matrix.txt.gz"))
clin1$sample_key <- gsub(" ", "", clin1$title)   # "NSCLC 990" → "NSCLC990" 对齐表达列
clin1$pfs       <- as.numeric(clin1[["progression-free survival (pfs)"]])
clin1$pfs_time  <- as.numeric(clin1[["pfs.time"]])

common1 <- intersect(colnames(mat1), clin1$sample_key)
idx1 <- match(common1, clin1$sample_key)
cl1 <- clin1[idx1, ]; rownames(cl1) <- cl1$sample_key
cl1 <- cl1[common1, , drop = FALSE]
sig1 <- sig1[common1]

# pfs=1 表示进展（事件），pfs=0 表示未进展/截尾（GEO 惯例，本数据 21 个 1 / 6 个 0 与此一致）
ok1 <- !is.na(sig1) & !is.na(cl1$pfs_time) & !is.na(cl1$pfs) & cl1$pfs_time > 0
cat("  可分析样本:", sum(ok1), " | PFS 事件:", sum(cl1$pfs[ok1] == 1), "\n")
t1 <- cl1$pfs_time[ok1]; e1 <- cl1$pfs[ok1]; s1 <- sig1[ok1]

# Cox（连续风险分 → PFS）
fit1 <- coxph(Surv(t1, e1) ~ s1)
hr1 <- exp(coef(fit1)); p1 <- summary(fit1)$coefficients[1, "Pr(>|z|)"]
c1 <- as.numeric(summary(fit1)$concordance[1])
cat("  Cox(连续): HR =", round(hr1, 3), "| P =", signif(p1, 3), "| C =", round(c1, 3), "\n")

# 中位二分 log-rank
grp1 <- factor(ifelse(s1 > median(s1), "High", "Low"), levels = c("Low", "High"))
lr1 <- survdiff(Surv(t1, e1) ~ grp1)$pvalue
hr1b <- exp(coef(coxph(Surv(t1, e1) ~ grp1)))
cat("  二分 log-rank: P =", signif(lr1, 3), "| HR(High/Low) =", round(hr1b, 3), "\n")

tab1 <- data.frame(cohort = "GSE135222", cancer = "NSCLC", outcome = "PFS",
                   N = sum(ok1), events = sum(e1),
                   HR_cont = round(hr1, 3), P_cont = signif(p1, 3), Cindex = round(c1, 3),
                   HR_highVSlow = round(hr1b, 3), logrank_P = signif(lr1, 3),
                   stringsAsFactors = FALSE)
write.csv(tab1, file.path(tab_dir, "icb_external_gse135222_pfs.csv"), row.names = FALSE)

# KM 图（base R）
fit_km <- survfit(Surv(t1, e1) ~ grp1)
png(file.path(fig_dir, "icb_pfs_km_gse135222.png"), width = 1400, height = 1100, res = 200)
plot(fit_km, col = c("#2E86AB", "#C0392B"), lwd = 2.5, mark.time = TRUE,
     xlab = "PFS (days)", ylab = "Progression-free survival probability",
     main = "GSE135222 (NSCLC, anti-PD-1/PD-L1)")
legend("topright", legend = c("Low risk", "High risk"),
       col = c("#2E86AB", "#C0392B"), lwd = 2.5)
dev.off()

# ============================================================
# 2. GSE91061（黑色素瘤，response，Pre 样本）
# ============================================================
cat("\n===== GSE91061（黑色素瘤 anti-PD-1，pre-treatment）=====\n")
dir2 <- file.path(base, "GSE91061")
exp2 <- read.csv(gzfile(file.path(dir2, "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz")),
                 row.names = 1, check.names = FALSE)
ent2 <- rownames(exp2)
sym2 <- mapIds(org.Hs.eg.db, ent2, "SYMBOL", "ENTREZID", multiVals = "first")
keep2 <- !is.na(sym2)
mat2 <- log2(as.matrix(exp2[keep2, , drop = FALSE]) + 1)
mat2 <- collapse_mean(mat2, sym2[keep2])
cat("  symbol 矩阵:", dim(mat2), "\n")

ps2 <- sapply(colnames(prog_df), function(p) score_program(mat2, prog_df[[p]]))
rownames(ps2) <- colnames(mat2)
sig2 <- signature_of(ps2, lc)

clin2 <- read_series_matrix(file.path(dir2, "GSE91061_series_matrix.txt.gz"))
clin2$response <- clin2[["response"]]
clin2$visit    <- clin2[["visit (pre or on treatment)"]]

common2 <- intersect(colnames(mat2), clin2$title)
idx2 <- match(common2, clin2$title)
cl2 <- clin2[idx2, ]; rownames(cl2) <- cl2$title
cl2 <- cl2[common2, , drop = FALSE]
s2 <- sig2[common2]

# 只保留 Pre + 有效应答（剔除 UNK/NA）
use2 <- cl2$visit == "Pre" & cl2$response %in% c("PRCR", "PD", "SD") & !is.na(s2)
cat("  Pre 样本:", sum(cl2$visit == "Pre"),
    " | Pre+有效应答:", sum(use2), "\n")
s2 <- s2[use2]
resp2 <- factor(ifelse(cl2$response[use2] == "PRCR", "Responder", "Non-responder"),
                levels = c("Non-responder", "Responder"))
cat("  应答分布:\n"); print(table(resp2))

wr2 <- wilcox.test(s2 ~ resp2)
n_r <- sum(resp2 == "Responder"); n_nr <- sum(resp2 == "Non-responder")
auc2 <- as.numeric(wr2$statistic) / (n_r * n_nr)   # >0.5 = 高分→不应答
cat("  Wilcoxon: P =", signif(wr2$p.value, 3),
    "| med(NonRsp/Rsp) =", round(median(s2[resp2 == "Non-responder"]), 3),
    "/", round(median(s2[resp2 == "Responder"]), 3),
    "| AUC(高分→不应答) =", round(auc2, 3), "\n")

tab2 <- data.frame(cohort = "GSE91061", cancer = "Melanoma", outcome = "response",
                   N = length(resp2), responders = sum(resp2 == "Responder"),
                   med_NonResponder = round(median(s2[resp2 == "Non-responder"]), 3),
                   med_Responder = round(median(s2[resp2 == "Responder"]), 3),
                   AUC_high_to_NonResponder = round(auc2, 3),
                   Wilcoxon_P = signif(wr2$p.value, 3),
                   stringsAsFactors = FALSE)
write.csv(tab2, file.path(tab_dir, "icb_external_gse91061_response.csv"), row.names = FALSE)

# 箱线图
df2 <- data.frame(value = s2, response = resp2)
p2 <- ggplot(df2, aes(x = response, y = value, fill = response)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.5, alpha = 0.55, color = "grey25") +
  scale_fill_manual(values = c("Non-responder" = "#C0392B", "Responder" = "#2E86AB"),
                    name = NULL) +
  labs(x = NULL, y = "Signature score (within-cohort)",
       title = "GSE91061 (Melanoma, anti-PD-1)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "icb_response_gse91061.png"), p2, width = 5, height = 5, dpi = 300)

# ---------- 汇总 ----------
cat("\n===== ICB 外部验证汇总 =====\n")
print(tab1, row.names = FALSE)
print(tab2, row.names = FALSE)
cat("\n说明：预后签名「分越高越差」。若签名能预测 ICB 应答，期望高分→不应答（GSE135222 HR>1、GSE91061 AUC>0.5）。\n")
cat("     脚本 31 发现层已得：7 程序签名不预测 ICB 应答（P=0.874）；本脚本为外部独立队列验证（预期亦为阴性，如实报告）。\n")
cat("已保存：tables/icb_external_gse135222_pfs.csv、tables/icb_external_gse91061_response.csv、\n")
cat("        figures/icb_pfs_km_gse135222.png、figures/icb_response_gse91061.png\n")
