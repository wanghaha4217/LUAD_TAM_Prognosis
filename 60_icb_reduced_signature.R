# ============================================================
# 60_icb_reduced_signature.R
# 目的：审稿人建议 —— 只用不良轴 P8/P9/P4（去掉保护轴）重测 ICB 预测
#   假设：若「预后≠预测」是保护轴稀释了不良轴的 ICB 信号，只留不良轴也许能预测 ICB。
#   预期仍是阴性（免疫冷），但这是审稿人要的对照，如实报告即可。
#   程序编号 = 原论文 NMF k=10（P8 缺氧 / P9 热休克 / P4 S100A8-9）。
#   签名分两种口径：(a) LASSO β 加权  (b) 等权 z 和。
# 依赖：org.Hs.eg.db / AnnotationDbi / survival / ggplot2
# 运行：本地 RStudio（约 1-2 分钟）
# ============================================================

suppressPackageStartupMessages({
  library(org.Hs.eg.db); library(AnnotationDbi)
  library(survival); library(ggplot2)
})

base    <- "."
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables"); fig_dir <- file.path(res_dir, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 工具（同脚本 39） ----------
collapse_mean <- function(mat, group) {
  mat <- as.matrix(mat); g <- as.character(group)
  keep <- !is.na(g); mat <- mat[keep, , drop = FALSE]; g <- g[keep]
  gf <- factor(g); s <- rowsum(mat, gf, reorder = TRUE)
  n <- as.vector(table(gf)[rownames(s)]); s / n
}
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
    v <- parts[[i]][-1]; v <- gsub('^"|"$', '', v)
    nm <- sub(":.*$", "", v[1]); out[[nm]] <- sub("^[^:]*:\\s*", "", v)
  }
  out
}
score_program <- function(mat_sym, genes) {
  g <- intersect(genes, rownames(mat_sym))
  if (length(g) < 3) return(rep(NA_real_, ncol(mat_sym)))
  z <- t(scale(t(mat_sym[g, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}
# 缩减签名：只取指定程序，β 加权 或 等权（z 后求和）
reduced_of <- function(score_mat, progs, betas = NULL) {
  s <- rep(0, nrow(score_mat)); names(s) <- rownames(score_mat)
  for (i in seq_along(progs)) {
    pp <- progs[i]
    if (!pp %in% colnames(score_mat)) next
    z <- scale(score_mat[, pp])[, 1]; z[is.na(z)] <- 0
    w <- if (is.null(betas)) 1 else betas[i]
    s <- s + w * z
  }
  s
}

# ---------- 程序基因 + LASSO β ----------
prog_df <- read.csv(file.path(base, "program_genes_top50.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
lc <- read.csv(file.path(tab_dir, "lasso_coefficients.csv"), stringsAsFactors = FALSE)
lc <- lc[lc$selected == TRUE, c("program", "beta")]

AXIS <- c("P8", "P9", "P4")                     # 不良轴
betas <- setNames(lc$beta[match(AXIS, lc$program)], AXIS)
betas[is.na(betas)] <- 0
cat("缩减签名系数（P8/P9/P4）:\n"); print(betas)

# ============================================================
# 1. GSE135222（NSCLC anti-PD-1/PD-L1，PFS）
# ============================================================
cat("\n===== GSE135222（NSCLC，PFS）=====\n")
dir1 <- file.path(base, "GSE135222")
exp1 <- read.delim(gzfile(file.path(dir1, "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz")),
                   row.names = 1, check.names = FALSE)
ens1 <- sub("\\..*$", "", rownames(exp1))
sym1 <- mapIds(org.Hs.eg.db, ens1, "SYMBOL", "ENSEMBL", multiVals = "first")
keep1 <- !is.na(sym1)
mat1 <- collapse_mean(log2(as.matrix(exp1[keep1, , drop = FALSE]) + 1), sym1[keep1])

ps1 <- sapply(AXIS, function(p) score_program(mat1, prog_df[[p]]))
rownames(ps1) <- colnames(mat1)
sig1b <- reduced_of(ps1, AXIS, betas)           # β 加权
sig1e <- reduced_of(ps1, AXIS, NULL)            # 等权

clin1 <- read_series_matrix(file.path(dir1, "GSE135222_series_matrix.txt.gz"))
clin1$sample_key <- gsub(" ", "", clin1$title)
clin1$pfs      <- as.numeric(clin1[["progression-free survival (pfs)"]])
clin1$pfs_time <- as.numeric(clin1[["pfs.time"]])
common1 <- intersect(colnames(mat1), clin1$sample_key)
idx1 <- match(common1, clin1$sample_key)
cl1 <- clin1[idx1, ]; rownames(cl1) <- cl1$sample_key; cl1 <- cl1[common1, , drop = FALSE]
sig1b <- sig1b[common1]; sig1e <- sig1e[common1]
ok1 <- !is.na(sig1b) & !is.na(cl1$pfs_time) & !is.na(cl1$pfs) & cl1$pfs_time > 0
t1 <- cl1$pfs_time[ok1]; e1 <- cl1$pfs[ok1]

report_cox <- function(tag, s) {
  fit <- coxph(Surv(t1, e1) ~ s[ok1])
  hr <- exp(coef(fit)); p <- summary(fit)$coefficients[1, "Pr(>|z|)"]
  ci <- as.numeric(summary(fit)$concordance[1])
  cat("  ", tag, ": HR =", round(hr, 3), "| P =", signif(p, 3), "| C =", round(ci, 3), "\n", sep = "")
  data.frame(cohort = "GSE135222", signature = tag, N = sum(ok1), events = sum(e1),
             HR_cont = round(hr, 3), P_cont = signif(p, 3), Cindex = round(ci, 3))
}
tab1 <- rbind(report_cox("P8P9P4_beta", sig1b), report_cox("P8P9P4_equal", sig1e))

# ============================================================
# 2. GSE91061（黑色素瘤 anti-PD-1，response，Pre）
# ============================================================
cat("\n===== GSE91061（黑色素瘤，response）=====\n")
dir2 <- file.path(base, "GSE91061")
exp2 <- read.csv(gzfile(file.path(dir2, "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz")),
                 row.names = 1, check.names = FALSE)
sym2 <- mapIds(org.Hs.eg.db, rownames(exp2), "SYMBOL", "ENTREZID", multiVals = "first")
keep2 <- !is.na(sym2)
mat2 <- collapse_mean(log2(as.matrix(exp2[keep2, , drop = FALSE]) + 1), sym2[keep2])

ps2 <- sapply(AXIS, function(p) score_program(mat2, prog_df[[p]]))
rownames(ps2) <- colnames(mat2)
sig2b <- reduced_of(ps2, AXIS, betas); sig2e <- reduced_of(ps2, AXIS, NULL)

clin2 <- read_series_matrix(file.path(dir2, "GSE91061_series_matrix.txt.gz"))
clin2$response <- clin2[["response"]]
clin2$visit    <- clin2[["visit (pre or on treatment)"]]
common2 <- intersect(colnames(mat2), clin2$title)
idx2 <- match(common2, clin2$title)
cl2 <- clin2[idx2, ]; rownames(cl2) <- cl2$title; cl2 <- cl2[common2, , drop = FALSE]
sig2b <- sig2b[common2]; sig2e <- sig2e[common2]
use2 <- cl2$visit == "Pre" & cl2$response %in% c("PRCR", "PD", "SD") & !is.na(sig2b)
resp2 <- factor(ifelse(cl2$response[use2] == "PRCR", "Responder", "Non-responder"),
                levels = c("Non-responder", "Responder"))

report_auc <- function(tag, s) {
  wr <- wilcox.test(s[use2] ~ resp2)
  n_r <- sum(resp2 == "Responder"); n_nr <- sum(resp2 == "Non-responder")
  auc <- as.numeric(wr$statistic) / (n_r * n_nr)   # >0.5 = 高分→不应答
  cat("  ", tag, ": P =", signif(wr$p.value, 3), "| AUC(高分→不应答) =", round(auc, 3), "\n", sep = "")
  data.frame(cohort = "GSE91061", signature = tag, N = length(resp2),
             responders = sum(resp2 == "Responder"),
             AUC_high_to_NonResponder = round(auc, 3), Wilcoxon_P = signif(wr$p.value, 3))
}
tab2 <- rbind(report_auc("P8P9P4_beta", sig2b), report_auc("P8P9P4_equal", sig2e))

# ---------- 汇总 ----------
cat("\n===== 缩减签名（P8/P9/P4）ICB 预测汇总 =====\n")
print(tab1, row.names = FALSE)
print(tab2, row.names = FALSE)
write.csv(tab1, file.path(tab_dir, "icb_reduced_signature.csv"), row.names = FALSE)
write.csv(tab2, file.path(tab_dir, "icb_reduced_signature_gse91061.csv"), row.names = FALSE)
cat("\n说明：若 β 加权与等权的 P8/P9/P4 缩减签名在 ICB 队列仍不显著（预期），\n")
cat("      则进一步坐实「预后签名≠预测签名」：不良轴单独也无 ICB 预测力，与免疫冷表型一致。\n")
