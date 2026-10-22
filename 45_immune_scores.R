# ============================================================
# 45_immune_scores.R
# 目的：补「免疫预测评分」——用已发表的免疫评分量化高/低风险组的免疫表型
#   ① CYT 细胞溶解活性（Rooney 2015，GZMA+PRF1）——手算，最稳
#   ② IFN-γ 6 基因签名（Ayers 2017：IDO1/CXCL10/CXCL9/HLA-DRA/STAT1/IFNG）——手算
#   ③ IPS 免疫表型分（Charoentong 2017，162 基因 4 类）——用 hacksig 包（有完整基因表）
#   TIDE 需 Python，本脚本不做（见说明）。
# 输出：4 bulk 队列 High vs Low 风险 + 2 ICB 队列 vs 结局
# 依赖：hacksig（可选，仅 IPS）；org.Hs.eg.db / AnnotationDbi / survival / ggplot2
# 运行：本地 RStudio
# ============================================================

suppressPackageStartupMessages({
  library(org.Hs.eg.db); library(AnnotationDbi); library(survival); library(ggplot2)
})

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

has_ips <- requireNamespace("hacksig", quietly = TRUE)
if (!has_ips) cat("【缺包】hacksig -> install.packages('hacksig')（IPS 将跳过，CYT/IFN-γ 照跑）\n")

# ---------- 工具 ----------
collapse_mean <- function(mat, group) {
  mat  <- as.matrix(mat); g <- as.character(group)
  keep <- !is.na(g); mat <- mat[keep, , drop = FALSE]; g <- g[keep]
  gf <- factor(g)
  s <- rowsum(mat, gf, reorder = TRUE)
  n <- as.vector(table(gf)[rownames(s)])
  s / n
}
read_gpl_soft <- function(path) {
  con <- gzfile(path, "rt"); lines <- readLines(con, warn = FALSE); close(con)
  tb <- which(lines == "!platform_table_begin")[1]
  te <- which(lines == "!platform_table_end")[1]
  if (is.na(tb) || is.na(te)) stop("找不到平台表标记：", path)
  tab <- read.delim(text = paste(lines[(tb + 1):(te - 1)], collapse = "\n"),
                    header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  id_col  <- grep("^ID$", colnames(tab), ignore.case = TRUE)[1]
  sym_col <- grep("symbol", colnames(tab), ignore.case = TRUE)[1]
  if (is.na(sym_col)) sym_col <- grep("gene", colnames(tab), ignore.case = TRUE)[1]
  id  <- as.character(tab[[id_col]]); sym <- as.character(tab[[sym_col]])
  sym <- trimws(sub("\\s*///.*$", "", sym))
  out <- sym; names(out) <- id; out[nzchar(out)]
}
# 评分：基因跨样本 z-score 后取均值（与脚本 38 同口径）
score <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < 2) return(NULL)
  z <- t(scale(t(mat[g, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}
CYT_genes  <- c("GZMA", "PRF1")
IFNG_genes <- c("IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG")

ips_of <- function(mat) {
  if (!requireNamespace("hacksig", quietly = TRUE)) return(NULL)
  tryCatch({
    zmat <- t(scale(t(mat)))                                  # 逐基因 z-score 去量纲
    res  <- hacksig::hack_immunophenoscore(expr_data = zmat)  # 返回 tibble(sample_id, raw_score, ips_score)
    setNames(as.numeric(res$raw_score), as.character(res$sample_id))[colnames(mat)]
  }, error = function(e) { cat("  IPS 失败：", conditionMessage(e), "\n"); NULL })
}

# ---------- 1) 4 bulk 队列 symbol 矩阵 ----------
cat("===== 构建 bulk symbol 矩阵 =====\n")
tcga <- readRDS(file.path(val_dir, "TCGA-LUAD.rds"))
cnt <- tcga$counts; rn <- rownames(cnt)
ens <- sub("\\..*$", "", rn)
sym <- mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first")
names(sym) <- rn; sym <- sym[!is.na(sym) & sym != ""]
keep <- !is.na(sym[rownames(cnt)])
cnt_sym <- collapse_mean(cnt[keep, , drop = FALSE], sym[rownames(cnt)[keep]])
tcga_sym <- log2(sweep(cnt_sym, 2, colSums(cnt_sym), "/") * 1e6 + 1)  # log2(CPM+1) 归一化
colnames(tcga_sym) <- substr(colnames(tcga_sym), 1, 15)

array_cohorts <- c(GSE68465 = "GPL96", GSE72094 = "GPL15048", GSE31210 = "GPL570")
geo_sym <- list()
for (g in names(array_cohorts)) {
  es <- readRDS(file.path(val_dir, paste0(g, "_eset.rds")))
  expr <- es$expr
  probe2sym <- read_gpl_soft(file.path(val_dir, paste0(array_cohorts[[g]], ".soft.gz")))
  keep <- !is.na(probe2sym[rownames(expr)])
  geo_sym[[g]] <- collapse_mean(expr[keep, , drop = FALSE], probe2sym[rownames(expr)[keep]])
}
bulk <- c(list("TCGA-LUAD" = tcga_sym), geo_sym)

sig <- readRDS(file.path(res_dir, "rds", "signature_results.rds"))
risk_list <- sig$risk

bulk_rows <- list()
for (g in names(bulk)) {
  mat <- bulk[[g]]
  risk <- risk_list[[g]]
  common <- intersect(colnames(mat), names(risk))
  mat <- mat[, common, drop = FALSE]; risk <- risk[common]
  grp <- factor(ifelse(risk > median(risk), "High", "Low"), levels = c("Low", "High"))
  cyt <- score(mat, CYT_genes); ifn <- score(mat, IFNG_genes); ips <- ips_of(mat)
  mk <- function(v, nm) {
    if (is.null(v)) return(c(NA, NA, NA))
    c(median(v[grp == "Low"]), median(v[grp == "High"]), wilcox.test(v ~ grp)$p.value)
  }
  r_cyt <- mk(cyt, "CYT"); r_ifn <- mk(ifn, "IFNG"); r_ips <- mk(ips, "IPS")
  bulk_rows[[g]] <- data.frame(cohort = g,
    CYT_Low = r_cyt[1], CYT_High = r_cyt[2], CYT_P = r_cyt[3],
    IFNg_Low = r_ifn[1], IFNg_High = r_ifn[2], IFNg_P = r_ifn[3],
    IPS_Low = r_ips[1], IPS_High = r_ips[2], IPS_P = r_ips[3],
    stringsAsFactors = FALSE)
}
bulk_tab <- do.call(rbind, bulk_rows)
cat("\n===== 免疫评分：High vs Low 风险（4 bulk）=====\n")
print(bulk_tab, row.names = FALSE)
write.csv(bulk_tab, file.path(tab_dir, "immune_scores_bulk_by_risk.csv"), row.names = FALSE)

# ---------- 2) 2 ICB 队列 symbol 矩阵 + 评分 vs 结局 ----------
cat("\n===== ICB 队列免疫评分 vs 结局 =====\n")
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
    out[[sub(":.*$", "", v[1])]] <- sub("^[^:]*:\\s*", "", v)
  }
  out
}

# GSE135222（NSCLC，PFS）
dir1 <- file.path(base, "GSE135222")
exp1 <- read.delim(gzfile(file.path(dir1, "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz")),
                   row.names = 1, check.names = FALSE)
sym1 <- mapIds(org.Hs.eg.db, sub("\\..*$", "", rownames(exp1)), "SYMBOL", "ENSEMBL", multiVals = "first")
keep1 <- !is.na(sym1)
mat1 <- log2(collapse_mean(as.matrix(exp1[keep1, , drop = FALSE]) + 1, sym1[keep1]) + 1)
clin1 <- read_series_matrix(file.path(dir1, "GSE135222_series_matrix.txt.gz"))
clin1$key <- gsub(" ", "", clin1$title)
clin1$pfs <- as.numeric(clin1[["progression-free survival (pfs)"]])
clin1$pfs_time <- as.numeric(clin1[["pfs.time"]])
cm1 <- intersect(colnames(mat1), clin1$key)
cl1 <- clin1[match(cm1, clin1$key), ]; rownames(cl1) <- cl1$key
m1 <- mat1[, cm1]
cyt1 <- score(m1, CYT_genes); ifn1 <- score(m1, IFNG_genes); ips1 <- ips_of(m1)
ok1 <- !is.na(cl1$pfs_time) & !is.na(cl1$pfs) & cl1$pfs_time > 0
icb1 <- data.frame(cohort = "GSE135222", score = c("CYT", "IFNg", "IPS"),
  med_event = NA, med_censored = NA, P = NA)
# 连续评分 → PFS Cox（事件=1 进展）
for (nm in c("CYT", "IFNg", "IPS")) {
  v <- switch(nm, CYT = cyt1, IFNg = ifn1, IPS = ips1)
  if (is.null(v)) next
  fit <- coxph(Surv(cl1$pfs_time[ok1], cl1$pfs[ok1]) ~ v[ok1])
  icb1[icb1$score == nm, "P"] <- signif(summary(fit)$coefficients[1, "Pr(>|z|)"], 3)
  icb1[icb1$score == nm, "med_event"] <- round(median(v[ok1][cl1$pfs[ok1] == 1]), 3)
  icb1[icb1$score == nm, "med_censored"] <- round(median(v[ok1][cl1$pfs[ok1] == 0]), 3)
}
cat("GSE135222（NSCLC，免疫评分 vs PFS，Cox P）:\n"); print(icb1, row.names = FALSE)

# GSE91061（黑色素瘤，response，Pre）
dir2 <- file.path(base, "GSE91061")
exp2 <- read.csv(gzfile(file.path(dir2, "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz")),
                 row.names = 1, check.names = FALSE)
sym2 <- mapIds(org.Hs.eg.db, rownames(exp2), "SYMBOL", "ENTREZID", multiVals = "first")
keep2 <- !is.na(sym2)
mat2 <- log2(collapse_mean(as.matrix(exp2[keep2, , drop = FALSE]) + 1, sym2[keep2]) + 1)
clin2 <- read_series_matrix(file.path(dir2, "GSE91061_series_matrix.txt.gz"))
clin2$response <- clin2[["response"]]; clin2$visit <- clin2[["visit (pre or on treatment)"]]
cm2 <- intersect(colnames(mat2), clin2$title)
cl2 <- clin2[match(cm2, clin2$title), ]; rownames(cl2) <- cl2$title
m2 <- mat2[, cm2]
cyt2 <- score(m2, CYT_genes); ifn2 <- score(m2, IFNG_genes); ips2 <- ips_of(m2)
use2 <- cl2$visit == "Pre" & cl2$response %in% c("PRCR", "PD", "SD")
resp2 <- factor(ifelse(cl2$response[use2] == "PRCR", "Responder", "Non-responder"),
                levels = c("Non-responder", "Responder"))
icb2 <- data.frame(cohort = "GSE91061", score = c("CYT", "IFNg", "IPS"),
  med_NonResp = NA, med_Resp = NA, P = NA)
for (nm in c("CYT", "IFNg", "IPS")) {
  v <- switch(nm, CYT = cyt2, IFNg = ifn2, IPS = ips2)
  if (is.null(v)) next
  vv <- v[cm2][use2]
  w <- wilcox.test(vv ~ resp2)
  icb2[icb2$score == nm, "P"] <- signif(w$p.value, 3)
  icb2[icb2$score == nm, "med_NonResp"] <- round(median(vv[resp2 == "Non-responder"]), 3)
  icb2[icb2$score == nm, "med_Resp"] <- round(median(vv[resp2 == "Responder"]), 3)
}
cat("\nGSE91061（黑色素瘤，免疫评分 vs 应答，Wilcoxon P）:\n"); print(icb2, row.names = FALSE)

write.csv(icb1, file.path(tab_dir, "immune_scores_icb_gse135222.csv"), row.names = FALSE)
write.csv(icb2, file.path(tab_dir, "immune_scores_icb_gse91061.csv"), row.names = FALSE)

cat("\n已保存 immune_scores_bulk_by_risk.csv / immune_scores_icb_gse135222.csv / immune_scores_icb_gse91061.csv\n")
cat("说明：CYT(Rooney)=GZMA+PRF1；IFN-γ(Ayers)=IDO1/CXCL10/CXCL9/HLA-DRA/STAT1/IFNG；IPS(Charoentong)用 hacksig。\n")
cat("     TIDE 未做（需 Python tidepy 或 TIDE 官网，见对话说明）；EaSIeR 基因表待核实，暂未加入。\n")
