# ============================================================
# 34_gene_level_lasso.R
# 补充分析：基因级 LASSO-Cox 签名（vs 程序级签名，公平对比）
#   思路：取 7 个选中程序(P2/P3/P4/P5/P8/P9/P10)的 top-50 基因并集 → TCGA 基因级 LASSO
#         → 选基因(带系数/符号) → 3 GEO 队列验证 → 与程序级并排对比
#   注意：基因必须保留方向(系数符号)——不能无符号平均(否则 P8 不良 / P3 保护 相互抵消)
# 依赖：survival / glmnet / org.Hs.eg.db / AnnotationDbi（GPL 注释读本地 .soft.gz，不联网）
# 运行：本地 RStudio（约 2-5 分钟）
# 输出：
#   results/tables/gene_level_lasso_coefficients.csv
#   results/tables/gene_level_validation.csv
#   results/rds/gene_level_lasso.rds
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 工具（与 20 脚本一致）----------
collapse_mean <- function(mat, group) {
  mat  <- as.matrix(mat); g <- as.character(group)
  keep <- !is.na(g); mat <- mat[keep, , drop = FALSE]; g <- g[keep]
  gf <- factor(g)
  s <- rowsum(mat, gf, reorder = TRUE)
  n <- as.vector(table(gf)[rownames(s)])
  s / n
}
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# 本地读 GPL .soft.gz（ID → symbol），不联网
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
  id  <- as.character(tab[[id_col]])
  sym <- as.character(tab[[sym_col]])
  sym <- trimws(sub("\\s*///.*$", "", sym))
  out <- sym; names(out) <- id
  out[nzchar(out)]
}

# ---------- 1) 7 选中程序的基因并集 ----------
prog_df <- read.csv(file.path(base, "program_genes_top50.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
lc <- read.csv(file.path(tab_dir, "lasso_coefficients.csv"), stringsAsFactors = FALSE)
sel_prog <- lc$program[lc$selected == TRUE]
sel_prog <- sel_prog[sel_prog %in% colnames(prog_df)]
cat("选中程序:", paste(sel_prog, collapse = ", "), "\n")

union_genes <- unique(unlist(lapply(sel_prog, function(p) as.character(prog_df[[p]]))))
union_genes <- union_genes[!is.na(union_genes) & union_genes != ""]
cat("并集基因数:", length(union_genes), "\n")

# ---------- 2) 加载 4 队列的 symbol 基因矩阵（仅保留并集基因）----------
d <- readRDS(file.path(val_dir, "bulk_program_scores.rds"))
surv <- d$surv

gene_mat <- list()

# --- TCGA（Ensembl → symbol）---
cat("\n===== TCGA-LUAD 基因矩阵 =====\n")
tcga <- readRDS(file.path(val_dir, "TCGA-LUAD.rds"))
cnt <- tcga$counts
rn <- rownames(cnt)
if (grepl("^ENSG", rn[1])) {
  ens <- sub("\\..*$", "", rn)
  sym <- mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first")
  names(sym) <- rn
  sym <- sym[!is.na(sym) & sym != ""]
} else sym <- setNames(rn, rn)
keep <- !is.na(sym[rownames(cnt)])
cnt <- cnt[keep, , drop = FALSE]; grp <- sym[rownames(cnt)]
tcga_sym <- collapse_mean(cnt, grp)
colnames(tcga_sym) <- substr(colnames(tcga_sym), 1, 15)
g_tcga <- intersect(union_genes, rownames(tcga_sym))
gene_mat[["TCGA-LUAD"]] <- tcga_sym[g_tcga, , drop = FALSE]
cat("  TCGA 可匹配并集基因:", length(g_tcga), "\n")

# --- GEO（探针 → symbol，本地 GPL）---
array_cohorts <- c(GSE68465 = "GPL96", GSE72094 = "GPL15048", GSE31210 = "GPL570")
for (g in names(array_cohorts)) {
  cat("\n===== ", g, " 基因矩阵 =====\n", sep = "")
  es <- readRDS(file.path(val_dir, paste0(g, "_eset.rds")))
  expr <- es$expr
  probe2sym <- read_gpl_soft(file.path(val_dir, paste0(array_cohorts[[g]], ".soft.gz")))
  keep <- !is.na(probe2sym[rownames(expr)])
  expr <- expr[keep, , drop = FALSE]; sym <- probe2sym[rownames(expr)]
  mat_sym <- collapse_mean(expr, sym)
  g_geo <- intersect(union_genes, rownames(mat_sym))
  gene_mat[[g]] <- mat_sym[g_geo, , drop = FALSE]
  cat("  ", g, " 可匹配并集基因:", length(g_geo), "\n", sep = "")
}

# ---------- 3) TCGA 基因级 LASSO ----------
cat("\n===== TCGA 基因级 LASSO =====\n")
tcga_surv <- surv[["TCGA-LUAD"]]
common <- intersect(colnames(gene_mat[["TCGA-LUAD"]]), rownames(tcga_surv))
X <- t(gene_mat[["TCGA-LUAD"]][, common, drop = FALSE])   # 样本 × 基因
Y <- Surv(tcga_surv[common, "time"], tcga_surv[common, "event"])
cat("  TCGA 训练样本:", length(common), " | 基因:", ncol(X), "\n")

mu <- colMeans(X); sd <- apply(X, 2, sd); sd[sd == 0 | is.na(sd)] <- 1
Xs <- sweep(sweep(X, 2, mu, "-"), 2, sd, "/")
set.seed(123)
cv <- cv.glmnet(Xs, Y, family = "cox", alpha = 1, nfolds = 10, standardize = FALSE)
fit <- glmnet(Xs, Y, family = "cox", alpha = 1, lambda = cv$lambda.min, standardize = FALSE)
beta <- as.vector(coef(fit)); names(beta) <- colnames(Xs)
sel_genes <- names(beta)[beta != 0]
cat("  lambda.min 选中基因数:", length(sel_genes), "\n")
cat("  系数（按 |β| 降序）:\n")
print(round(sort(beta[sel_genes], decreasing = TRUE), 4))

# ---------- 4) 风险分 + 4 队列评估 ----------
risk_of <- function(mat_gene, beta, mu, sd) {
  Xs <- t(mat_gene)
  cmn <- intersect(names(beta), colnames(Xs))
  if (length(cmn) == 0) return(NULL)
  Xs <- Xs[, cmn, drop = FALSE]                              # 先子集到 cmn 基因
  Xs <- sweep(sweep(Xs, 2, mu[cmn], "-"), 2, sd[cmn], "/")   # 用 cmn 对应的 mean/sd
  drop(Xs %*% beta[cmn])                                     # drop 保留样本行名（关键）
}

eval_cohort <- function(g, beta, mu, sd) {
  r <- risk_of(gene_mat[[g]], beta, mu, sd)
  if (is.null(r)) return(NULL)
  s <- surv[[g]]
  common <- intersect(names(r), rownames(s))
  r <- r[common]; s <- s[common, , drop = FALSE]
  yy <- Surv(s$time, s$event)
  fc <- coxph(yy ~ r)
  ci <- as.numeric(summary(fc)$concordance[1])
  hr_c <- exp(coef(fc)); p_c <- summary(fc)$coefficients[1, "Pr(>|z|)"]
  med <- median(r)
  grp <- factor(ifelse(r > med, "High", "Low"), levels = c("Low", "High"))
  hr_b <- NA_real_; lr_p <- NA_real_
  if (length(unique(grp)) == 2) {
    lr_p <- survdiff(yy ~ grp)$pvalue
    hr_b <- exp(coef(coxph(yy ~ grp)))
  }
  data.frame(cohort = g, N = length(common), events = sum(s$event == 1),
             Cindex = round(ci, 3), HR_cont = round(hr_c, 3), P_cont = signif(p_c, 3),
             HR_highVSlow = round(hr_b, 3), logrank_P = signif(lr_p, 3),
             stringsAsFactors = FALSE)
}

gene_res <- do.call(rbind, lapply(names(gene_mat), function(g)
  eval_cohort(g, beta, mu, sd)))
cat("\n===== 基因级签名：4 队列验证 =====\n"); print(gene_res, row.names = FALSE)

# ---------- 5) 与程序级签名对比 ----------
prog_ref <- read.csv(file.path(tab_dir, "signature_validation.csv"), stringsAsFactors = FALSE)
cmp <- data.frame(
  cohort = gene_res$cohort,
  Cindex_gene = gene_res$Cindex,
  Cindex_prog = prog_ref$Cindex[match(gene_res$cohort, prog_ref$cohort)],
  HR_gene = gene_res$HR_highVSlow,
  HR_prog = prog_ref$HR_highVSlow[match(gene_res$cohort, prog_ref$cohort)],
  logrankP_gene = gene_res$logrank_P,
  logrankP_prog = prog_ref$logrank_P[match(gene_res$cohort, prog_ref$cohort)],
  stringsAsFactors = FALSE)
cat("\n===== 基因级 vs 程序级（C-index / HR 高对低 / log-rank P）=====\n")
print(cmp, row.names = FALSE)

# ---------- 6) 保存 ----------
coef_df <- data.frame(gene = sel_genes, beta = round(beta[sel_genes], 4),
                      stringsAsFactors = FALSE)
coef_df <- coef_df[order(-abs(coef_df$beta)), ]
write.csv(coef_df, file.path(tab_dir, "gene_level_lasso_coefficients.csv"), row.names = FALSE)
write.csv(gene_res, file.path(tab_dir, "gene_level_validation.csv"), row.names = FALSE)
saveRDS(list(union_genes = union_genes, sel_genes = sel_genes, beta = beta,
             mu = mu, sd = sd, validation = gene_res, comparison = cmp),
        file.path(res_dir, "rds", "gene_level_lasso.rds"))

cat("\n已保存: gene_level_lasso_coefficients.csv / gene_level_validation.csv / gene_level_lasso.rds\n")
cat("\n提示：基因级 LASSO 选择可能不稳定（不同 CV 折选出的基因会漂），务必以 3 外部队列验证为准判断，而非只看 TCGA 表观 C。\n")
