# ============================================================
# 35_gene_level_bootstrap.R
# 基因级签名：bootstrap 乐观校正 + per-SD HR（修 34 的 HR_cont 尺度伪影）
#   1) 重新加载基因级矩阵（复用 34 逻辑，末尾存 gene_level_matrices.rds 供复用）
#   2) 1000× bootstrap 乐观校正（整条管线每次重跑：z-score → cv.glmnet → 选基因 → 风险分）
#      → 校正 C，与程序级 0.618 公平比
#   3) per-SD HR（风险分在队列内标准化后 Cox，尺度不变，4 队列可比；基因级 + 程序级并排）
# 依赖：survival / glmnet / org.Hs.eg.db / AnnotationDbi
# 运行：本地 RStudio（约 15-25 分钟，B=1000）
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

B <- 1000

# ---------- 工具（与 34 一致）----------
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
  tab <- read.delim(text = paste(lines[(tb + 1):(te - 1)], collapse = "\n"),
                    header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  id_col  <- grep("^ID$", colnames(tab), ignore.case = TRUE)[1]
  sym_col <- grep("symbol", colnames(tab), ignore.case = TRUE)[1]
  if (is.na(sym_col)) sym_col <- grep("gene", colnames(tab), ignore.case = TRUE)[1]
  id  <- as.character(tab[[id_col]]); sym <- as.character(tab[[sym_col]])
  sym <- trimws(sub("\\s*///.*$", "", sym))
  out <- sym; names(out) <- id; out[nzchar(out)]
}

# ---------- 1) 并集基因 ----------
prog_df <- read.csv(file.path(base, "program_genes_top50.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
lc <- read.csv(file.path(tab_dir, "lasso_coefficients.csv"), stringsAsFactors = FALSE)
sel_prog <- lc$program[lc$selected == TRUE]
sel_prog <- sel_prog[sel_prog %in% colnames(prog_df)]
union_genes <- unique(unlist(lapply(sel_prog, function(p) as.character(prog_df[[p]]))))
union_genes <- union_genes[!is.na(union_genes) & union_genes != ""]
cat("并集基因数:", length(union_genes), "\n")

# ---------- 2) 加载 4 队列基因矩阵 ----------
d <- readRDS(file.path(val_dir, "bulk_program_scores.rds"))
surv <- d$surv
gene_mat <- list()

cat("\n===== TCGA-LUAD =====\n")
tcga <- readRDS(file.path(val_dir, "TCGA-LUAD.rds"))
cnt <- tcga$counts; rn <- rownames(cnt)
if (grepl("^ENSG", rn[1])) {
  ens <- sub("\\..*$", "", rn)
  sym <- mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first")
  names(sym) <- rn; sym <- sym[!is.na(sym) & sym != ""]
} else sym <- setNames(rn, rn)
keep <- !is.na(sym[rownames(cnt)])
cnt <- cnt[keep, , drop = FALSE]; grp <- sym[rownames(cnt)]
tcga_sym <- collapse_mean(cnt, grp)
colnames(tcga_sym) <- substr(colnames(tcga_sym), 1, 15)
gene_mat[["TCGA-LUAD"]] <- tcga_sym[intersect(union_genes, rownames(tcga_sym)), , drop = FALSE]
cat("  TCGA 匹配并集基因:", nrow(gene_mat[["TCGA-LUAD"]]), "\n")

array_cohorts <- c(GSE68465 = "GPL96", GSE72094 = "GPL15048", GSE31210 = "GPL570")
for (g in names(array_cohorts)) {
  es <- readRDS(file.path(val_dir, paste0(g, "_eset.rds")))
  expr <- es$expr
  probe2sym <- read_gpl_soft(file.path(val_dir, paste0(array_cohorts[[g]], ".soft.gz")))
  keep <- !is.na(probe2sym[rownames(expr)])
  expr <- expr[keep, , drop = FALSE]; sym <- probe2sym[rownames(expr)]
  mat_sym <- collapse_mean(expr, sym)
  gene_mat[[g]] <- mat_sym[intersect(union_genes, rownames(mat_sym)), , drop = FALSE]
  cat("  ", g, " 匹配并集基因:", nrow(gene_mat[[g]]), "\n", sep = "")
}
saveRDS(gene_mat, file.path(val_dir, "gene_level_matrices.rds"))

# ---------- 3) 基因级 bootstrap（TCGA）----------
X <- gene_mat[["TCGA-LUAD"]]
s <- surv[["TCGA-LUAD"]]
common <- intersect(colnames(X), rownames(s))
X <- X[, common, drop = FALSE]                      # 基因 × 样本
Y <- Surv(s[common, "time"], s[common, "event"])
n <- ncol(X)
cat("\n训练样本:", n, " | 基因:", nrow(X), "\n")

fit_pipeline <- function(Xtr, Ytr) {                # Xtr: 基因 × 样本
  Xt <- t(Xtr)
  mu <- colMeans(Xt); sd <- apply(Xt, 2, sd); sd[sd == 0 | is.na(sd)] <- 1
  Xs <- sweep(sweep(Xt, 2, mu, "-"), 2, sd, "/")
  cv <- suppressWarnings(cv.glmnet(Xs, Ytr, family = "cox", alpha = 1,
                                   nfolds = 10, standardize = FALSE))
  fit <- suppressWarnings(glmnet(Xs, Ytr, family = "cox", alpha = 1,
                                 lambda = cv$lambda.min, standardize = FALSE))
  beta <- as.vector(coef(fit)); names(beta) <- colnames(Xs)
  list(mu = mu, sd = sd, beta = beta, nsel = sum(beta != 0))
}
c_index <- function(Xnew, Ynew, model) {
  if (model$nsel == 0) return(0.5)
  Xt <- t(Xnew)
  cmn <- intersect(names(model$beta), colnames(Xt))
  if (length(cmn) == 0) return(0.5)
  Xs <- sweep(sweep(Xt[, cmn, drop = FALSE], 2, model$mu[cmn], "-"), 2, model$sd[cmn], "/")
  risk <- drop(Xs %*% model$beta[cmn])
  if (sd(risk) == 0) return(0.5)
  fc <- coxph(Ynew ~ risk)
  as.numeric(summary(fc)$concordance[1])
}

set.seed(123)
m_full <- fit_pipeline(X, Y)
apparent <- c_index(X, Y, m_full)
cat("表观 C-index =", round(apparent, 4), " | 全样本选中基因数 =", m_full$nsel, "\n")

idx_all <- seq_len(n)
boot <- vector("list", B)
set.seed(20260908)
for (b in seq_len(B)) {
  ib <- sample(idx_all, n, replace = TRUE)
  Xb <- X[, ib, drop = FALSE]; Yb <- Y[ib]
  m <- fit_pipeline(Xb, Yb)
  c_boot <- c_index(Xb, Yb, m)
  c_orig <- c_index(X,  Y,  m)
  boot[[b]] <- c(boot_app = c_boot, orig = c_orig,
                 optimism = c_boot - c_orig, nsel = m$nsel)
  if (b %% 100 == 0) cat("  ...", b, "/", B, "\n")
}
res <- do.call(rbind, boot)

opt_mean <- mean(res[, "optimism"]); sd_opt <- sd(res[, "optimism"])
c_corrected <- apparent - opt_mean
c_med <- median(res[, "orig"])
c_lo <- as.numeric(quantile(res[, "orig"], 0.025))
c_hi <- as.numeric(quantile(res[, "orig"], 0.975))
n_null <- sum(res[, "nsel"] == 0)

cat("\n===== 基因级 bootstrap 乐观校正 =====\n")
cat("表观 C-index        =", round(apparent, 4), "\n")
cat("平均乐观度          =", round(opt_mean, 4), " (SD =", round(sd_opt, 4), ")\n")
cat("乐观校正 C-index    =", round(c_corrected, 4), "\n")
cat("校正后 95% CI(±1.96·SD) = [", round(c_corrected - 1.96 * sd_opt, 4), ",",
    round(c_corrected + 1.96 * sd_opt, 4), "]\n")
cat("bootstrap 检验 C 分布: 中位 =", round(c_med, 4),
    " | 2.5% =", round(c_lo, 4), " | 97.5% =", round(c_hi, 4), "\n")
cat("bootstrap 模型未坍缩(nsel≥1)比例 =", round(mean(res[, "nsel"] >= 1), 3), "\n")
cat("选中基因数分布:\n"); print(table(res[, "nsel"]))

# ---------- 4) per-SD HR（基因级 + 程序级）----------
risk_gene <- function(mat_gene, model) {
  Xt <- t(mat_gene)
  cmn <- intersect(names(model$beta), colnames(Xt))
  if (length(cmn) == 0) return(NULL)
  Xs <- sweep(sweep(Xt[, cmn, drop = FALSE], 2, model$mu[cmn], "-"), 2, model$sd[cmn], "/")
  drop(Xs %*% model$beta[cmn])
}
hr_per_sd <- function(risk, s) {
  common <- intersect(names(risk), rownames(s))
  risk <- risk[common]; s <- s[common, , drop = FALSE]
  z <- (risk - mean(risk)) / sd(risk)
  fit <- coxph(Surv(s$time, s$event) ~ z)
  ci <- exp(confint(fit))
  data.frame(HR_perSD = as.numeric(exp(coef(fit))), lower = ci[1], upper = ci[2],
             P = summary(fit)$coefficients[1, "Pr(>|z|)"])
}

gene_hr <- do.call(rbind, lapply(names(gene_mat), function(g) {
  r <- risk_gene(gene_mat[[g]], m_full)
  if (is.null(r)) return(NULL)
  data.frame(cohort = g, hr_per_sd(r, surv[[g]]))
}))

sig <- readRDS(file.path(res_dir, "rds", "signature_results.rds"))
prog_hr <- do.call(rbind, lapply(names(sig$risk), function(g) {
  data.frame(cohort = g, hr_per_sd(sig$risk[[g]], surv[[g]]))
}))

cmp_hr <- data.frame(
  cohort = gene_hr$cohort,
  HR_perSD_gene = round(gene_hr$HR_perSD, 2),
  HR_perSD_prog = round(prog_hr$HR_perSD[match(gene_hr$cohort, prog_hr$cohort)], 2),
  P_gene = signif(gene_hr$P, 3),
  P_prog = signif(prog_hr$P[match(gene_hr$cohort, prog_hr$cohort)], 3),
  stringsAsFactors = FALSE)
cat("\n===== per-SD HR（基因级 vs 程序级，尺度不变可比）=====\n")
print(cmp_hr, row.names = FALSE)

# ---------- 5) 保存 ----------
out <- list(apparent = apparent, optimism_mean = opt_mean, corrected = c_corrected,
            boot_res = res, HR_perSD = cmp_hr,
            gene_hr = gene_hr, prog_hr = prog_hr)
saveRDS(out, file.path(res_dir, "rds", "gene_level_bootstrap.rds"))
write.csv(as.data.frame(res), file.path(tab_dir, "gene_level_bootstrap.csv"), row.names = FALSE)
write.csv(cmp_hr, file.path(tab_dir, "gene_level_perSD_HR.csv"), row.names = FALSE)

cat("\n已保存: gene_level_bootstrap.rds / gene_level_bootstrap.csv / gene_level_perSD_HR.csv\n")
