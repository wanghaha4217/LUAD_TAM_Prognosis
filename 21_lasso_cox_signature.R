# ============================================================
# 21_lasso_cox_signature.R
# 目的：LASSO-Cox 建 TAM 程序预后签名，并在 3 个独立队列验证
#   训练：TCGA-LUAD（9 个程序打分 → LASSO 选程序 → 线性风险分）
#   验证：GSE68465 / GSE72094 / GSE31210（用 TCGA 的系数算同一风险分）
# 依赖：survival / glmnet（若报"没有 glmnet 包"，先 install.packages("glmnet")）
# 运行：本地 RStudio（依赖 20 脚本输出的 bulk_program_scores.rds）
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
})

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
dir.create(file.path(res_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(res_dir, "rds"),    recursive = TRUE, showWarnings = FALSE)

d <- readRDS(file.path(val_dir, "bulk_program_scores.rds"))
scores   <- d$scores
surv     <- d$surv
coverage <- d$coverage

# ============================================================
# 1. 训练：TCGA-LUAD LASSO-Cox
# ============================================================
cat("===== 1. 训练（TCGA-LUAD）=====\n")
tr_scores <- as.matrix(scores[["TCGA-LUAD"]])
tr_surv   <- surv[["TCGA-LUAD"]]
tr_y <- Surv(tr_surv$time, tr_surv$event)

# 用「训练集的 mean/sd」把各程序分标准化（跨队列统一尺度），再喂给 glmnet
tr_mean <- colMeans(tr_scores)
tr_sd   <- apply(tr_scores, 2, sd)

scale_scores <- function(df, mu, sd) {
  df  <- as.matrix(df)
  cmn <- intersect(colnames(df), names(mu))
  sweep(sweep(df[, cmn, drop = FALSE], 2, mu[cmn], "-"), 2, sd[cmn], "/")
}
tr_x <- scale_scores(tr_scores, tr_mean, tr_sd)

set.seed(123)
cv <- cv.glmnet(tr_x, tr_y, family = "cox", alpha = 1, nfolds = 10,
                standardize = FALSE)   # 已在上面手动标准化，这里关掉内部标准化
cat("lambda.min =", round(cv$lambda.min, 4),
    " | lambda.1se =", round(cv$lambda.1se, 4), "\n")

coef_at <- function(lam) {
  fit <- glmnet(tr_x, tr_y, family = "cox", alpha = 1, lambda = lam,
                standardize = FALSE)
  b <- as.vector(coef(fit)); names(b) <- colnames(tr_x); b
}
coef_min <- coef_at(cv$lambda.min)
coef_1se <- coef_at(cv$lambda.1se)
cat("lambda.min 选中:", paste(names(coef_min)[coef_min != 0], collapse = ", "), "\n")
cat("lambda.1se 选中:", paste(names(coef_1se)[coef_1se != 0], collapse = ", "), "\n")

# 主签名用 lambda.min；若全被压成 0 则退回 1se
coef_final <- coef_min
if (sum(coef_final != 0) == 0) coef_final <- coef_1se
sel <- names(coef_final)[coef_final != 0]
cat("\n最终签名程序（", length(sel), " 个）:", paste(sel, collapse = ", "), "\n", sep = "")
print(round(coef_final[sel], 4))

# ============================================================
# 2. 风险分 + 各队列评估
# ============================================================
# 风险分 = Σ β_i × (标准化后的程序分)；所有队列都用 TCGA 的 mean/sd 标准化
risk_of <- function(score_df, coefs, mu, sd) {
  xx <- scale_scores(score_df, mu, sd)
  cmn <- intersect(names(coefs), colnames(xx))
  as.numeric(xx[, cmn, drop = FALSE] %*% coefs[cmn])
}

eval_cohort <- function(name, score_df, surv_df, coefs, mu, sd) {
  risk <- risk_of(score_df, coefs, mu, sd)
  ok <- !is.na(risk) & !is.na(surv_df$time) & !is.na(surv_df$event)
  risk <- risk[ok]; s <- surv_df[ok, , drop = FALSE]
  yy <- Surv(s$time, s$event)
  n <- length(risk); ev <- sum(s$event == 1)

  # 连续风险分 HR（每 +1 单位）；C-index 直接取 coxph 自带的（方向与 HR 一致）
  fc  <- coxph(yy ~ risk)
  ci  <- as.numeric(summary(fc)$concordance[1])   # Harrell's C
  hr_c <- exp(coef(fc)); p_c <- summary(fc)$coefficients[1, "Pr(>|z|)"]

  # 中位二分（High/Low）+ log-rank
  med <- median(risk)
  grp <- factor(ifelse(risk > med, "High", "Low"), levels = c("Low", "High"))
  hr_b <- NA_real_; p_b <- NA_real_; lr_p <- NA_real_
  if (length(unique(grp)) == 2) {
    lr_p <- survdiff(yy ~ grp)$pvalue
    fb   <- coxph(yy ~ grp)
    hr_b <- exp(coef(fb)); p_b <- summary(fb)$coefficients[1, "Pr(>|z|)"]
  }

  data.frame(cohort = name, N = n, events = ev,
             Cindex = round(ci, 3),
             HR_cont = round(hr_c, 3), P_cont = signif(p_c, 3),
             HR_highVSlow = round(hr_b, 3), P_binary = signif(p_b, 3),
             logrank_P = signif(lr_p, 3),
             stringsAsFactors = FALSE)
}

res <- do.call(rbind, lapply(names(scores), function(g)
  eval_cohort(g, scores[[g]], surv[[g]], coef_final, tr_mean, tr_sd)))
cat("\n===== 2. 各队列验证（TCGA 为训练集内评估，偏乐观）=====\n")
print(res, row.names = FALSE)

# ============================================================
# 3. 保存
# ============================================================
risk_list <- lapply(names(scores), function(g) {
  r <- risk_of(scores[[g]], coef_final, tr_mean, tr_sd)
  names(r) <- rownames(scores[[g]]); r
})
names(risk_list) <- names(scores)

coef_df <- data.frame(program = names(coef_final),
                      beta    = round(coef_final, 4),
                      selected = names(coef_final) %in% sel,
                      stringsAsFactors = FALSE)

write.csv(coef_df, file.path(res_dir, "tables", "lasso_coefficients.csv"), row.names = FALSE)
write.csv(res,     file.path(res_dir, "tables", "signature_validation.csv"),  row.names = FALSE)
saveRDS(list(coef = coef_final, sel = sel, risk = risk_list, validation = res,
             tr_mean = tr_mean, tr_sd = tr_sd,
             lambda.min = cv$lambda.min, lambda.1se = cv$lambda.1se),
        file.path(res_dir, "rds", "signature_results.rds"))
cat("\n已保存:", res_dir, "\n")
