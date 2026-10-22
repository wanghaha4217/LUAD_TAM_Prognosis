# ============================================================
# 27_lasso_figure.R
# LASSO-Cox 方法图（B）：三张面板
#   Fig2a 交叉验证曲线（deviance vs log(lambda)，标 lambda.min / lambda.1se）
#   Fig2b 系数路径（各程序系数随 lambda 收缩到 0）
#   Fig2c 最终 7 程序系数柱状图（红=不良，蓝=保护）
# 依赖：survival / glmnet / ggplot2
# 运行：本地 RStudio（复用 21 脚本的训练设置，seed=123 完全一致）
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(ggplot2)
})

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
fig_dir <- file.path(res_dir, "figures")

d   <- readRDS(file.path(val_dir, "bulk_program_scores.rds"))
sig <- readRDS(file.path(res_dir, "rds", "signature_results.rds"))

# ---- 重跑 cv.glmnet（与 21 脚本完全一致，快速）----
tr_scores <- as.matrix(d$scores[["TCGA-LUAD"]])
tr_surv   <- d$surv[["TCGA-LUAD"]]
tr_y      <- Surv(tr_surv$time, tr_surv$event)
tr_mean   <- colMeans(tr_scores)
tr_sd     <- apply(tr_scores, 2, sd)

scale_scores <- function(df, mu, sd) {
  df  <- as.matrix(df)
  cmn <- intersect(colnames(df), names(mu))
  sweep(sweep(df[, cmn, drop = FALSE], 2, mu[cmn], "-"), 2, sd[cmn], "/")
}
tr_x <- scale_scores(tr_scores, tr_mean, tr_sd)

set.seed(123)
cv <- cv.glmnet(tr_x, tr_y, family = "cox", alpha = 1, nfolds = 10,
                standardize = FALSE)

# ---- Fig2a 交叉验证曲线（标准 LASSO 样式：顶轴=非零程序数，图例标注两条 lambda 线）----
fit      <- cv$glmnet.fit
coef_mat <- as.matrix(coef(fit))          # (截距 + 9 程序) × length(lambda)
nzero    <- colSums(coef_mat[-1, , drop = FALSE] != 0)   # 每个 lambda 的非零程序数
lx       <- log(cv$lambda)

cvdf <- data.frame(lambda = lx, cvm = cv$cvm,
                   cvlo = cv$cvlo, cvup = cv$cvup, nzero = nzero)

# 顶轴刻度：每个不同程序数取其中位 log(lambda) 位置（左=程序多 → 右=程序少）
uniq_nz <- sort(unique(nzero), decreasing = TRUE)
brk_lx  <- sapply(uniq_nz, function(z) median(lx[nzero == z]))

lab_min <- sprintf("lambda.min = %.3g (%d programs)",
                   cv$lambda.min, sum(coef(cv, s = "lambda.min") != 0))
lab_1se <- sprintf("lambda.1se = %.3g (%d programs)",
                   cv$lambda.1se, sum(coef(cv, s = "lambda.1se") != 0))

p_a <- ggplot(cvdf, aes(lambda, cvm)) +
  geom_point(size = 1.2) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), width = 0, alpha = 0.4) +
  geom_vline(aes(xintercept = log(cv$lambda.min),
                 color = "lambda.min", linetype = "lambda.min")) +
  geom_vline(aes(xintercept = log(cv$lambda.1se),
                 color = "lambda.1se", linetype = "lambda.1se")) +
  scale_color_manual(name = NULL,
                     values = c("lambda.min" = "firebrick",
                                "lambda.1se" = "grey40"),
                     labels = c("lambda.min" = lab_min,
                                "lambda.1se" = lab_1se)) +
  scale_linetype_manual(name = NULL,
                        values = c("lambda.min" = "dashed",
                                   "lambda.1se" = "dotted"),
                        labels = c("lambda.min" = lab_min,
                                   "lambda.1se" = lab_1se)) +
  labs(x = "log(lambda)", y = "Partial likelihood deviance",
       title = "LASSO-Cox cross-validation") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)),
                     sec.axis = sec_axis(~ ., name = "Number of nonzero coefficients",
                                         breaks = brk_lx, labels = uniq_nz)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  theme_classic() +
  theme(legend.position = "top")
ggsave(file.path(fig_dir, "method_2a_cv_curve.png"), p_a, width = 6, height = 4.5, dpi = 300)

# ---- Fig2b 系数路径 ----
coef_long <- data.frame(
  lambda  = rep(log(fit$lambda), each = nrow(coef_mat)),
  program = rep(rownames(coef_mat), times = ncol(coef_mat)),
  beta    = as.vector(coef_mat)
)
p_b <- ggplot(coef_long, aes(lambda, beta, color = program)) +
  geom_line(size = 1) +
  geom_vline(xintercept = log(cv$lambda.min), linetype = "dashed", color = "firebrick") +
  labs(x = "log(lambda)", y = "Coefficient", title = "Coefficient path") +
  theme_classic() +
  theme(legend.title = element_blank())
ggsave(file.path(fig_dir, "method_2b_coef_path.png"), p_b, width = 6.5, height = 4.5, dpi = 300)

# ---- Fig2c 最终系数柱状图 ----
coef_sel <- sig$coef[sig$sel]
coef_df  <- data.frame(program = names(coef_sel), beta = coef_sel,
                       stringsAsFactors = FALSE)
coef_df$program <- factor(coef_df$program, levels = coef_df$program[order(coef_df$beta)])
coef_df$sign    <- ifelse(coef_df$beta > 0, "Adverse", "Protective")

p_c <- ggplot(coef_df, aes(program, beta, fill = sign)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c(Adverse = "#C84630", Protective = "#2E86AB"),
                    guide = "none") +
  labs(x = NULL, y = "LASSO coefficient",
       title = "TAM program signature coefficients") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12))
ggsave(file.path(fig_dir, "method_2c_coef_bar.png"), p_c, width = 6, height = 4, dpi = 300)

cat("三张方法图已保存到:", fig_dir, "\n")
cat("lambda.min =", round(cv$lambda.min, 4),
    " | lambda.1se 选中程序数:", sum(coef(cv, s = "lambda.1se") != 0), "\n")
