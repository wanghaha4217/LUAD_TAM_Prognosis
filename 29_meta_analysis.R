# ============================================================
# 29_meta_analysis.R  (v3 — 双菱形(固定+随机) + 星号 + 加粗无效线)
# 随机效应 Meta 分析：4 队列 HR(高/低) 合并 + 异质性 I²
#   统计用 meta::metagen()，画图用 forestplot（与 28 统一）
#   meta_forest.png   森林图（4 队列点 + common/random 两个菱形 + 右侧 P 值列带星号）
# 依赖：survival / meta / forestplot —— 先 install.packages(c("meta","forestplot"))
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(meta)
  library(forestplot)
})

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
fig_dir <- file.path(res_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

star  <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
fmt_p <- function(p) {
  base <- ifelse(p < 0.001, "P<0.001", sprintf("P=%.3f", p))
  s <- star(p)
  ifelse(nzchar(s), paste0(base, " ", s), base)
}

BOX    <- "#2E86AB"
XLAB   <- "Hazard ratio (95% CI, log scale)"
XTICKS <- c(0.5, 1, 2, 4, 8)
CLIP   <- c(0.5, 8)

# ---- 每队列 HR/CI/P（中位数切分，Wald）----
d   <- readRDS(file.path(val_dir, "bulk_program_scores.rds"))
sig <- readRDS(file.path(res_dir, "rds", "signature_results.rds"))
risk <- sig$risk
surv <- d$surv

fr <- do.call(rbind, lapply(names(risk), function(g) {
  r <- risk[[g]]; s <- surv[[g]]
  common <- intersect(names(r), rownames(s))
  df <- data.frame(time = s[common, "time"], event = s[common, "event"], risk = r[common])
  df <- df[complete.cases(df), , drop = FALSE]
  med <- median(df$risk)
  df$grp <- factor(ifelse(df$risk > med, "High", "Low"), levels = c("Low", "High"))
  fit <- coxph(Surv(time, event) ~ grp, data = df)
  ci  <- summary(fit)$conf.int[1, c(1, 3, 4)]
  wp  <- summary(fit)$coefficients[1, "Pr(>|z|)"]
  data.frame(cohort = g, N = nrow(df), HR = ci[1], lower = ci[2], upper = ci[3],
             P = wp, stringsAsFactors = FALSE)
}))

# ---- Meta（统计计算，不画图）----
m <- metagen(TE    = log(fr$HR),
             seTE  = (log(fr$upper) - log(fr$lower)) / (2 * qnorm(0.975)),
             studlab = fr$cohort,
             sm    = "HR")

# pooled 效应（兼容 meta 7.x 的 "common" 与 6.x 的 "fixed"）
if ("TE.common" %in% names(m)) {
  c_HR <- exp(m$TE.common); c_lo <- exp(m$lower.common); c_up <- exp(m$upper.common); c_P <- m$pval.common
} else if ("TE.fixed" %in% names(m)) {
  c_HR <- exp(m$TE.fixed); c_lo <- exp(m$lower.fixed); c_up <- exp(m$upper.fixed); c_P <- m$pval.fixed
} else {
  c_HR <- c_lo <- c_up <- c_P <- NA
}
r_HR <- exp(m$TE.random); r_lo <- exp(m$lower.random); r_up <- exp(m$upper.random); r_P <- m$pval.random
I2_pct <- ifelse(m$I2 > 1, m$I2, 100 * m$I2)
title_meta <- sprintf("Random-effects meta-analysis: pooled HR %.2f (%.2f–%.2f), I² = %.1f%%, Cochran's Q P = %.3f",
                      r_HR, r_lo, r_up, I2_pct, signif(m$pval.Q, 3))

# ---- 行数据：4 队列 + 2 合并（common + random，画菱形）----
cohort_lab <- c(sprintf("%s (N=%d)", fr$cohort, fr$N), "Common effect", "Random effects")
hr_val <- c(fr$HR, c_HR, r_HR)
lo_val <- c(fr$lower, c_lo, r_lo)
up_val <- c(fr$upper, c_up, r_up)
p_val  <- c(fmt_p(fr$P), fmt_p(c_P), fmt_p(r_P))

labeltext_meta <- list(
  c("Cohort", cohort_lab),
  c("HR (95% CI)", sprintf("%.2f (%.2f-%.2f)", hr_val, lo_val, up_val)),
  c("P value", p_val)
)

png(file.path(fig_dir, "meta_forest.png"), width = 8, height = 5.5, units = "in", res = 300)
forestplot(
  labeltext = labeltext_meta,
  title    = title_meta,
  mean  = c(NA, hr_val),
  lower = c(NA, lo_val),
  upper = c(NA, up_val),
  is.summary = c(TRUE, rep(FALSE, nrow(fr)), TRUE, TRUE),   # 表头 + 两个 pooled = 菱形
  zero   = 1,
  xlog   = TRUE,
  xticks = XTICKS,
  clip   = CLIP,
  xlab   = XLAB,
  col    = fpColors(box = BOX, zero = "grey20", summary = "#1F3864"),
  boxsize = 0.25,
  lwd.ci = 1.5,
  lwd.zero = 1.5,
  lineheight = unit(1.6, "lines"),
  graph.pos = 2,
  txt_gp = fpTxtGp(label = gpar(cex = 0.9),
                   ticks = gpar(cex = 0.85),
                   xlab  = gpar(cex = 0.9))
)
dev.off()

# ---- 关键数字（写进图注）----
cat("图已保存到:", file.path(fig_dir, "meta_forest.png"), "\n")
cat("\n===== 关键数字（写进图注）=====\n")
cat("Common effect HR =", round(c_HR, 2), " [", round(c_lo, 2), "-", round(c_up, 2), "] P =", signif(c_P, 3), "\n")
cat("Random effects HR =", round(r_HR, 2), " [", round(r_lo, 2), "-", round(r_up, 2), "] P =", signif(r_P, 3), "\n")
cat("I2 =", I2_pct, "% | Q test P =", signif(m$pval.Q, 3), " | tau2 =", round(m$tau2, 4), "\n")
