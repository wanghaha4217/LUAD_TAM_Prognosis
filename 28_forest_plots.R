# ============================================================
# 28_forest_plots.R  (v5 — 加显著性星号 + 加粗线 + 深色无效线)
# forestplot 多列布局：左=变量/队列(N)，中=HR点+CI，右=HR(95%CI)+P值(带星号)
#   统一：同色、同坐标、同列结构；标题放图注
#   forest_multivariate.png  多因素 Cox（签名 + stage + age + gender）
#   forest_validation.png    4 队列 HR(高/低)
# 依赖：survival / forestplot —— 先 install.packages("forestplot")
# ============================================================

suppressPackageStartupMessages({
  library(survival)
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

# ---- 多因素 Cox ----
mv <- read.csv(file.path(res_dir, "tables", "multivariate_cox.csv"), stringsAsFactors = FALSE)
relabel <- c(risk = "Signature (risk)", stageLate = "Stage (late vs early)",
             age = "Age (per year)", gendermale = "Gender (male vs female)")
mv$variable <- relabel[mv$variable]
ord_mv <- c("Signature (risk)", "Stage (late vs early)", "Age (per year)", "Gender (male vs female)")
mv <- mv[order(match(mv$variable, ord_mv)), ]

labeltext_mv <- list(
  c("Variable", mv$variable),
  c("HR (95% CI)", sprintf("%.2f (%.2f-%.2f)", mv$HR, mv$lower, mv$upper)),
  c("P value", fmt_p(mv$P))
)

png(file.path(fig_dir, "forest_multivariate.png"), width = 8, height = 4, units = "in", res = 300)
forestplot(
  labeltext = labeltext_mv,
  mean  = c(NA, mv$HR),
  lower = c(NA, mv$lower),
  upper = c(NA, mv$upper),
  is.summary = c(TRUE, rep(FALSE, nrow(mv))),
  zero   = 1,
  xlog   = TRUE,
  xticks = XTICKS,
  clip   = CLIP,
  xlab   = XLAB,
  col    = fpColors(box = BOX, zero = "grey20"),
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

# ---- 4 队列 HR(高/低) ----
# 切分点 = 各队列内部中位数（非最优阈值 → 无阈值选择乐观偏差）
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
# 固定队列顺序（与 meta_forest 一致）：TCGA(训练) 在最上

labeltext_val <- list(
  c("Cohort", sprintf("%s (N=%d)", fr$cohort, fr$N)),
  c("HR (95% CI)", sprintf("%.2f (%.2f-%.2f)", fr$HR, fr$lower, fr$upper)),
  c("P value", fmt_p(fr$P))
)

png(file.path(fig_dir, "forest_validation.png"), width = 8, height = 4.5, units = "in", res = 300)
forestplot(
  labeltext = labeltext_val,
  mean  = c(NA, fr$HR),
  lower = c(NA, fr$lower),
  upper = c(NA, fr$upper),
  is.summary = c(TRUE, rep(FALSE, nrow(fr))),
  zero   = 1,
  xlog   = TRUE,
  xticks = XTICKS,
  clip   = CLIP,
  xlab   = XLAB,
  col    = fpColors(box = BOX, zero = "grey20"),
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

cat("森林图已保存到:", fig_dir, "\n")
cat("\n===== 跨队列 HR（含 P 值）=====\n"); print(fr)
cat("\n===== 多因素 HR（含 P 值）=====\n"); print(mv)
