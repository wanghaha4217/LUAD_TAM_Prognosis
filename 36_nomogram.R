# ============================================================
# 36_nomogram.R
# 目的：7 程序签名 + 临床（stage/age/gender）→ Nomogram，
#       并补三件套：校准曲线 / 决策曲线(DCA) / 时间依赖 ROC-AUC
#   建模队列：TCGA-LUAD（唯一有完整临床字段的队列）
# 依赖：rms（nomogram+校准）、timeROC（时间依赖 AUC），缺则跳过并提示
#       DCA 用自写净获益函数（3 年时点二值化近似，免第三方包）
# 运行：本地 RStudio；首次需 install.packages(c("rms","timeROC"))
# ============================================================

suppressPackageStartupMessages({ library(survival) })

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

has_rms  <- requireNamespace("rms",     quietly = TRUE)
has_troc <- requireNamespace("timeROC", quietly = TRUE)
if (!has_rms)  cat("【缺包】rms     -> install.packages('rms')（nomogram/校准）\n")
if (!has_troc) cat("【缺包】timeROC -> install.packages('timeROC')（时间依赖 ROC）\n")

# ---------- 数据 ----------
d    <- readRDS(file.path(val_dir, "bulk_program_scores.rds"))
sig  <- readRDS(file.path(res_dir, "rds", "signature_results.rds"))
tcga <- readRDS(file.path(val_dir, "TCGA-LUAD.rds"))

surv   <- d$surv[["TCGA-LUAD"]]
risk   <- sig$risk[["TCGA-LUAD"]]
common <- intersect(rownames(surv), names(risk))

clin     <- tcga$clinical
key_clin <- substr(as.character(clin$sample), 1, 15)
idx      <- match(common, key_clin)

stage_num <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^STAGE[ ]*", "", x)
  x <- gsub("III", "3", x); x <- gsub("II", "2", x)
  x <- gsub("IV", "4", x);  x <- gsub("I", "1", x)
  suppressWarnings(as.integer(sub("([0-9]).*", "\\1", x)))
}
stg   <- stage_num(clin$ajcc_pathologic_stage.diagnoses[idx])
stage <- factor(ifelse(stg %in% c(1, 2), "Early",
                ifelse(stg %in% c(3, 4), "Late", NA)),
                levels = c("Early", "Late"))
age    <- suppressWarnings(as.numeric(clin$age_at_earliest_diagnosis_in_years.diagnoses.xena_derived[idx]))
gender <- factor(tolower(clin$gender.demographic[idx]))

df <- data.frame(time  = surv[common, "time"],
                 event = surv[common, "event"],
                 risk  = risk[common],
                 stage = stage, age = age, gender = gender)
df <- df[complete.cases(df), , drop = FALSE]
cat("建模样本:", nrow(df), " | 事件:", sum(df$event), "\n")

# ---------- 1. 多因素 Cox（与脚本 23 同口径，供核对）----------
fit_cox <- coxph(Surv(time, event) ~ risk + stage + age + gender, data = df)
mv <- data.frame(
  variable = rownames(summary(fit_cox)$coefficients),
  HR = round(summary(fit_cox)$conf.int[, 1], 3),
  lo = round(summary(fit_cox)$conf.int[, 3], 3),
  hi = round(summary(fit_cox)$conf.int[, 4], 3),
  P  = signif(summary(fit_cox)$coefficients[, "Pr(>|z|)"], 3))
cat("\n===== 多因素 Cox（nomogram 变量）=====\n")
print(mv, row.names = FALSE)
write.csv(mv, file.path(tab_dir, "nomogram_multivariate_cox.csv"), row.names = FALSE)

# ---------- 2. Nomogram + 校准（rms）----------
if (has_rms) {
  suppressPackageStartupMessages(library(rms))
  dd <- datadist(df); options(datadist = "dd")
  f <- cph(Surv(time, event) ~ risk + stage + age + gender, data = df,
           x = TRUE, y = TRUE, surv = TRUE, time.inc = 365)
  survfun <- Survival(f)

  nom <- nomogram(f,
                  fun = list(function(x) survfun(365,  x),
                             function(x) survfun(1095, x),
                             function(x) survfun(1825, x)),
                  funlabel = c("1-yr OS", "3-yr OS", "5-yr OS"),
                  fun.at = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  png(file.path(fig_dir, "nomogram.png"), width = 2000, height = 1300, res = 200)
  plot(nom)
  dev.off()
  cat("已保存 figures/nomogram.png\n")

  # 校准（1/3/5 年各一张，KM 校准）
  for (u in c(365, 1095, 1825)) {
    cal <- calibrate(f, cmethod = "KM", u = u, m = 100, B = 200)
    png(file.path(fig_dir, sprintf("calibration_%dy.png", round(u / 365))),
        width = 1200, height = 1000, res = 200)
    plot(cal, xlab = "Predicted survival", ylab = "Observed survival", subtitles = FALSE)
    title(main = paste(round(u / 365), "-year calibration"))
    dev.off()
  }
  cat("已保存 figures/calibration_{1,3,5}y.png\n")
}

# ---------- 3. 时间依赖 ROC-AUC（timeROC）----------
if (has_troc) {
  suppressPackageStartupMessages(library(timeROC))
  troc <- timeROC(T = df$time, delta = df$event, marker = df$risk,
                  cause = 1, weighting = "marginal",
                  times = c(365, 1095, 1825), iid = TRUE)
  auc_tab <- data.frame(
    horizon = c("1-year", "3-year", "5-year"),
    AUC = round(as.numeric(troc$AUC), 3))
  # confint 返回的是百分比（×100），除 100 换成与 AUC 同单位的 0-1 小数
  ci_tab <- tryCatch(round(confint(troc, level = 0.95)$CI_AUC / 100, 3),
                     error = function(e) NULL)
  if (!is.null(ci_tab) && ncol(ci_tab) >= 2) {
    auc_tab$lo <- ci_tab[, 1]; auc_tab$hi <- ci_tab[, 2]
  }
  cat("\n===== 时间依赖 ROC-AUC（签名风险分）=====\n")
  print(auc_tab, row.names = FALSE)
  write.csv(auc_tab, file.path(tab_dir, "timeROC_auc.csv"), row.names = FALSE)

  for (tt in c(365, 1095, 1825)) {
    png(file.path(fig_dir, sprintf("timeROC_%dy.png", round(tt / 365))),
        width = 1000, height = 1000, res = 200)
    plot(troc, time = tt, col = "#E41A1C", lwd = 2, title = FALSE)
    dev.off()
  }
  cat("已保存 figures/timeROC_{1,3,5}y.png\n")
}

# ---------- 4. DCA（自写净获益，3 年时点）----------
# 二值化近似：事件发生在 3 年内=阳性，3 年内仍存活且未截尾=阴性，
# 3 年内先截尾者结局未知→剔除。审稿若挑截尾处理，可换 dcurves 包做 KM 加权版。
dca_surv <- function(time, event, risk, horizon = 1095, thresholds = seq(0, 0.5, 0.02)) {
  pos  <- event == 1 & time <= horizon
  neg  <- event == 0 & time >  horizon
  keep <- pos | neg
  y    <- as.integer(pos[keep]); p <- risk[keep]; N <- length(y)
  prev <- mean(y)
  nb_model <- sapply(thresholds, function(pt) {
    w  <- pt / (1 - pt)
    pp <- p >= pt
    (sum(pp & y == 1) - w * sum(pp & y == 0)) / N
  })
  nb_all <- prev - (1 - prev) * thresholds / (1 - thresholds)
  data.frame(threshold = thresholds, model = nb_model, all = nb_all, none = 0)
}
dca <- dca_surv(df$time, df$event, df$risk, horizon = 1095)
write.csv(dca, file.path(tab_dir, "dca_net_benefit.csv"), row.names = FALSE)

png(file.path(fig_dir, "dca.png"), width = 1200, height = 900, res = 200)
yl <- range(c(dca$model, dca$all), na.rm = TRUE)
plot(dca$threshold, dca$model, type = "l", col = "#E41A1C", lwd = 2.5,
     xlab = "Threshold probability", ylab = "Net benefit", ylim = yl)
lines(dca$threshold, dca$all,  col = "grey40", lwd = 2, lty = 2)
lines(dca$threshold, dca$none, col = "grey70", lwd = 2, lty = 3)
abline(h = 0, col = "grey85", lty = 1)
legend("topright", legend = c("Signature", "Treat all", "Treat none"),
       col = c("#E41A1C", "grey40", "grey70"), lty = c(1, 2, 3), lwd = 2)
dev.off()
cat("已保存 figures/dca.png\n")

cat("\n===== 完成 =====\n")
