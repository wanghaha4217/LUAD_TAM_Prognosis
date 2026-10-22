# ============================================================
# 22_km_curves.R
# 目的：4 个队列的 KM 生存曲线（高/低风险组），标注 log-rank P 和 HR(95%CI)
# 依赖：survival / survminer（若报"没有 survminer"，install.packages("survminer")）
# 运行：本地 RStudio（依赖 21 脚本输出的 signature_results.rds）
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
})

base    <- "."
res_dir <- file.path(base, "results")
fig_dir <- file.path(res_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

sig <- readRDS(file.path(res_dir, "rds", "signature_results.rds"))
d   <- readRDS(file.path(base, "validation_cohorts", "bulk_program_scores.rds"))
surv <- d$surv
risk <- sig$risk            # list: cohort -> 命名风险分向量

cohorts <- names(risk)
plot_list <- list()

for (g in cohorts) {
  r <- risk[[g]]
  s <- surv[[g]]
  common <- intersect(names(r), rownames(s))
  df <- data.frame(
    time  = s[common, "time"] / 365.25,          # 天 → 年
    event = s[common, "event"],
    risk  = r[common]
  )
  df <- df[!is.na(df$time) & !is.na(df$event) & !is.na(df$risk), ]

  med <- median(df$risk)
  df$group <- factor(ifelse(df$risk > med, "High", "Low"), levels = c("Low", "High"))

  fit  <- survfit(Surv(time, event) ~ group, data = df)
  lr_p <- survdiff(Surv(time, event) ~ group, data = df)$pvalue
  cf   <- summary(coxph(Surv(time, event) ~ group, data = df))$conf.int[1, c(1, 3, 4)]
  hr_txt <- sprintf("HR = %.2f (%.2f-%.2f)", cf[1], cf[2], cf[3])

  p <- ggsurvplot(fit, data = df,
                  pval = FALSE,
                  risk.table = TRUE,
                  palette = c("#2E86AB", "#C84630"),
                  xlab = "Time (years)", ylab = "Overall survival",
                  legend.labs = c("Low risk", "High risk"),
                  legend.title = "Group",
                  title = g,
                  ggtheme = theme_classic2(),
                  font.legend = 10,
                  break.time.by = 2)

  # 手动加 log-rank P 和 HR
  p$plot <- p$plot +
    annotate("text", x = 0, y = 0.10, hjust = 0,
             label = sprintf("log-rank P = %.3g\n%s", lr_p, hr_txt), size = 3.5)

  # 单张图：只存 p$plot（真 ggplot，含曲线+注释）。risk table 放合并图里，单张不放，
  # 因为 ggsave 直接吃 ggsurvplot 列表会存成黑图
  ggsave(file.path(fig_dir, paste0("KM_", g, ".png")),
         plot = p$plot, width = 6.5, height = 5, dpi = 300)
  plot_list[[g]] <- p
  cat(g, ": 完成\n", sep = "")
}

# ---- 合并成 2×2 一张图（论文用） ----
if (length(plot_list) == 4) {
  combo <- arrange_ggsurvplots(plot_list, ncol = 2, nrow = 2)
  ggsave(file.path(fig_dir, "KM_combined.png"),
         plot = combo, width = 12, height = 10, dpi = 300)
  cat("合并图已保存\n")
}

cat("\n全部完成。输出目录:", fig_dir, "\n")
