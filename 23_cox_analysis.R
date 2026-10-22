# ============================================================
# 23_cox_analysis.R
# 1. 单因素 Cox（TCGA）：9 程序各自 + 签名的 HR / P / C-index
# 2. 多因素 Cox（TCGA）：签名 + stage + age + gender（证明独立预后）
# 依赖：survival（依赖 20/21 脚本的输出）
# 运行：本地 RStudio
# ============================================================

suppressPackageStartupMessages({ library(survival) })

base    <- "."
res_dir <- file.path(base, "results")

d    <- readRDS(file.path(base, "validation_cohorts", "bulk_program_scores.rds"))
sig  <- readRDS(file.path(res_dir, "rds", "signature_results.rds"))
tcga <- readRDS(file.path(base, "validation_cohorts", "TCGA-LUAD.rds"))

scores    <- d$scores[["TCGA-LUAD"]]
surv_tcga <- d$surv[["TCGA-LUAD"]]
risk      <- sig$risk[["TCGA-LUAD"]]
common    <- intersect(rownames(scores), names(risk))

# ---------- 工具：单因素 Cox ----------
univ <- function(x, s, label) {
  df <- data.frame(time = s$time, event = s$event, x = x)
  df <- df[complete.cases(df), , drop = FALSE]
  fit <- coxph(Surv(time, event) ~ x, data = df)
  ci  <- summary(fit)$conf.int[1, c(1, 3, 4)]
  data.frame(variable = label, N = nrow(df),
             HR = round(ci[1], 3), lower = round(ci[2], 3), upper = round(ci[3], 3),
             P = signif(summary(fit)$coefficients[1, "Pr(>|z|)"], 3),
             Cindex = round(summary(fit)$concordance[1], 3),
             stringsAsFactors = FALSE)
}

cat("===== 1. 单因素 Cox（TCGA，9 程序 + 签名）=====\n")
uni_tab <- do.call(rbind, c(
  lapply(colnames(scores), function(p)
    univ(scores[common, p], surv_tcga[common, , drop = FALSE], p)),
  list(univ(risk[common], surv_tcga[common, , drop = FALSE], "Signature"))
))
print(uni_tab, row.names = FALSE)

# ---------- 2. 多因素 Cox ----------
# 临床对齐：clinical 的 sample 是完整 barcode，取前 15 位
clin     <- tcga$clinical
key_clin <- substr(as.character(clin$sample), 1, 15)
idx      <- match(common, key_clin)
cat("\n临床对齐：总样本", length(common), " | 临床匹配", sum(!is.na(idx)), "\n")

stage_raw  <- clin$ajcc_pathologic_stage.diagnoses[idx]
age_raw    <- clin$age_at_earliest_diagnosis_in_years.diagnoses.xena_derived[idx]
gender_raw <- clin$gender.demographic[idx]

# stage → early/late
stage_num <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^STAGE[ ]*", "", x)
  x <- gsub("III", "3", x); x <- gsub("II", "2", x)
  x <- gsub("IV", "4", x);  x <- gsub("I", "1", x)
  suppressWarnings(as.integer(sub("([0-9]).*", "\\1", x)))
}
stg <- stage_num(stage_raw)
earlylate <- factor(ifelse(stg %in% c(1, 2), "Early",
                    ifelse(stg %in% c(3, 4), "Late", NA)),
                    levels = c("Early", "Late"))
cat("\n阶段分布（原始 → 早/晚）:\n")
print(table(stage_raw, earlylate, useNA = "ifany"))

age    <- suppressWarnings(as.numeric(age_raw))
gender <- factor(tolower(gender_raw))
cat("\nage 范围:", range(age, na.rm = TRUE),
    " | gender:", paste(names(table(gender)), collapse = "/"), "\n")

mv <- data.frame(time  = surv_tcga[common, "time"],
                 event = surv_tcga[common, "event"],
                 risk  = risk[common],
                 stage = earlylate,
                 age   = age,
                 gender = gender)
mv <- mv[complete.cases(mv), , drop = FALSE]
cat("多因素 Cox 样本数:", nrow(mv), " | 事件:", sum(mv$event), "\n")

fit_mv <- coxph(Surv(time, event) ~ risk + stage + age + gender, data = mv)
cat("\n===== 2. 多因素 Cox（签名 + stage + age + gender）=====\n")
mv_tab <- data.frame(
  variable = rownames(summary(fit_mv)$coefficients),
  HR       = round(summary(fit_mv)$conf.int[, 1], 3),
  lower    = round(summary(fit_mv)$conf.int[, 3], 3),
  upper    = round(summary(fit_mv)$conf.int[, 4], 3),
  P        = signif(summary(fit_mv)$coefficients[, "Pr(>|z|)"], 3),
  stringsAsFactors = FALSE
)
print(mv_tab, row.names = FALSE)

# ---------- 保存 ----------
write.csv(uni_tab, file.path(res_dir, "tables", "univariate_cox.csv"), row.names = FALSE)
write.csv(mv_tab,  file.path(res_dir, "tables", "multivariate_cox.csv"), row.names = FALSE)
saveRDS(list(univariate = uni_tab, multivariate = mv_tab),
        file.path(res_dir, "rds", "cox_analysis.rds"))
cat("\n已保存到:", res_dir, "\n")
