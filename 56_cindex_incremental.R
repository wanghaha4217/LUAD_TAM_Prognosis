# ============================================================
# 56_cindex_incremental.R
# 目的：签名 vs 纯临床模型的 C-index 增量 + 似然比检验（LRT）
#       审稿人几乎必问「签名比 stage/age/sex 多多少」——此脚本正面回答。
#   主分析：4 队列统一 stage + age + sex
#   TCGA 额外：stage + age + sex + smoking（敏感性，pack-years 派生 Never/Ever）
# 输出：results/tables/cindex_incremental.csv
# 运行：本地 RStudio，秒级（每队列 2 个 Cox 模型）
# ============================================================
suppressPackageStartupMessages({ library(survival) })

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

d   <- readRDS(file.path(val_dir, "bulk_program_scores.rds"))
sig <- readRDS(file.path(res_dir, "rds", "signature_results.rds"))
tcga <- readRDS(file.path(val_dir, "TCGA-LUAD.rds"))

# ---------- 工具 ----------
stage_num <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^STAGE[ ]*", "", x)
  x <- gsub("III", "3", x); x <- gsub("II", "2", x)
  x <- gsub("IV", "4", x);  x <- gsub("I", "1", x)
  suppressWarnings(as.integer(sub("([0-9]).*", "\\1", x)))
}
earlylate <- function(s) factor(ifelse(s %in% c(1, 2), "Early",
                                ifelse(s %in% c(3, 4), "Late", NA)),
                                levels = c("Early", "Late"))
stage_from_tnm <- function(x) {
  x <- as.character(x)
  pT <- suppressWarnings(as.integer(sub(".*pT([1-4]).*", "\\1", x)))
  pN <- suppressWarnings(as.integer(sub(".*pN([0-3]).*", "\\1", x)))
  stg <- rep(NA_integer_, length(x))
  stg[pT %in% 1:2 & pN == 0] <- 1L
  stg[(pT == 3 & pN == 0) | (pT %in% 1:2 & pN == 1)] <- 2L
  stg[(pT == 4 & pN %in% 0:1) | (pT %in% 1:3 & pN == 2) | (pN == 3)] <- 3L
  stg
}
norm_sex <- function(x) {
  x <- trimws(tolower(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[grepl("^(m|male|1)$", x)] <- "male"
  out[grepl("^(f|female|2)$", x)] <- "female"
  out
}

# 单队列：临床模型 vs 临床+签名模型
run_cohort <- function(g, time, event, risk, stage, age, sex, smoke = NULL) {
  df <- data.frame(time = time, event = event, risk = risk,
                   stage = stage, age = age, sex = sex)
  if (!is.null(smoke)) df$smoke <- smoke
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 20) return(NULL)

  form_clin <- if (is.null(smoke)) "Surv(time,event) ~ stage + age + sex"
               else "Surv(time,event) ~ stage + age + sex + smoke"
  fit_clin <- coxph(as.formula(form_clin), data = df)
  fit_full <- coxph(as.formula(paste0(form_clin, " + risk")), data = df)

  C_clin <- summary(fit_clin)$concordance[1]
  C_full <- summary(fit_full)$concordance[1]

  # 手动 LRT（嵌套模型，避免 anova 列名版本差异）
  chisq <- 2 * (as.numeric(logLik(fit_full)) - as.numeric(logLik(fit_clin)))
  df_x  <- length(coef(fit_full)) - length(coef(fit_clin))
  lrt_p <- pchisq(chisq, df_x, lower.tail = FALSE)

  data.frame(cohort = g,
             model  = ifelse(is.null(smoke), "stage+age+sex", "stage+age+sex+smoke"),
             N = nrow(df), events = sum(df$event),
             C_clin = round(C_clin, 3), C_full = round(C_full, 3),
             deltaC = round(C_full - C_clin, 3),
             LRT_P  = signif(lrt_p, 3),
             stringsAsFactors = FALSE)
}

out <- list()

# ---------- TCGA ----------
surv_tcga <- d$surv[["TCGA-LUAD"]]
risk_tcga <- sig$risk[["TCGA-LUAD"]]
common <- intersect(rownames(surv_tcga), names(risk_tcga))
clin <- tcga$clinical
key  <- substr(as.character(clin$sample), 1, 15)
idx  <- match(common, key)

stg <- earlylate(stage_num(clin$ajcc_pathologic_stage.diagnoses[idx]))
age <- suppressWarnings(as.numeric(clin$age_at_earliest_diagnosis_in_years.diagnoses.xena_derived[idx]))
sex <- factor(norm_sex(clin$gender.demographic[idx]), levels = c("male", "female"))

py_col <- grep("pack_year", colnames(clin), ignore.case = TRUE, value = TRUE)[1]
smoke  <- NULL
if (!is.na(py_col)) {
  py <- suppressWarnings(as.numeric(clin[[py_col]][idx]))
  smoke <- factor(ifelse(is.na(py) | py <= 0, "Never", "Ever"), levels = c("Never", "Ever"))
}

out[["TCGA-3var"]] <- run_cohort("TCGA-LUAD", surv_tcga[common, "time"], surv_tcga[common, "event"],
                                 risk_tcga[common], stg, age, sex, NULL)
out[["TCGA-4var"]] <- run_cohort("TCGA-LUAD", surv_tcga[common, "time"], surv_tcga[common, "event"],
                                 risk_tcga[common], stg, age, sex, smoke)

# ---------- GEO（3 队列，均有 stage/age/sex）----------
geo_clin <- list(
  GSE68465 = list(age = "age",            sex = "Sex",    stage = "disease_stage",      mode = "tnm"),
  GSE72094 = list(age = "age_at_diagnosis", sex = "gender", stage = "Stage",             mode = "roman"),
  GSE31210 = list(age = "age (years)",   sex = "gender", stage = "pathological stage", mode = "roman")
)
for (g in names(geo_clin)) {
  es <- readRDS(file.path(val_dir, paste0(g, "_eset.rds")))
  pdata <- es$pdata
  s <- d$surv[[g]]; r <- sig$risk[[g]]
  common <- intersect(rownames(s), names(r))
  cfg <- geo_clin[[g]]
  pd  <- pdata[common, , drop = FALSE]

  age <- suppressWarnings(as.numeric(pd[[cfg$age]]))
  sex <- factor(norm_sex(pd[[cfg$sex]]), levels = c("male", "female"))
  stg_raw <- if (cfg$mode == "tnm") stage_from_tnm(pd[[cfg$stage]]) else stage_num(pd[[cfg$stage]])
  stg <- earlylate(stg_raw)

  out[[g]] <- run_cohort(g, s[common, "time"], s[common, "event"], r[common], stg, age, sex, NULL)
}

res <- do.call(rbind, out)
print(res, row.names = FALSE)
write.csv(res, file.path(tab_dir, "cindex_incremental.csv"), row.names = FALSE)
cat("\n解读：deltaC>0 = 加签名后 C-index 提升；LRT_P<0.05 = 签名显著改善模型（independent information）。\n")
cat("已保存:", file.path(tab_dir, "cindex_incremental.csv"), "\n")
