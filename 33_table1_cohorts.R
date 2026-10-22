# ============================================================
# 33_table1_cohorts.R
# Table 1：4 队列临床特征汇总（N / 事件 / 平台 / 中位随访 / 年龄 / 性别 / 分期）
#   生存 = bulk_program_scores.rds 的 harmonized surv（time=天, event=0/1）
#   临床 = TCGA-LUAD.rds$clinical（Xena） + GSE*_eset.rds$pdata（series_matrix 元数据，grep 自动匹配列）
# 输出：results/tables/table1_cohorts.csv + 控制台打印（含每队列用到的临床列名，便于核对）
# 依赖：survival
# 运行：本地 RStudio
# ============================================================

suppressPackageStartupMessages({ library(survival) })

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

platform <- c("TCGA-LUAD" = "RNA-seq (STAR)",
              "GSE68465"  = "Affymetrix HG-U133A (GPL96)",
              "GSE72094"  = "Affymetrix Gene 1.0 ST (GPL15048)",
              "GSE31210"  = "Affymetrix HG-U133 Plus 2.0 (GPL570)")

# ---------- 工具 ----------
followup_months <- function(time, event) {           # 反向 KM 中位随访（天 → 月）
  fit <- survfit(Surv(time, 1 - event) ~ 1)
  m <- as.numeric(summary(fit)$table["median"])
  if (is.na(m)) m <- max(time, na.rm = TRUE)
  m / 30.44
}
med_iqr <- function(x) {                             # "65 [58-72]"
  x <- as.numeric(x); x <- x[!is.na(x)]
  if (length(x) < 3) return(NA_character_)
  sprintf("%.0f [%.0f-%.0f]", median(x), quantile(x, 0.25), quantile(x, 0.75))
}
norm_sex <- function(x) {
  x <- trimws(tolower(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[grepl("^(m|male|1)$", x)] <- "male"
  out[grepl("^(f|female|2)$", x)] <- "female"
  out
}
sex_str <- function(x) {                             # "238 (47%)"
  x <- norm_sex(x); nm <- sum(x == "male", na.rm = TRUE); nf <- sum(x == "female", na.rm = TRUE)
  if (nm + nf == 0) return(NA_character_)
  sprintf("%d (%d%%)", nm, round(100 * nm / (nm + nf)))
}
stage_num <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^STAGE[ ]*", "", x)
  x <- gsub("III", "3", x); x <- gsub("II", "2", x)
  x <- gsub("IV", "4", x);  x <- gsub("I", "1", x)
  suppressWarnings(as.integer(sub("([0-9]).*", "\\1", x)))
}
stage_str <- function(x) {                           # "I=180; II=120; III=88; IV=26"
  s <- stage_num(x); s <- s[s %in% 1:4]
  if (length(s) == 0) return(NA_character_)
  tb <- table(factor(s, levels = 1:4))
  paste(sprintf("I=%d; II=%d; III=%d; IV=%d", tb[1], tb[2], tb[3], tb[4]), collapse = "")
}
find_col <- function(pdata, pats) {
  hit <- unique(unlist(lapply(pats, function(p)
    grep(p, colnames(pdata), ignore.case = TRUE, value = TRUE))))
  if (length(hit) == 0) NA_character_ else hit[1]
}

# ---------- 生存（N / 事件 / 随访）----------
d <- readRDS(file.path(val_dir, "bulk_program_scores.rds"))
surv <- d$surv

rows <- list()
for (g in names(platform)) {
  s <- surv[[g]]
  rows[[g]] <- data.frame(
    cohort = g, platform = platform[g],
    N = nrow(s), events = sum(s$event == 1, na.rm = TRUE),
    followup_months = round(followup_months(s$time, s$event), 1),
    stringsAsFactors = FALSE)
}
tab <- do.call(rbind, rows)

# ---------- 临床（年龄/性别/分期）----------
age_c <- sex_c <- stage_c <- rep(NA_character_, nrow(tab))
names(age_c) <- names(sex_c) <- names(stage_c) <- names(platform)
src <- character(0)

fmt_stage <- function(s) {                          # s: 1-4 数值 → "I=..; II=..; III=..; IV=.."
  s <- s[s %in% 1:4]
  if (length(s) == 0) return(NA_character_)
  tb <- table(factor(s, levels = 1:4))
  paste(sprintf("I=%d; II=%d; III=%d; IV=%d", tb[1], tb[2], tb[3], tb[4]), collapse = "")
}
# GSE68465 的 disease_stage 是 TNM（"pN0pT1"），解析成 1-3 期（简化版 AJCC，无 T 亚分期）
stage_from_tnm <- function(x) {
  x  <- as.character(x)
  pT <- suppressWarnings(as.integer(sub(".*pT([1-4]).*", "\\1", x)))
  pN <- suppressWarnings(as.integer(sub(".*pN([0-3]).*", "\\1", x)))
  stg <- rep(NA_integer_, length(x))
  stg[pT %in% 1:2 & pN == 0] <- 1L
  stg[(pT == 3 & pN == 0) | (pT %in% 1:2 & pN == 1)] <- 2L
  stg[(pT == 4 & pN %in% 0:1) | (pT %in% 1:3 & pN == 2) | (pN == 3)] <- 3L
  stg
}

# TCGA（Xena clinical，字段名已知）
tcga <- readRDS(file.path(val_dir, "TCGA-LUAD.rds"))
clin <- tcga$clinical
key_clin <- substr(as.character(clin$sample), 1, 15)
idx <- match(rownames(surv[["TCGA-LUAD"]]), key_clin)
age_c["TCGA-LUAD"]   <- med_iqr(clin[["age_at_earliest_diagnosis_in_years.diagnoses.xena_derived"]][idx])
sex_c["TCGA-LUAD"]   <- sex_str(clin[["gender.demographic"]][idx])
stage_c["TCGA-LUAD"] <- fmt_stage(stage_num(clin[["ajcc_pathologic_stage.diagnoses"]][idx]))
src <- c(src, "TCGA-LUAD: age/sex/stage 用 Xena clinical 已知字段")

# GEO（字段名已确认，显式指定；GSE31210 为早期切除队列，III/IV 被预后排除过滤）
geo_clin <- list(
  GSE68465 = list(age = "age",            sex = "Sex",    stage = "disease_stage",      mode = "tnm"),
  GSE72094 = list(age = "age_at_diagnosis", sex = "gender", stage = "Stage",             mode = "roman"),
  GSE31210 = list(age = "age (years)",   sex = "gender", stage = "pathological stage", mode = "roman")
)
for (g in names(geo_clin)) {
  es <- readRDS(file.path(val_dir, paste0(g, "_eset.rds")))
  pdata <- es$pdata[rownames(surv[[g]]), , drop = FALSE]
  cfg <- geo_clin[[g]]
  if (cfg$age %in% colnames(pdata))  age_c[g] <- med_iqr(pdata[[cfg$age]])
  if (cfg$sex %in% colnames(pdata))  sex_c[g] <- sex_str(pdata[[cfg$sex]])
  if (cfg$stage %in% colnames(pdata)) {
    stg <- if (cfg$mode == "tnm") stage_from_tnm(pdata[[cfg$stage]]) else stage_num(pdata[[cfg$stage]])
    stage_c[g] <- fmt_stage(stg)
  }
  src <- c(src, sprintf("%s: age=%s sex=%s stage=%s(%s)", g, cfg$age, cfg$sex, cfg$stage, cfg$mode))
}

tab$age    <- age_c
tab$sex_male_n_pct <- sex_c
tab$stage  <- stage_c

cat("===== 各队列临床字段来源（请核对 grep 是否抓对列）=====\n")
cat(paste(src, collapse = "\n"), "\n\n")

cat("===== Table 1 =====\n")
print(tab, row.names = FALSE)

write.csv(tab, file.path(tab_dir, "table1_cohorts.csv"), row.names = FALSE)
cat("\n已保存:", file.path(tab_dir, "table1_cohorts.csv"), "\n")
