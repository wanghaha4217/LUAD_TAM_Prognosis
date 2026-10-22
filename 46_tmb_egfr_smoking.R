# ============================================================
# 46_tmb_egfr_smoking.R
# 目的：TMB + EGFR 突变 + 吸烟状态 与 7 程序签名风险 的关系（TCGA-LUAD）
#   连接逻辑链：免疫冷(脚本38/45) → 低 TMB → 抗 PD-1 不敏感(脚本39 阴性)
#   并做 EGFR 驱动亚型 / 吸烟史 的风险分层
# 输入：
#   TCGA-LUAD_MAF/data_mutations.txt（cBioPortal 下载，用户已下好）
#   validation_cohorts/TCGA-LUAD.rds（临床 + counts）
#   validation_cohorts/bulk_program_scores.rds（生存 time/event）
#   results/rds/signature_results.rds（风险分）
# 输出：
#   tables/tmb_egfr_smoking_by_risk.csv（组比较汇总）
#   tables/tmb_risk_spearman.csv（TMB 与风险相关）
#   figures/tmb_km.png / egfr_km.png / smoking_km.png / risk_by_egfr_smoking.png
# 依赖：data.table（读 232MB MAF 快；缺则 read.delim 慢读）; survival; ggplot2
# 运行：本地 RStudio
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
})

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures")
maf_dir <- file.path(base, "TCGA-LUAD_MAF")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- MAF 文件定位 ----------
maf_file <- list.files(maf_dir, pattern = "data_mutations", full.names = TRUE)
if (length(maf_file) == 0) stop("找不到 MAF，请把 data_mutations.txt 放到 ", maf_dir)
maf_path <- maf_file[1]
cat("MAF 文件:", maf_path, "\n")

has_dt <- requireNamespace("data.table", quietly = TRUE)
if (!has_dt) cat("【缺包】data.table -> install.packages('data.table')（可选，快读大文件）\n")

# 只读需要的列（省内存）
con <- file(maf_path, "r"); hdr <- readLines(con, n = 1); close(con)
all_cols <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
need <- c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")
have_hgvsp <- "HGVSp_Short" %in% all_cols
sel <- intersect(c(need, if (have_hgvsp) "HGVSp_Short" else "HGVSp"), all_cols)
cat("读取 MAF 列:", paste(sel, collapse = ", "), "\n")

if (has_dt) {
  maf <- as.data.frame(data.table::fread(maf_path, select = sel, sep = "\t",
                                         header = TRUE, quote = "", showProgress = TRUE))
} else {
  maf <- read.delim(maf_path, sep = "\t", header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "")
  maf <- maf[, intersect(sel, colnames(maf)), drop = FALSE]
}
cat("MAF 行数:", nrow(maf), "\n")

if (!have_hgvsp && "HGVSp" %in% colnames(maf)) maf$HGVSp_Short <- maf$HGVSp
if (!"HGVSp_Short" %in% colnames(maf)) maf$HGVSp_Short <- NA_character_

# 非沉默变异（TMB / EGFR 都用这个口径）
nonsyn_cls <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins",
                "In_Frame_Del","In_Frame_Ins","Splice_Site","Translation_Start_Site",
                "Nonstop_Mutation","Start_Codon_Del","Start_Codon_Ins","Stop_Codon_Del","Stop_Codon_Ins")
maf$is_nonsyn <- maf$Variant_Classification %in% nonsyn_cls
maf$sid       <- substr(maf$Tumor_Sample_Barcode, 1, 15)

# ---------- 风险 + 生存 ----------
d    <- readRDS(file.path(val_dir, "bulk_program_scores.rds"))
sig  <- readRDS(file.path(res_dir, "rds", "signature_results.rds"))
tcga <- readRDS(file.path(val_dir, "TCGA-LUAD.rds"))

surv   <- d$surv[["TCGA-LUAD"]]
risk   <- sig$risk[["TCGA-LUAD"]]
common <- intersect(rownames(surv), names(risk))
time   <- surv[common, "time"]; event <- surv[common, "event"]; rk <- risk[common]
cat("风险+生存样本:", length(common), " | 事件:", sum(event), "\n")

# ---------- TMB ----------
nsc <- table(maf$sid[maf$is_nonsyn])
tmb_cnt <- setNames(rep(0L, length(common)), common)
mm <- intersect(names(nsc), names(tmb_cnt))   # 只匹配 common 里的样本，避免 MAF 多出的样本拉长向量
tmb_cnt[mm] <- as.integer(nsc[mm])
tmb <- tmb_cnt / 38   # mut/Mb，WES 外显子按 38Mb（Samstein 口径；只影响绝对值，不影响分层）
cat("TMB 中位:", round(median(tmb), 2), "mut/Mb | 计数中位:", median(tmb_cnt),
    " | >=10 mut/Mb:", sum(tmb >= 10), "\n")

sp <- cor.test(tmb, rk, method = "spearman")
cat("\n===== TMB vs 风险（Spearman）=====\n")
cat("rho =", round(sp$estimate, 3), "| P =", signif(sp$p.value, 3), "\n")
write.csv(data.frame(rho = sp$estimate, P = sp$p.value),
          file.path(tab_dir, "tmb_risk_spearman.csv"), row.names = FALSE)

# ---------- EGFR 突变 ----------
egfr_rows <- maf[maf$Hugo_Symbol == "EGFR" & maf$is_nonsyn, , drop = FALSE]
egfr_sid  <- unique(egfr_rows$sid)
egfr_mut  <- common %in% egfr_sid
cat("\nEGFR 突变样本:", sum(egfr_mut), "/", length(common), "\n")
if (nrow(egfr_rows)) {
  cat("EGFR 突变类型分布（HGVSp_Short）:\n")
  print(sort(table(egfr_rows$HGVSp_Short), decreasing = TRUE))
}

# ---------- 吸烟 ----------
clin <- tcga$clinical
key_clin <- substr(as.character(clin$sample), 1, 15)
idx <- match(common, key_clin)

smoke_cols <- grep("smok|tobacco|cigarette|pack", colnames(clin), ignore.case = TRUE, value = TRUE)
cat("\n吸烟相关字段:", paste(smoke_cols, collapse = ", "), "\n")
smoke <- NULL

# 1) 优先分类字段（值含 smoker / non-smoker / never）
cat_col <- NULL
for (sc in smoke_cols) {
  v <- as.character(clin[[sc]][idx])
  if (any(grepl("smoker|non.?smoker|never", v, ignore.case = TRUE), na.rm = TRUE)) {
    cat_col <- sc; break
  }
}
if (!is.null(cat_col)) {
  raw <- as.character(clin[[cat_col]][idx])
  cat("吸烟分类字段:", cat_col, "\n"); print(table(raw, useNA = "ifany"))
  smoke <- factor(ifelse(grepl("never|non.?smoker", raw, ignore.case = TRUE), "Never",
                  ifelse(grepl("smoker|smoking", raw, ignore.case = TRUE), "Ever", NA)),
                  levels = c("Never", "Ever"))
} else {
  # 2) 无分类字段 → 用 pack_years 派生：NA/0=Never，>0=Ever（TCGA 惯例，缺失多为从不吸烟）
  py_col <- grep("pack_year", smoke_cols, ignore.case = TRUE, value = TRUE)
  if (length(py_col)) {
    py <- suppressWarnings(as.numeric(clin[[py_col[1]]][idx]))
    cat("吸烟字段(派生自):", py_col[1], "→ Never=pack-years NA/0，Ever=pack-years>0\n")
    cat("  pack-years: NA/0 =", sum(is.na(py) | py <= 0), " | >0 =", sum(py > 0, na.rm = TRUE), "\n")
    smoke <- factor(ifelse(is.na(py) | py <= 0, "Never", "Ever"), levels = c("Never", "Ever"))
  } else {
    cat("【警告】无 pack_years 字段，吸烟分层跳过。\n")
  }
}
if (!is.null(smoke)) cat("吸烟 Never/Ever:", sum(smoke == "Never"), "/", sum(smoke == "Ever"), "\n")

# ---------- 组比较工具 ----------
compare2 <- function(g, rk, time, event, name) {
  keep <- !is.na(g)
  g <- g[keep]; rk <- rk[keep]; time <- time[keep]; event <- event[keep]
  lv <- levels(g)
  if (length(unique(g)) < 2) return(NULL)
  w  <- wilcox.test(rk ~ g)$p.value
  lr <- survdiff(Surv(time, event) ~ g)$pvalue
  hr <- exp(coef(coxph(Surv(time, event) ~ g)))
  data.frame(analysis = name,
             group = paste0(lv[2], " vs ", lv[1], "(ref)"),
             n_ref = sum(g == lv[1]), n_grp = sum(g == lv[2]),
             med_risk_ref = round(median(rk[g == lv[1]]), 4),
             med_risk_grp = round(median(rk[g == lv[2]]), 4),
             Wilcox_P = signif(w, 3),
             logrank_P = signif(lr, 3),
             HR_grp_vs_ref = round(hr, 3),
             stringsAsFactors = FALSE)
}

tmb_grp  <- factor(ifelse(tmb > median(tmb), "TMB-High", "TMB-Low"),
                   levels = c("TMB-Low", "TMB-High"))
egfr_grp <- factor(ifelse(egfr_mut, "EGFR-mut", "EGFR-WT"),
                   levels = c("EGFR-WT", "EGFR-mut"))

res <- rbind(
  compare2(tmb_grp,  rk, time, event, "TMB median-split"),
  compare2(egfr_grp, rk, time, event, "EGFR mutation")
)

if (!is.null(smoke)) {
  res <- rbind(res, compare2(smoke, rk, time, event, "Smoking"))
  ok_smoke <- !is.na(smoke)
  driven <- factor(ifelse(egfr_mut[ok_smoke] & smoke[ok_smoke] == "Never",
                          "EGFR-mut+Never", "Other"), levels = c("Other", "EGFR-mut+Never"))
  res <- rbind(res, compare2(driven, rk[ok_smoke], time[ok_smoke], event[ok_smoke],
                             "EGFR-mut & never-smoker"))
}

cat("\n===== 组比较汇总 =====\n")
print(res, row.names = FALSE)
write.csv(res, file.path(tab_dir, "tmb_egfr_smoking_by_risk.csv"), row.names = FALSE)

# ---------- KM 图 ----------
km_plot <- function(g, rk, time, event, fname, ttl, cols) {
  keep <- !is.na(g)
  g <- g[keep]; time <- time[keep]; event <- event[keep]
  if (length(unique(g)) < 2) { cat("  [跳过KM]", ttl, "组数<2\n"); return(invisible(NULL)) }
  fit <- survfit(Surv(time, event) ~ g)
  png(file.path(fig_dir, fname), width = 1400, height = 1100, res = 200)
  plot(fit, col = cols, lwd = 2.5, mark.time = TRUE,
       xlab = "Overall survival (days)", ylab = "Survival probability", main = ttl)
  legend("topright", legend = levels(g), col = cols, lwd = 2.5)
  dev.off()
  cat("已保存 figures/", fname, "\n", sep = "")
}
km_plot(tmb_grp,  rk, time, event, "tmb_km.png",  "TCGA-LUAD: TMB median-split",  c("#2E86AB", "#C0392B"))
km_plot(egfr_grp, rk, time, event, "egfr_km.png", "TCGA-LUAD: EGFR mutation",    c("#2E86AB", "#C0392B"))
if (!is.null(smoke)) km_plot(smoke, rk, time, event, "smoking_km.png", "TCGA-LUAD: Smoking", c("#7B8B4F", "#B8860B"))

# ---------- 箱线图：风险 vs EGFR / 吸烟 ----------
panel1 <- data.frame(risk = rk, group = egfr_grp, var = "EGFR")
dfb <- panel1
if (!is.null(smoke)) dfb <- rbind(dfb, data.frame(risk = rk, group = smoke, var = "Smoking"))
dfb$var <- factor(dfb$var, levels = c("EGFR", "Smoking"))
p <- ggplot(dfb, aes(x = group, y = risk, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 0.5, alpha = 0.35, color = "grey30") +
  facet_wrap(~ var, scales = "free_x") +
  scale_fill_manual(values = c("#2E86AB", "#C0392B", "#7B8B4F", "#B8860B"), name = NULL) +
  labs(x = NULL, y = "Signature risk score") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", color = NA))
ggsave(file.path(fig_dir, "risk_by_egfr_smoking.png"), p, width = 8, height = 4.5, dpi = 300)

cat("\n已保存 tables/tmb_egfr_smoking_by_risk.csv + tmb_risk_spearman.csv + 4 张图\n")
cat("说明：TMB=非沉默突变数/38Mb；EGFR-mut=任意非沉默 EGFR 突变；吸烟 Never=never/non-smoker，其余=Ever。\n")
cat("     HR>1=该组风险更高/预后更差；方向如实读，勿预设。\n")
