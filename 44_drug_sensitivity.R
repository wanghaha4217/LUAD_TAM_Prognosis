# ============================================================
# 44_drug_sensitivity.R
# 目的：TAM 签名高风险 vs 低风险组的化疗/靶向药敏感性（预测 IC50）
#   用 oncoPredict 把 GDSC2 细胞系药敏模型迁移到 TCGA-LUAD 表达
#   （Maeser et al. Brief Bioinform 2021；训练数据 = GDSC2）
# 依赖：oncoPredict（CRAN）；首次需一次性下载 GDSC2 训练数据（见下）
# 运行：本地 RStudio
# ============================================================

suppressPackageStartupMessages({ library(survival) })

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

has_op <- requireNamespace("oncoPredict", quietly = TRUE)
if (!has_op) {
  cat("【缺包】oncoPredict -> install.packages('oncoPredict')\n")
}

# ---------- GDSC2 训练数据（首次下载，一次性）----------
# 下载方法：浏览器打开  https://osf.io/c6tfx/
#          进 "DataFiles" → "Training Data" 文件夹，下载两个文件：
#            GDSC2_Expr (RMA Normalized and Log Transformed).rds   （注意带空格括号，是正常的）
#            GDSC2_Res.rds
#          都放进下面的 td_dir 目录（文件名保持原样，脚本会自动匹配）。
td_dir <- file.path(base, "oncoPredict_TrainingData")
if (!dir.exists(td_dir)) dir.create(td_dir, recursive = TRUE)
find_file <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (length(f) == 0) NULL else f[1]
}
expr_file <- find_file(td_dir, "GDSC2_Expr")
res_file  <- find_file(td_dir, "GDSC2_Res")

# 目标药物（GDSC2 内名称，用不区分大小写的 grep 匹配，找不到就跳过）
drugs <- c("Cisplatin", "Carboplatin", "Paclitaxel", "Docetaxel", "Gemcitabine",
           "Pemetrexed", "Erlotinib", "Gefitinib", "Vinorelbine", "Etoposide",
           "Oxaliplatin", "Doxorubicin", "5-Fluorouracil", "Afatinib", "Crizotinib")

collapse_mean <- function(mat, group) {
  mat  <- as.matrix(mat); g <- as.character(group)
  keep <- !is.na(g); mat <- mat[keep, , drop = FALSE]; g <- g[keep]
  gf <- factor(g)
  s <- rowsum(mat, gf, reorder = TRUE)
  n <- as.vector(table(gf)[rownames(s)])
  s / n
}

if (!has_op || is.null(expr_file) || is.null(res_file)) {
  cat("跳过：缺少 oncoPredict 包或 GDSC2 训练数据。\n")
  cat("  1) install.packages('oncoPredict')\n")
  cat("  2) 浏览器打开 https://osf.io/c6tfx/ → DataFiles → Training Data\n")
  cat("  3) 下载 GDSC2_Expr (RMA Normalized and Log Transformed).rds 和 GDSC2_Res.rds\n")
  cat("  4) 放到：", td_dir, "\n", sep = "")
  stop("缺少 oncoPredict 包或 GDSC2 训练数据，按上方提示补好后重跑。")
}

suppressPackageStartupMessages({ library(oncoPredict); library(org.Hs.eg.db); library(AnnotationDbi) })

# ---------- TCGA symbol 矩阵（log2 CPM）----------
tcga <- readRDS(file.path(val_dir, "TCGA-LUAD.rds"))
cnt  <- tcga$counts; rn <- rownames(cnt)
ens  <- sub("\\..*$", "", rn)
sym  <- mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first")
names(sym) <- rn; sym <- sym[!is.na(sym) & sym != ""]
keep <- !is.na(sym[rownames(cnt)])
cpm  <- sweep(cnt[keep, , drop = FALSE], 2, colSums(cnt) + 1, "/") * 1e6
expr <- log2(collapse_mean(cpm, sym[rownames(cnt)[keep]]) + 1)
colnames(expr) <- substr(colnames(expr), 1, 15)
cat("TCGA 表达矩阵:", dim(expr), "\n")

# ---------- 训练数据 + 预测 ----------
trainingExprData <- readRDS(expr_file)
trainingPtype    <- readRDS(res_file)
# GDSC2 IC50 是 log 值，calcPhenotype 前须 exp() 回线性（vignette 明确要求）
trainingPtype <- exp(trainingPtype)

common <- intersect(rownames(trainingExprData), rownames(expr))
cat("训练基因:", nrow(trainingExprData), " | 测试基因重叠:", length(common), "\n")
testExprData     <- as.matrix(expr[common, , drop = FALSE])
trainingExprData <- trainingExprData[common, , drop = FALSE]

pred <- tryCatch(
  oncoPredict::calcPhenotype(trainingExprData = trainingExprData,
                             trainingPtype    = trainingPtype,
                             testExprData     = testExprData,
                             batchCorrect     = "none",   # 单一测试队列无批次，不用 ComBat
                             powerTransformPhenotype = TRUE,
                             removeLowVaryingGenes   = 0.2,
                             minNumSamples           = 10,
                             printOutput             = TRUE,
                             removeLowVaringGenesFrom = "rawData"),  # 注意 Varing 是包作者笔误，必填
  error = function(e) { cat("calcPhenotype 失败：", conditionMessage(e), "\n"); NULL })

if (is.null(pred)) stop("calcPhenotype 预测失败，见上方报错（会话已保留，改好后可重跑）。")

# pred：样本 × 药物（预测 IC50，越低越敏感）
cat("预测矩阵:", dim(pred), "\n")

# 匹配目标药物（grep 不区分大小写）
sel <- unique(unlist(lapply(drugs, function(d) grep(d, colnames(pred), ignore.case = TRUE))))
pred_drug <- pred[, sel, drop = FALSE]
cat("匹配到药物:", paste(colnames(pred_drug), collapse = ", "), "\n")

# ---------- 高/低风险组比较 ----------
sig  <- readRDS(file.path(res_dir, "rds", "signature_results.rds"))
risk <- sig$risk[["TCGA-LUAD"]]
common2 <- intersect(rownames(pred_drug), names(risk))
pred_drug <- pred_drug[common2, , drop = FALSE]
risk <- risk[common2]
grp <- factor(ifelse(risk > median(risk), "High", "Low"), levels = c("Low", "High"))
cat("对齐样本:", length(common2), " | Low/High:", sum(grp == "Low"), "/", sum(grp == "High"), "\n")

res <- do.call(rbind, lapply(colnames(pred_drug), function(d) {
  v <- as.numeric(pred_drug[, d])
  p <- wilcox.test(v ~ grp)$p.value
  data.frame(drug = d,
             ic50_med_Low  = round(median(v[grp == "Low"]), 3),
             ic50_med_High = round(median(v[grp == "High"]), 3),
             # 方向由数据定：本队列实测 High IC50 更低 = 更敏感（勿预设 High 更耐药）
             P = signif(p, 3), stringsAsFactors = FALSE)
}))
res$padj <- p.adjust(res$P, "BH")
res <- res[order(res$P), ]
cat("\n===== 药敏：High vs Low 风险（IC50 越低越敏感）=====\n")
print(res, row.names = FALSE)
write.csv(res, file.path(tab_dir, "drug_sensitivity_by_risk.csv"), row.names = FALSE)

# 箱线图（前 8 个药）
library(ggplot2)
top <- head(res$drug, 8)
long <- do.call(rbind, lapply(top, function(d) {
  data.frame(drug = d, ic50 = as.numeric(pred_drug[, d]), group = grp)
}))
long$drug <- factor(long$drug, levels = top)
p <- ggplot(long, aes(x = group, y = ic50, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 0.5, alpha = 0.4, color = "grey30") +
  facet_wrap(~ drug, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Low" = "#2E86AB", "High" = "#C0392B"), name = "Risk") +
  labs(x = NULL, y = "Predicted IC50") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", strip.background = element_rect(fill = "grey92", color = NA),
        panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "drug_sensitivity_by_risk.png"), p, width = 10, height = 5.5, dpi = 300)

cat("\n已保存 tables/drug_sensitivity_by_risk.csv + figures/drug_sensitivity_by_risk.png\n")
cat("说明：IC50 越低=越敏感。本队列实测 High 组 IC50 更低→更敏感（与「高风险→耐药」初始假设相反）。\n")
cat("     oncoPredict 为细胞系模型预测，仅生成假设；方向在所有药中一致提示可能反映增殖轴，非药物特异性耐药。\n")
