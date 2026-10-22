# ============================================================
# 47_tmb_cold_quadrant.R
# 目的：把「高 TMB + 免疫冷」的分离做成图（2×2 象限：TMB 高/低 × CYT 冷/热）
#   验证高风险组富集在「高 TMB × 免疫冷」象限——呼应讨论：TMB 单独不足以预测 ICB，
#   须结合免疫表型（基因型 TMB × 免疫表型 CYT 分离）
# 输入：TCGA-LUAD MAF（TMB）+ TCGA-LUAD 表达（CYT）+ 风险分
# 输出：figures/tmb_cyt_quadrant.png（散点）+ tmb_cyt_quadrant_tile.png（象限瓦片）
#       tables/tmb_cyt_quadrant_summary.csv
# 依赖：data.table / org.Hs.eg.db / AnnotationDbi / ggplot2
# 运行：本地 RStudio
# ============================================================

suppressPackageStartupMessages({
  library(org.Hs.eg.db); library(AnnotationDbi); library(ggplot2)
})

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures")
maf_dir <- file.path(base, "TCGA-LUAD_MAF")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

has_dt <- requireNamespace("data.table", quietly = TRUE)
if (!has_dt) cat("【缺包】data.table -> install.packages('data.table')（可选，快读大文件）\n")

# ---------- 1) TMB（与脚本 46 同）----------
maf_file <- list.files(maf_dir, pattern = "data_mutations", full.names = TRUE)
if (length(maf_file) == 0) stop("找不到 MAF，请把 data_mutations.txt 放到 ", maf_dir)
con <- file(maf_file[1], "r"); hdr <- readLines(con, n = 1); close(con)
all_cols <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
sel <- intersect(c("Hugo_Symbol","Tumor_Sample_Barcode","Variant_Classification"), all_cols)
if (has_dt) {
  maf <- as.data.frame(data.table::fread(maf_file[1], select = sel, sep = "\t",
                                         header = TRUE, quote = "", showProgress = TRUE))
} else {
  maf <- read.delim(maf_file[1], sep = "\t", header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "")
  maf <- maf[, intersect(sel, colnames(maf)), drop = FALSE]
}
nonsyn_cls <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins",
                "In_Frame_Del","In_Frame_Ins","Splice_Site","Translation_Start_Site",
                "Nonstop_Mutation","Start_Codon_Del","Start_Codon_Ins","Stop_Codon_Del","Stop_Codon_Ins")
maf$is_nonsyn <- maf$Variant_Classification %in% nonsyn_cls
maf$sid       <- substr(maf$Tumor_Sample_Barcode, 1, 15)

# ---------- 2) 风险 ----------
d   <- readRDS(file.path(val_dir, "bulk_program_scores.rds"))
sig <- readRDS(file.path(res_dir, "rds", "signature_results.rds"))
surv <- d$surv[["TCGA-LUAD"]]; risk <- sig$risk[["TCGA-LUAD"]]
common <- intersect(rownames(surv), names(risk))
rk <- risk[common]

nsc <- table(maf$sid[maf$is_nonsyn])
tmb_cnt <- setNames(rep(0L, length(common)), common)
mm <- intersect(names(nsc), names(tmb_cnt)); tmb_cnt[mm] <- as.integer(nsc[mm])
tmb <- tmb_cnt / 38

# ---------- 3) CYT（与脚本 45 同口径：GZMA+PRF1 z-score 均值）----------
collapse_mean <- function(mat, group) {
  mat <- as.matrix(mat); g <- as.character(group)
  k <- !is.na(g); mat <- mat[k, , drop = FALSE]; g <- g[k]
  gf <- factor(g); s <- rowsum(mat, gf, reorder = TRUE)
  n <- as.vector(table(gf)[rownames(s)]); s / n
}
score <- function(mat, genes) {
  g <- intersect(genes, rownames(mat)); if (length(g) < 2) return(NULL)
  z <- t(scale(t(mat[g, , drop = FALSE]))); colMeans(z, na.rm = TRUE)
}
tcga <- readRDS(file.path(val_dir, "TCGA-LUAD.rds"))
cnt <- tcga$counts; rn <- rownames(cnt)
ens <- sub("\\..*$", "", rn)
sym <- mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first")
names(sym) <- rn; sym <- sym[!is.na(sym) & sym != ""]
keep <- !is.na(sym[rownames(cnt)])
tcga_sym <- collapse_mean(cnt[keep, , drop = FALSE], sym[rownames(cnt)[keep]])
colnames(tcga_sym) <- substr(colnames(tcga_sym), 1, 15)
cyt <- score(tcga_sym, c("GZMA", "PRF1"))[common]

# ---------- 4) 象限 ----------
df <- data.frame(sample = common, tmb = tmb, cyt = cyt, risk = rk)
df$risk_grp <- factor(ifelse(df$risk > median(df$risk), "High", "Low"), levels = c("Low", "High"))
df$tmb_q <- factor(ifelse(df$tmb > median(df$tmb), "TMB-high", "TMB-low"),
                   levels = c("TMB-low", "TMB-high"))
df$cyt_q <- factor(ifelse(df$cyt > median(df$cyt), "CYT-hot", "CYT-cold"),
                   levels = c("CYT-cold", "CYT-hot"))
df$quadrant <- paste(df$tmb_q, df$cyt_q, sep = " / ")

# Fisher：高风险 vs 低风险 × 「TMB-high + CYT-cold」象限 vs 其他
target <- df$tmb_q == "TMB-high" & df$cyt_q == "CYT-cold"
ft <- fisher.test(table(df$risk_grp, target))
cat("Fisher（高/低风险 × TMB-high+CYT-cold 象限）: P =", signif(ft$p.value, 3),
    "| OR =", round(ft$estimate, 3), "\n")

qsum <- do.call(rbind, lapply(split(df, df$quadrant), function(x) {
  data.frame(quadrant = x$quadrant[1], n = nrow(x),
             high_frac = round(mean(x$risk_grp == "High"), 3),
             med_risk = round(median(x$risk), 3),
             med_tmb = round(median(x$tmb), 2), med_cyt = round(median(x$cyt), 3))
}))
qsum <- qsum[order(qsum$med_risk, decreasing = TRUE), ]
cat("\n===== 四象限汇总 =====\n"); print(qsum, row.names = FALSE)
write.csv(qsum, file.path(tab_dir, "tmb_cyt_quadrant_summary.csv"), row.names = FALSE)

# ---------- 5) 图 ----------
p1 <- ggplot(df, aes(x = tmb, y = cyt, color = risk_grp)) +
  geom_hline(yintercept = median(df$cyt), linetype = 2, color = "grey50") +
  geom_vline(xintercept = median(df$tmb), linetype = 2, color = "grey50") +
  geom_point(size = 1.6, alpha = 0.7) +
  scale_color_manual(values = c("Low" = "#2E86AB", "High" = "#C0392B"), name = "Risk") +
  annotate("text", x = max(df$tmb), y = max(df$cyt), hjust = 1, vjust = 1, label = "TMB-high · hot", size = 3) +
  annotate("text", x = max(df$tmb), y = min(df$cyt), hjust = 1, vjust = 0, label = "TMB-high · cold", size = 3) +
  annotate("text", x = min(df$tmb), y = max(df$cyt), hjust = 0, vjust = 1, label = "TMB-low · hot", size = 3) +
  annotate("text", x = min(df$tmb), y = min(df$cyt), hjust = 0, vjust = 0, label = "TMB-low · cold", size = 3) +
  labs(x = "TMB (mut/Mb)", y = "Cytolytic score (CYT, GZMA+PRF1)",
       title = "TCGA-LUAD: high-risk concentrates in TMB-high × immune-cold") +
  theme_bw(base_size = 12) + theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "tmb_cyt_quadrant.png"), p1, width = 7, height = 6, dpi = 300)

tile <- expand.grid(tmb_q = c("TMB-low", "TMB-high"), cyt_q = c("CYT-cold", "CYT-hot"))
tile$med_risk  <- sapply(seq_len(nrow(tile)), function(i)
  median(df$risk[df$tmb_q == tile$tmb_q[i] & df$cyt_q == tile$cyt_q[i]]))
tile$high_frac <- sapply(seq_len(nrow(tile)), function(i)
  mean(df$risk_grp[df$tmb_q == tile$tmb_q[i] & df$cyt_q == tile$cyt_q[i]] == "High"))
p2 <- ggplot(tile, aes(x = cyt_q, y = tmb_q, fill = med_risk)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("med risk %.2f\nHigh %.0f%%", med_risk, high_frac * 100)),
            color = "white", size = 3.4, fontface = "bold") +
  scale_fill_gradient2(low = "#2E86AB", mid = "white", high = "#C0392B",
                       midpoint = 0, name = "Median risk") +
  labs(x = "Immune phenotype", y = "TMB", title = "Median signature risk by quadrant") +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "tmb_cyt_quadrant_tile.png"), p2, width = 6, height = 4.5, dpi = 300)

cat("\n已保存 tmb_cyt_quadrant.png + tmb_cyt_quadrant_tile.png + tmb_cyt_quadrant_summary.csv\n")
cat("说明：CYT 冷=低细胞溶解活性；若高风险富集于「高 TMB × 冷」象限，即实证『TMB 高但免疫冷』的分离。\n")
