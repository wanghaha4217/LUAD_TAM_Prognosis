# ============================================================
# 74_spp1_prognosis.R
# 目的：回答「SPP1 对预后有没有作用」——补测 SPP1 单基因（而非 P1 程序）
#       在 4 个 bulk 队列的单变量 Cox。
#       背景：正文签名用的是 P1 程序（SPP1+C1q+细胞周期 前50基因等权均值），
#             P1 因单变量无预后 + 与 P2 共线被排除；但这不能推出 SPP1 单基因无预后。
# 测什么（每队列单变量 Cox，表达 z-score 后 = 每 +1 SD 的 HR）：
#   (a) SPP1 单基因
#   (b) C1QA / C1QB / C1QC（P1 的补体组分，检验是否稀释了 P1 信号）
#   (c) CXCL9 单基因 + SPP1−CXCL9 极性差值（方向一撞车文献的 CXCL9/SPP1 轴，
#       越高越偏 SPP1 化/免疫排斥）
# 输入：validation_cohorts/bulk_symbol_matrices.rds（4 队列符号矩阵 + 生存）
# 输出：tables/spp1_prognosis.csv（透明可查）
# 运行：本地 RStudio
# ============================================================

suppressPackageStartupMessages({ library(survival) })

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
tab_dir <- file.path(base, "results", "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

bulk <- readRDS(file.path(val_dir, "bulk_symbol_matrices.rds"))
mat_list  <- bulk$mat
surv_list <- bulk$surv
cohorts   <- intersect(names(mat_list), names(surv_list))
cat("可用队列:", paste(cohorts, collapse = ", "), "\n\n")

# ---------- 工具：单变量 Cox（表达已 z-score，HR = 每 +1 SD） ----------
univ_cox <- function(x, surv, cohort, label) {
  df <- data.frame(time = surv$time, event = surv$event, x = x)
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 10 || sd(df$x) == 0) return(NULL)
  fit <- tryCatch(coxph(Surv(time, event) ~ x, data = df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ci <- summary(fit)$conf.int[1, c(1, 3, 4)]
  data.frame(cohort = cohort, feature = label, N = nrow(df),
             HR = ci[1], lower = ci[2], upper = ci[3],
             P = summary(fit)$coefficients[1, "Pr(>|z|)"],
             Cindex = summary(fit)$concordance[1],
             stringsAsFactors = FALSE)
}

zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

rows <- list()
for (g in cohorts) {
  m <- mat_list[[g]]
  sv <- surv_list[[g]]

  # 样本对齐：表达矩阵（列名）与生存数据（行名）取交集
  common <- intersect(colnames(m), rownames(sv))
  cat(sprintf("[%s] 表达样本 %d，生存样本 %d，交集 %d\n",
              g, ncol(m), nrow(sv), length(common)))
  if (length(common) < 10) { cat("  交集样本过少，跳过\n"); next }
  m  <- m[, common, drop = FALSE]
  sv <- sv[common, , drop = FALSE]

  gsym <- rownames(m)

  # 基因可用性
  avail <- function(gene) if (gene %in% gsym) gene else NA
  genes <- c("SPP1", "C1QA", "C1QB", "C1QC", "CXCL9")
  expr <- lapply(genes, function(gene) {
    if (is.na(avail(gene))) return(NULL)
    zscore(as.numeric(m[gene, , drop = TRUE]))
  })
  names(expr) <- genes

  # (a) SPP1 单基因
  if (!is.null(expr$SPP1)) {
    rows[[length(rows) + 1]] <- univ_cox(expr$SPP1, sv, g, "SPP1")
  }
  # (b) 补体组分
  for (gene in c("C1QA", "C1QB", "C1QC")) {
    if (!is.null(expr[[gene]])) {
      rows[[length(rows) + 1]] <- univ_cox(expr[[gene]], sv, g, gene)
    }
  }
  # (c) CXCL9 单基因 + SPP1−CXCL9 极性差值
  if (!is.null(expr$CXCL9)) {
    rows[[length(rows) + 1]] <- univ_cox(expr$CXCL9, sv, g, "CXCL9")
    if (!is.null(expr$SPP1)) {
      pol <- expr$SPP1 - expr$CXCL9   # 越高越 SPP1 化
      rows[[length(rows) + 1]] <- univ_cox(pol, sv, g, "SPP1_minus_CXCL9")
    }
  }
  cat("[", g, "] 完成\n", sep = "")
}

out <- do.call(rbind, rows)
out <- out[order(out$cohort, out$feature), ]
print(out, row.names = FALSE)

write.csv(out, file.path(tab_dir, "spp1_prognosis.csv"), row.names = FALSE)
cat("\n已写出: results/tables/spp1_prognosis.csv\n")
