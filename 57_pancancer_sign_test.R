# ============================================================
# 57_pancancer_sign_test.R
# 目的：泛癌「方向一致性」符号检验（二项检验），固化正文 4.9 引用的 P 值。
#   正文 4.9 写的 P8/P9=0.0037、P4=0.036 此前是手算，无脚本无 CSV——
#   审稿人要可复现，本脚本把 9 个程序的 HR>1 计数 + 双侧二项检验固化成 CSV。
# 输入：results/tables/pancancer_program_hr.csv（脚本 52 已生成）
# 输出：results/tables/pancancer_sign_test.csv
# 运行：本地 RStudio，秒级（纯读 CSV + binom.test，不碰表达矩阵）
# ============================================================
suppressPackageStartupMessages({})

base    <- "."
tab_dir <- file.path(base, "results", "tables")

d <- read.csv(file.path(tab_dir, "pancancer_program_hr.csv"),
              stringsAsFactors = FALSE)

progs <- sort(unique(d$program))
res <- data.frame(program = progs,
                  n_cancers = NA_integer_,
                  n_HR_gt_1 = NA_integer_,
                  n_HR_lt_1 = NA_integer_,
                  sign_test_P = NA_real_,
                  direction = NA_character_,
                  stringsAsFactors = FALSE)

for (p in progs) {
  h <- d$HR[d$program == p]
  n_up   <- sum(h > 1)
  n_down <- sum(h < 1)
  n      <- n_up + n_down
  bt <- binom.test(n_up, n, p = 0.5, alternative = "two.sided")
  res[res$program == p, c("n_cancers", "n_HR_gt_1", "n_HR_lt_1", "sign_test_P")] <-
    c(n, n_up, n_down, bt$p.value)
}

res$direction <- ifelse(res$n_HR_gt_1 > res$n_HR_lt_1, "adverse", "protective")
res <- res[order(res$sign_test_P), ]
res$sign_test_P <- signif(res$sign_test_P, 3)

write.csv(res, file.path(tab_dir, "pancancer_sign_test.csv"), row.names = FALSE)
print(res, row.names = FALSE)

cat("\n说明：sign_test_P < 0.05 = 该程序在 28 癌种中 HR 方向显著一致（倾向 >1 或 <1）。\n")
cat("正文 4.9 只引用「方向显著不良」的程序（n_HR_gt_1 占多数 且 P<0.05）。\n")
cat("预期（与正文一致）：P8、P9 各 22/28 不良 P=0.0037；P4 20/28 P=0.036。\n")
cat("已保存:", file.path(tab_dir, "pancancer_sign_test.csv"), "\n")
