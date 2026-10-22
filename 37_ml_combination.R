# ============================================================
# 37_ml_combination.R  (v3 — 方法改为 LASSO + RSF + XGBoost，与正文 4.8 / Fig 7E 一致)
# 目的：方法敏感性 —— 随机生存森林(RSF) 与 XGBoost(Cox) 与 LASSO-Cox 横向比，
#       输入同为「9 程序打分（P2-P10）」，训练 TCGA → 3 外部队列 C-index
# v3 修改：CoxBoost 换成 XGBoost。正文 4.8 + Fig 7E 图注均写「LASSO vs RSF vs XGBoost」，
#         且正文已引用 XGBoost 的具体数字（来自原 59_xgboost_ml.R）。CoxBoost 从未跑通
#         （包未装，旧 CSV 里没有 CoxBoost 列），故弃 CoxBoost、改用 XGBoost 对齐正文。
#         XGBoost 逻辑从原 59_xgboost_ml.R 并入。
# v2 修改：ci() 返回 C + se(C)，柱上标数值 + 显著性星号（单侧 C>0.5）。
# 依赖：randomForestSRC、xgboost（缺则跳过对应模型）；glmnet 已装
# 运行：本地 RStudio；首次需 install.packages(c("randomForestSRC","xgboost"))
# ============================================================

suppressPackageStartupMessages({ library(survival) })

base    <- "."
val_dir <- file.path(base, "validation_cohorts")
res_dir <- file.path(base, "results")
tab_dir <- file.path(res_dir, "tables")
fig_dir <- file.path(res_dir, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

has_rsf <- requireNamespace("randomForestSRC", quietly = TRUE)
has_xgb <- requireNamespace("xgboost", quietly = TRUE)
if (!has_rsf) cat("【缺包】randomForestSRC -> install.packages('randomForestSRC')（跳过 RSF）\n")
if (!has_xgb) cat("【缺包】xgboost -> install.packages('xgboost')（跳过 XGBoost）\n")

d       <- readRDS(file.path(val_dir, "bulk_program_scores.rds"))
sig     <- readRDS(file.path(res_dir, "rds", "signature_results.rds"))
scores  <- d$scores
surv    <- d$surv
cohorts <- names(scores)

# 特征 = 9 程序（P2-P10，与训练一致，让各方法自选）
feats   <- colnames(scores[["TCGA-LUAD"]])
tr_mean <- sig$tr_mean
tr_sd   <- sig$tr_sd

scale_scores <- function(df, mu, sd) {
  df  <- as.matrix(df)
  cmn <- intersect(colnames(df), names(mu))
  sweep(sweep(df[, cmn, drop = FALSE], 2, mu[cmn], "-"), 2, sd[cmn], "/")
}
X <- lapply(scores, function(s) scale_scores(s, tr_mean, tr_sd))

x_tr <- X[["TCGA-LUAD"]]
t_tr <- surv[["TCGA-LUAD"]]$time
e_tr <- surv[["TCGA-LUAD"]]$event
ok_tr <- complete.cases(x_tr) & !is.na(t_tr) & !is.na(e_tr)
x_tr <- x_tr[ok_tr, , drop = FALSE]; t_tr <- t_tr[ok_tr]; e_tr <- e_tr[ok_tr]
cat("训练样本:", nrow(x_tr), " | 事件:", sum(e_tr), "\n")

# ---------- 工具：Harrell's C（coxph 自带，方向正确）+ se + C>0.5 单侧 P ----------
ci_full <- function(risk_named, surv_df) {
  rn <- rownames(surv_df)
  r  <- as.numeric(risk_named[rn])
  t  <- surv_df$time
  e  <- surv_df$event
  ok <- !is.na(r) & !is.na(t) & !is.na(e)
  if (sum(ok) < 10) {
    cat("  [ci] 可对齐样本只有", sum(ok), "个（<10），返回 NA\n", sep = "")
    return(c(C = NA_real_, se = NA_real_, p = NA_real_))
  }
  r <- r[ok]; t <- t[ok]; e <- e[ok]
  if (length(unique(r)) < 2) {
    cat("  [ci] 风险分无差异（unique<2），返回 NA\n")
    return(c(C = NA_real_, se = NA_real_, p = NA_real_))
  }
  cc <- summary(coxph(Surv(t, e) ~ r))$concordance
  C  <- as.numeric(cc[1]); se <- as.numeric(cc[2])
  c(C = C, se = se, p = pnorm((C - 0.5) / se, lower.tail = FALSE))  # 单侧 C>0.5
}

# ---------- 各方法风险向量（method -> cohort -> named risk）----------
RISK <- list(LASSO = sig$risk)

# ---------- 2. 随机生存森林（RSF）----------
if (has_rsf) {
  library(randomForestSRC)
  rsf_res <- tryCatch({
    tr_df <- data.frame(time = t_tr, event = e_tr, x_tr, check.names = FALSE)
    set.seed(123)
    rf <- rfsrc(Surv(time, event) ~ ., data = tr_df, ntree = 1000, importance = "none")
    out <- list("TCGA-LUAD" = setNames(as.numeric(rf$predicted), rownames(x_tr)))
    for (g in setdiff(cohorts, "TCGA-LUAD")) {
      vd <- data.frame(X[[g]], check.names = FALSE)
      out[[g]] <- setNames(as.numeric(predict(rf, newdata = vd)$predicted), rownames(X[[g]]))
    }
    out
  }, error = function(e) { cat("RSF 失败：", conditionMessage(e), "\n"); NULL })
  if (!is.null(rsf_res)) RISK$RSF <- rsf_res
}

# ---------- 3. XGBoost（survival:cox，与原 59_xgboost_ml.R 一致）----------
if (has_xgb) {
  library(xgboost)
  cat("xgboost 版本:", as.character(packageVersion("xgboost")), "\n")
  xgb_res <- tryCatch({
    dtr <- xgb.DMatrix(x_tr, label = e_tr)
    setinfo(dtr, "label_lower_bound", t_tr)
    setinfo(dtr, "label_upper_bound", t_tr)
    params <- list(objective = "survival:cox", eval_metric = "cox-nloglik",
                   max_depth = 3, eta = 0.05, subsample = 0.8, colsample_bytree = 0.8)
    nb <- 100
    cv <- tryCatch({
      set.seed(123)
      xgb.cv(data = dtr, nrounds = 200, nfold = 5, params = params, verbose = 0)
    }, error = function(e) { cat("  xgb.cv 失败（用默认 100 轮）:", conditionMessage(e), "\n"); NULL })
    if (!is.null(cv)) {
      ev <- cv$evaluation_log
      tc <- grep("^test_.*mean$", names(ev), value = TRUE)
      if (length(tc) > 0 && nrow(ev) > 0) {
        it   <- ev[[1]]; v <- ev[[tc[1]]]; imin <- which.min(v)
        nb   <- max(as.integer(it[imin]), 30L)
        cat("  测试损失最小在 iter", it[imin], "(", signif(v[imin], 5), ")\n")
      }
    }
    nb <- as.integer(nb[1]); cat("  采用 nrounds:", nb, "\n")
    set.seed(123)
    bst <- xgb.train(params = params, data = dtr, nrounds = nb, verbose = 0)
    lp <- lapply(X, function(xx)
      setNames(as.numeric(predict(bst, xgb.DMatrix(xx))), rownames(xx)))
    list(lp = lp, nrounds = nb)
  }, error = function(e) { cat("XGBoost 失败：", conditionMessage(e), "\n"); NULL })
  if (!is.null(xgb_res)) {
    RISK$XGBoost <- xgb_res$lp
    cat("XGBoost 最优 nrounds:", xgb_res$nrounds, "\n")
    ext <- sapply(setdiff(cohorts, "TCGA-LUAD"), function(g) {
      rr <- RISK$XGBoost[[g]]; s <- surv[[g]]
      rn <- rownames(s); r <- as.numeric(rr[rn])
      ok <- !is.na(r) & !is.na(s$time) & !is.na(s$event)
      if (sum(ok) < 10) return(NA_real_)
      as.numeric(summary(coxph(Surv(s$time[ok], s$event[ok]) ~ r[ok]))$concordance[1])
    })
    if (all(!is.na(ext)) && median(ext, na.rm = TRUE) < 0.5)
      cat("【提示】XGBoost 外部队列 C 中位 <0.5，预测方向可能相反，请检查。\n")
  }
}

# ---------- 输出：C-index + se + 显著性（长表）----------
methods <- intersect(c("LASSO", "RSF", "XGBoost"), names(RISK))
res_list <- lapply(methods, function(m) {
  do.call(rbind, lapply(cohorts, function(g) {
    r <- ci_full(RISK[[m]][[g]], surv[[g]])
    data.frame(cohort = g, method = m, Cindex = r["C"], se = r["se"],
               p_vs_05 = r["p"], stringsAsFactors = FALSE)
  }))
})
long <- do.call(rbind, res_list); rownames(long) <- NULL
long$star <- ifelse(long$p_vs_05 < 0.001, "***",
                    ifelse(long$p_vs_05 < 0.01, "**",
                           ifelse(long$p_vs_05 < 0.05, "*", "ns")))
long$lab  <- sprintf("%.3f%s", long$Cindex, ifelse(long$star == "ns", "", long$star))

# 宽表 CSV（cohort × method，覆盖旧版仅 LASSO/RSF 两列）
wide <- data.frame(cohort = cohorts, stringsAsFactors = FALSE)
for (i in seq_along(methods)) wide[[methods[i]]] <- res_list[[i]]$Cindex
write.csv(wide, file.path(tab_dir, "ml_combination_cindex.csv"), row.names = FALSE)

cat("\n===== 方法对比 C-index（TCGA 训练 + 3 外部队列；* = C 显著>0.5）=====\n")
print(long, row.names = FALSE)

suppressPackageStartupMessages(library(ggplot2))
long$cohort <- factor(long$cohort, levels = cohorts)
ymax <- max(long$Cindex, na.rm = TRUE) * 1.15
p <- ggplot(long, aes(cohort, Cindex, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = lab), position = position_dodge(0.8), vjust = -0.4, size = 3.2) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey50") +
  scale_fill_manual(values = c(LASSO = "#E41A1C", RSF = "#4DAF4A", XGBoost = "#377EB8")) +
  labs(x = NULL, y = "Harrell's C-index", fill = NULL,
       caption = "* P<0.05, ** P<0.01, *** P<0.001 (one-sided test C > 0.5)") +
  coord_cartesian(ylim = c(0.4, ymax)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        axis.text  = element_text(size = 11),
        axis.title = element_text(size = 13))
ggsave(file.path(fig_dir, "ml_combination_cindex.png"), p, width = 8, height = 5, dpi = 300)

cat("\n已保存 tables/ml_combination_cindex.csv + figures/ml_combination_cindex.png\n")
cat("说明：三种方法输入同为 9 程序打分；LASSO 内部选 7 程序，RSF/XGBoost 自行非参/梯度提升。\n")
cat("      TCGA 一栏为训练集内评估（偏乐观），看 3 个外部队列的一致性才是重点。\n")
