# ============================================================
# 40_program_enrichment.R
# 对 NMF 10 程序的 top 基因做功能富集(Hallmark + GO:BP),
# 用数据驱动方式给 P1~P10 打功能标签(尤其敲定 P5/P6)。
# 依赖(均已装): clusterProfiler / org.Hs.eg.db / msigdbr
# ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(dplyr)
  library(ggplot2)
})

base    <- "."
tab_dir <- file.path(base, "results", "tables")
fig_dir <- file.path(base, "results", "figures")

# ---- 1) 读 top-50 基因(宽表:列 P1~P10,每列 = 该程序 top-50 基因) ----
gdf   <- read.csv(file.path(base, "program_genes_top50.csv"), check.names = FALSE, stringsAsFactors = FALSE)
progs <- colnames(gdf)
cat("程序:", paste(progs, collapse = ", "), "\n")

# ---- 2) Hallmark 基因集(symbol 版 TERM2GENE) ----
hm     <- msigdbr(species = "Homo sapiens", collection = "H")
hm_t2g <- hm[, c("gs_name", "gene_symbol")]
cat("Hallmark 基因集数:", length(unique(hm_t2g$gs_name)), "\n")

# ---- 3) 逐程序富集 ----
all_res <- list()
for (p in progs) {
  genes <- gdf[[p]]
  genes <- genes[!is.na(genes) & genes != ""]
  cat("\n===== ", p, " (", length(genes), " 基因) =====\n", sep = "")

  # Hallmark ORA(用 symbol 直接 enricher)
  hk <- tryCatch(
    enricher(gene = genes, TERM2GENE = hm_t2g, pvalueCutoff = 0.05,
             pAdjustMethod = "BH", minGSSize = 5, maxGSSize = 500),
    error = function(e) { cat("  Hallmark 出错:", e$message, "\n"); NULL })
  if (!is.null(hk) && nrow(as.data.frame(hk)) > 0) {
    hk_df <- as.data.frame(hk)
    cat("  [Hallmark] top 富集:\n")
    print(hk_df[1:min(6, nrow(hk_df)), c("Description", "p.adjust", "Count")])
    all_res[[paste0(p, "__H")]] <- cbind(program = p, source = "Hallmark", hk_df)
  } else {
    cat("  [Hallmark] 无显著富集\n")
  }

  # GO:BP(转 entrez 再 enrichGO)
  eg <- tryCatch(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
                 error = function(e) NULL)
  if (!is.null(eg) && nrow(eg) > 0) {
    go <- tryCatch(
      enrichGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
               pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE),
      error = function(e) { cat("  GO 出错:", e$message, "\n"); NULL })
    if (!is.null(go) && nrow(as.data.frame(go)) > 0) {
      go_df <- as.data.frame(go)
      cat("  [GO:BP] top 富集:\n")
      print(go_df[1:min(5, nrow(go_df)), c("Description", "p.adjust", "Count")])
      all_res[[paste0(p, "__GO")]] <- cbind(program = p, source = "GO:BP", go_df)
    } else {
      cat("  [GO:BP] 无显著富集\n")
    }
  }
}

# ---- 4) 合并保存 ----
if (length(all_res) > 0) {
  out <- do.call(rbind, lapply(all_res, function(df) {
    keep <- intersect(c("program", "source", "ID", "Description", "GeneRatio",
                        "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"),
                      names(df))
    df[, keep, drop = FALSE]
  }))
  rownames(out) <- NULL
  write.csv(out, file.path(tab_dir, "program_enrichment.csv"), row.names = FALSE)
  cat("\n已保存:", file.path(tab_dir, "program_enrichment.csv"), "(", nrow(out), " 行)\n")

  # ---- 5) dotplot:每程序 top-3 Hallmark 术语 ----
  hk_only <- out[out$source == "Hallmark", ]
  if (nrow(hk_only) > 0) {
    top3 <- hk_only %>%
      group_by(program) %>%
      slice_min(order_by = p.adjust, n = 3) %>%
      ungroup() %>%
      mutate(Description = sub("^HALLMARK_", "", Description))
    top3$program <- factor(top3$program, levels = progs)
    top3$Description <- factor(top3$Description, levels = rev(unique(top3$Description)))
    g <- ggplot(top3, aes(x = program, y = Description, size = Count, color = -log10(p.adjust))) +
      geom_point() +
      scale_color_gradient(low = "steelblue", high = "firebrick") +
      scale_size_continuous(range = c(2, 9)) +
      labs(x = "NMF program", y = NULL, color = "-log10(p.adjust)",
           size = "Gene count",
           title = "Program functional annotation (top-3 Hallmark terms)") +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_blank())
    ggsave(file.path(fig_dir, "program_enrichment_dotplot.png"), g,
           width = 11, height = 7, dpi = 300)
    cat("已保存:", file.path(fig_dir, "program_enrichment_dotplot.png"), "\n")
  }
} else {
  cat("无任何显著富集结果\n")
}
