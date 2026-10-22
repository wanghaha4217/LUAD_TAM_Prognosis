# ============================================================
# 62_cellchat_prep.R
# 目的：从完整 GSE243013 提取 LUAD 免疫细胞对象，供 CellChat（脚本 63）用
#   - 巨噬保留 11 亚型（Mφ_*），非巨噬合并到 ~15 类（Treg/CD4T/CD8T/CD8T_Tex/NK/gdT/B/Plasma/cDC/Neu/Mast/ILC）
#   - 带 pathological_response（non-MPR vs MPR/pCR）分组
# 运行：云端 RStudio（内存 ~30-40GB，读 7.1GB mtx）
# 输出：./GSE243013_LUAD_immune_cellchat.rds
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

GEO <- "."
DIS <- "."
out <- file.path(DIS, "GSE243013_LUAD_immune_cellchat.rds")

# ---------- 1. 读 metadata，筛 LUAD ----------
meta <- fread(file.path(GEO, "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz"))
meta_luad <- meta[cancer_type == "LUAD"]
cat("LUAD 免疫细胞数:", nrow(meta_luad), "（应 336685）\n")
cat("major_cell_type 分布:\n"); print(table(meta_luad$major_cell_type))

# ---------- 2. 细胞类型注释 ----------
# 巨噬 11 亚型：从巨噬 rds 的 cellID -> sub_cell_type
so_mac <- readRDS(file.path(DIS, "GSE243013_LUAD_TAM_clustered.rds"))
mac_map <- setNames(so_mac@meta.data$sub_cell_type, rownames(so_mac@meta.data))
rm(so_mac); gc()

celltype <- ifelse(meta_luad$cellID %in% names(mac_map),
                   as.character(mac_map[meta_luad$cellID]), NA_character_)

# 非巨噬：sub_cell_type 前缀合并（注意顺序，Treg 先于 CD4T、Tex 先于 CD8T）
idx <- is.na(celltype)
s <- meta_luad$sub_cell_type[idx]
mapped <- fcase(
  grepl("CD4T_Treg", s), "Treg",
  grepl("^CD4T", s), "CD4T",
  grepl("CD8T_NK-like", s), "NK",
  grepl("Tex", s), "CD8T_Tex",
  grepl("^CD8T", s), "CD8T",
  grepl("^NK", s), "NK",
  grepl("^T_gdT", s), "gdT",
  grepl("^(Bm|Bn)_", s), "B cell",
  grepl("^Plasma", s), "Plasma",
  grepl("^cDC", s), "cDC",
  grepl("^Neu", s), "Neutrophil",
  grepl("^Mast", s), "Mast",
  grepl("^ILC", s), "ILC",
  default = s
)
celltype[idx] <- mapped
meta_luad$celltype <- celltype
cat("\ncelltype 分布:\n"); print(sort(table(celltype), decreasing = TRUE))

# ---------- 3. 读 counts mtx，筛 LUAD 细胞 ----------
barcodes <- readLines(gzfile(file.path(GEO, "GSE243013_barcodes.csv.gz")))
barcodes <- barcodes[-1]                  # 去掉表头 "barcode"
stopifnot(length(barcodes) == 1254749)
# barcodes 可能是 cellID（P266-xxx-1）或纯 barcode（xxx-1），两种都兼容
if (any(barcodes %in% meta_luad$cellID)) {
  keep_rows <- match(meta_luad$cellID, barcodes)          # 直接按 cellID 对齐
} else {
  bc <- sub("^[^-]+-", "", meta_luad$cellID)               # 去 sampleID 前缀
  keep_rows <- match(bc, barcodes)
}
stopifnot(!anyNA(keep_rows))
cat("\n匹配到的 LUAD 细胞行数:", length(keep_rows), "（应 = 336685）\n")

mtx <- readMM(gzfile(file.path(GEO, "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz")))
cat("完整 mtx 维度:", dim(mtx), "（细胞 x 基因）\n")
mtx <- mtx[keep_rows, , drop = FALSE]   # 筛 LUAD 细胞行
mtx <- t(mtx)                           # 转置成 基因 x 细胞

genes <- readLines(gzfile(file.path(GEO, "GSE243013_genes.csv(1).gz")))
genes <- genes[-1]                        # 去掉表头 "geneSymbol"
stopifnot(length(genes) == 31831)
if (anyDuplicated(genes)) {
  cat("geneSymbol 有重复，make.unique 处理\n")
  genes <- make.unique(genes)
}
rownames(mtx) <- genes
colnames(mtx) <- meta_luad$cellID       # 列顺序 = meta_luad 顺序，直接对齐
rm(barcodes); gc()

# ---------- 4. 构建 Seurat 对象 ----------
md <- as.data.frame(meta_luad)
rownames(md) <- md$cellID
so <- CreateSeuratObject(counts = mtx, meta.data = md, min.cells = 0, min.features = 0)
cat("\n对象维度:", nrow(so), "基因 x", ncol(so), "细胞\n")
cat("Idents 备用 celltype 列，pathological_response 分组:\n")
print(table(so$celltype, so$pathological_response))

# ---------- 5. 存 ----------
saveRDS(so, out)
cat("\n已保存:", out, "\n")
cat("object.size =", format(object.size(so), units = "auto"), "\n")
