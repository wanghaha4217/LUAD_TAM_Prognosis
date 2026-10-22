# ============================================================
# 10_build_discovery_seurat.R
# 目的：从 GSE243013 的 7.1G 大矩阵里，只取「LUAD 巨噬」细胞，建 Seurat 对象
# 依赖：Matrix / Seurat / data.table
# 内存：readMM 全读需要大内存，请在云服务器跑（本地 31.6G 会 OOM）
# 运行：本地/云端 RStudio，全选运行
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(data.table)
})

base    <- "."
out_dir <- file.path(base, "discovery")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mtx_file     <- file.path(base, "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz")
gene_file    <- file.path(base, "GSE243013_genes.csv(1).gz")   # 若你已改名，改这里
barcode_file <- file.path(base, "GSE243013_barcodes.csv.gz")
meta_file    <- file.path(base, "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz")

# ---------- 1. 读 metadata，锁定目标细胞（LUAD 巨噬） ----------
meta <- fread(meta_file)          # data.table 自动解压 .gz
cat("metadata:", nrow(meta), "细胞 ×", ncol(meta), "字段\n")

is_luad  <- meta$cancer_type == "LUAD"
is_macro <- grepl("φ", meta$sub_cell_type, fixed = TRUE)   # 作者用 Mφ_* 标注巨噬亚型
target_cells <- meta$cellID[is_luad & is_macro]

cat("LUAD 巨噬细胞:", length(target_cells),
    " | 患者:", uniqueN(meta$sampleID[is_luad & is_macro]), "\n")
cat("--- LUAD 巨噬亚型（sub_cell_type）分布 ---\n")
print(table(meta$sub_cell_type[is_luad & is_macro]))

if (length(target_cells) == 0) {
  stop("没匹配到任何巨噬细胞。先跑 table(meta$sub_cell_type[meta$cancer_type=='LUAD']) 看实际取值，把结果贴给我。")
}

# ---------- 2. barcodes（= 矩阵行标签）与基因名 ----------
barcodes <- fread(barcode_file, header = TRUE)$barcode
genes    <- fread(gene_file, header = TRUE)$geneSymbol
cat("barcodes:", length(barcodes), " | 基因:", length(genes), "\n")

target_rows <- match(target_cells, barcodes)   # 矩阵行 = barcodes 顺序
if (anyNA(target_rows)) {
  stop("有 ", sum(is.na(target_rows)), " 个 cellID 在 barcodes.csv 里找不到。把不匹配的 cellID 贴给我看。")
}
cat("目标矩阵行号:", length(target_rows), "（全部匹配成功）\n")

# ---------- 3. 读大矩阵（cells × genes），只取目标行 ----------
# 注意：该矩阵和标准 10x 转置了 —— 头部 1254749 31831 表示 行=细胞、列=基因
# 所以下面先按行子集（取细胞），再 t() 转回 基因×细胞
cat("读取 7.1G 矩阵（1,254,749 × 31,831，非零元 ~20 亿）... 需要大内存，请耐心\n")
M <- readMM(gzfile(mtx_file))       # dgCMatrix，行=细胞，列=基因
cat("全矩阵维度:", dim(M), "\n")

M <- M[target_rows, , drop = FALSE]
gc()

# 转置成 基因 × 细胞（Seurat 惯例）
M <- t(M)
rownames(M) <- make.unique(genes)   # 若基因名有重复，会加 .1 后缀；一般不会
colnames(M) <- target_cells
cat("子集矩阵（基因 × 细胞）:", dim(M), "\n")

# ---------- 4. 建 Seurat 对象 ----------
sub_meta <- as.data.frame(meta[is_luad & is_macro])
rownames(sub_meta) <- sub_meta$cellID

so <- CreateSeuratObject(counts = M, meta.data = sub_meta,
                         project = "GSE243013_LUAD_TAM")
cat("Seurat 对象:", nrow(so), "基因 ×", ncol(so), "细胞\n")

# ---------- 5. 保存 ----------
saveRDS(so, file.path(out_dir, "GSE243013_LUAD_TAM_seurat.rds"))
cat("已保存:", file.path(out_dir, "GSE243013_LUAD_TAM_seurat.rds"), "\n")
