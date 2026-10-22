# Script manifest

Each row maps a script to the figure / table / section it produces.
Tier: core = main figures & tables; supp = supplementary figures cited in the main text; novelty = SPP1/regulon differentiation analyses.

| Script | Tier | Produces |
|---|---|---|
| `10_build_discovery_seurat.R` | 核心 · 主图主表 | Build the discovery Seurat object (GSE243013, 38,211 TAM) |
| `26_program_subtype_annotation.R` | 核心 · 主图主表 | Annotate 11 macrophage subtypes |
| `42b_nmf_one_k.R` | 核心 · 主图主表 | NMF k=10 factorization -> 10 programs + top-50 genes |
| `25_discovery_figures.R` | 核心 · 主图主表 | Fig 1: UMAP + NMF heatmap + program-subtype DotPlot |
| `40_program_enrichment.R` | 补充 · 主文引用的 Fig S | Fig S5: functional enrichment of the 10 programs |
| `00_fetch_validation_cohorts.R` | 核心 · 主图主表 | Fetch TCGA-LUAD / GEO validation data |
| `01_build_validation_cohorts.R` | 核心 · 主图主表 | Build the 3 validation cohort expression/survival matrices |
| `20_bulk_program_scores.R` | 核心 · 主图主表 | Score the 10 programs in every bulk cohort |
| `21_lasso_cox_signature.R` | 核心 · 主图主表 | LASSO-Cox -> 7-program signature (fixed coefficients) |
| `27_lasso_figure.R` | 核心 · 主图主表 | Fig 2: CV curve + coefficient paths + coefficient bar |
| `76_program_correlation.R` | 补充 · 主文引用的 Fig S | Fig S11: program-program correlation heatmap |
| `22_km_curves.R` | 核心 · 主图主表 | Fig 3: Kaplan-Meier by risk group in 4 cohorts |
| `23_cox_analysis.R` | 核心 · 主图主表 | Fig 4A + Table 2: multivariable Cox + C-index/HR |
| `29_meta_analysis.R` | 核心 · 主图主表 | Fig 4B: random-effects meta-analysis |
| `28_forest_plots.R` | 补充 · 主文引用的 Fig S | Fig S2: per-cohort HR forest plots |
| `33_table1_cohorts.R` | 核心 · 主图主表 | Table 1: cohort characteristics |
| `49_scissor.R` | 核心 · 主图主表 | Fig 5A: Scissor poor-prognosis macrophage enrichment |
| `50_viper.R` | 核心 · 主图主表 | Fig 5B: VIPER TF activity (HIF1A/EPAS1/RELA up, IFN down) |
| `64_regulon_enrichment.R` | 补充 · 主文引用的 Fig S | Fig S10 + Table S11: DoRothEA regulon enrichment (NF-kB master regulator) |
| `45_immune_scores.R` | 核心 · 主图主表 | Compute CYT / IFN-g / IPS immune scores |
| `46_immune_scores_plot.R` | 核心 · 主图主表 | Fig 6B: immune-score boxplots |
| `38_checkpoint_hla.R` | 核心 · 主图主表 | Fig 6A: checkpoint + HLA scores |
| `48_metabolism.R` | 核心 · 主图主表 | Fig 6C: glycolysis / hypoxia / OXPHOS scores |
| `47_tmb_cold_quadrant.R` | 核心 · 主图主表 | Fig 6D: TMB x cytolytic quadrant |
| `46_tmb_egfr_smoking.R` | 核心 · 主图主表 | Fig 6 + S4: TMB / EGFR / smoking stratifications |
| `51_macrophage_subtype_bulk.R` | 核心 · 主图主表 | Fig 6E + S6: macrophage-subtype bulk scores |
| `36_nomogram.R` | 核心 · 主图主表 | Fig 7A/B: nomogram + timeROC + calibration + DCA |
| `44_drug_sensitivity.R` | 核心 · 主图主表 | Fig 7C: oncoPredict IC50 by risk group |
| `31_icb_response.R` | 核心 · 主图主表 | Fig 7 + S3: ICB response (GSE243013) |
| `39_icb_external.R` | 核心 · 主图主表 | S3: ICB external cohorts (GSE135222 / GSE91061) |
| `62_cellchat_prep.R` | 核心 · 主图主表 | Fig 7D: CellChat object preparation |
| `63_cellchat_communication.R` | 核心 · 主图主表 | Fig 7D: differential cell-cell communication |
| `70_spp1_program_response.R` | 核心 · 主图主表 | Fig 7D: SPP1 / P1 program in non-responders |
| `69_regulon_ligand_check.R` | 核心 · 主图主表 | Regulon ligand check: SPP1 = AP-1 (JUN) target, CXCL8/IL1B = NF-kB targets |
| `37_ml_combination.R` | 核心 · 主图主表 | Fig 7E: LASSO vs RSF vs XGBoost comparison |
| `56_cindex_incremental.R` | 核心 · 主图主表 | Section 4.5: incremental C-index + likelihood-ratio test |
| `52_pancancer.R` | 核心 · 主图主表 | Fig 8: pan-cancer program HR across 28 cancer types |
| `57_pancancer_sign_test.R` | 核心 · 主图主表 | Fig 8: binomial sign test for direction consistency |
| `34_gene_level_lasso.R` | 补充 · 主文引用的 Fig S | Fig S7: gene-level LASSO (vs program-level) |
| `35_gene_level_bootstrap.R` | 补充 · 主文引用的 Fig S | Fig S7: gene-level bootstrap optimism correction |
| `41_nmf_rank_selection.R` | 补充 · 主文引用的 Fig S | Fig S8: NMF rank survey (cophenetic/silhouette) |
| `42_nmf_robustness_across_k.R` | 补充 · 主文引用的 Fig S | Fig S9: NMF re-run across k (cloud workflow) |
| `61_crossk_figures.R` | 补充 · 主文引用的 Fig S | Fig S9: cross-k + cross-algorithm overlap heatmaps |
| `74_spp1_prognosis.R` | 差异化 · SPP1/regulon 深挖 | SPP1 single-gene univariate Cox in 4 cohorts |
| `75_spp1_component_forest.R` | 差异化 · SPP1/regulon 深挖 | SPP1 vs C1q component forest (P1 decomposition) |
| `60_icb_reduced_signature.R` | 差异化 · SPP1/regulon 深挖 | Reduced adverse-axis (P8/P9/P4) ICB test (Section 4.8) |

Run order follows the list top-to-bottom (discovery -> bulk scoring -> signature -> validation -> mechanistic -> clinical -> pan-cancer).
