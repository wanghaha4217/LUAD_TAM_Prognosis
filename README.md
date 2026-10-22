# Code for: A macrophage functional-program prognostic signature in lung adenocarcinoma

R scripts that reproduce the single-cell atlas, NMF program decomposition, LASSO-Cox signature, validation, mechanistic (Scissor / VIPER / regulon), clinical (nomogram / drug / ICB / CellChat) and pan-cancer analyses.

## Repository layout

- `R/` — analysis scripts (run top-to-bottom; see `MANIFEST.md` for the figure/table each script produces).
- `MANIFEST.md` — script-to-figure/table mapping.

## How to run

1. Open the repository root in RStudio (File > Open Project, or set the working directory to the repo root).
2. Each script sets its base path as `"."` (the project root) and builds all downstream paths from it with `file.path()`.
   If you run a script from a different directory, replace `"."` in that script with the absolute path to the repo root.
3. Place the input data where the scripts expect it:
   - Discovery scRNA-seq (GEO GSE243013): `GSE243013_LUAD_TAM_clustered.rds` (and the `GSE243013/` raw folder) at the repo root.
   - Validation cohorts (GSE68465 / GSE72094 / GSE31210) under `validation_cohorts/`.
   - TCGA-LUAD expression / survival / mutation files at the repo root (see `00_fetch_validation_cohorts.R`).
4. Run the scripts in the order given in `MANIFEST.md`; later scripts read the objects saved by earlier ones.

## Data availability

Single-cell data: GEO GSE243013. Bulk validation: TCGA-LUAD and GEO (GSE68465 / GSE72094 / GSE31210). Mutation annotation: cBioPortal PanCancer Atlas. Pan-cancer expression: Xena PanCanAtlas.
