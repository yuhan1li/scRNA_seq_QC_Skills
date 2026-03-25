---
name: single-cell-basic-analysis
description: Standard Seurat-based single-cell RNA-seq analysis for raw count matrices or 10x folders. Use when Codex needs to build or adapt a reusable workflow for QC filtering, normalization, PCA/UMAP, clustering, Harmony batch correction, SingleR cell-type annotation, marker discovery, and cell composition bar plots across samples, groups, or conditions.
---

# Single-Cell Basic Analysis

Use `scripts/run_scrna_basic_pipeline.R` as the default implementation for baseline scRNA-seq analysis.

## Quick start

- Expect a sample sheet TSV with at least `sample_id` and `data_dir`.
- Prefer adding `group` and `batch` columns. If they are missing, the script falls back to `sample_id`.
- Read [references/input-format.md](references/input-format.md) only when the input table needs to be created or repaired.
- Run the pipeline with `Rscript scripts/run_scrna_basic_pipeline.R --sample-sheet /abs/path/sample_sheet.tsv --outdir /abs/path/results`.
- Use `--figure-dir` and `--bundle-dir` when figures must be separated from data and code outputs.

## What the script does

1. Read one or more 10x count folders from a sample sheet.
2. Build and merge Seurat objects with consistent sample metadata.
3. Compute mitochondrial ratio and run QC filtering.
4. Run normalization, variable-feature selection, scaling, PCA, neighbors, clustering, and UMAP.
5. Run Harmony when the batch column contains more than one batch.
6. Run cluster-level SingleR annotation with a celldex reference.
7. Export marker tables, cluster heatmap, UMAP plots, metadata tables, and cell composition bar plots.

## Default command

```bash
Rscript scripts/run_scrna_basic_pipeline.R \
  --sample-sheet /abs/path/sample_sheet.tsv \
  --outdir /abs/path/scrna_results \
  --figure-dir /abs/path/scrna_results/figures \
  --bundle-dir /abs/path/scrna_results/data_and_code \
  --species human \
  --reference hpca \
  --dims 1:30 \
  --resolution 0.6 \
  --min-features 300 \
  --max-features 7500 \
  --max-mt 20 \
  --batch-column batch \
  --group-column group
```

## Practical guidance

- Use `hpca` first for broad human cell classes.
- Use `blueprint` or `monaco` when the dataset is immune-heavy and finer immune annotation is needed.
- Keep `group` as the biological comparison used for composition plots.
- Keep `batch` as the technical batch used for Harmony correction.
- Review `annotation/singler_cluster_annotations.tsv` before trusting labels blindly. SingleR is a baseline, not the final biological truth.
- If the user already has a Seurat object rather than raw 10x folders, adapt the script instead of rewriting the whole workflow from scratch.

## Expected outputs

- `rds/seurat_integrated.rds`
- `metadata/cell_metadata.tsv`
- `figures/*.pdf`
- `markers/all_markers.tsv`
- `annotation/singler_cluster_annotations.tsv`
- `composition/celltype_counts.tsv`
- `composition/celltype_proportions.tsv`
