# Input Format

Create a tab-delimited sample sheet with one row per sample.

## Required columns

- `sample_id`: Short unique sample name used in merged cell barcodes.
- `data_dir`: Absolute path to the 10x matrix directory.

## Recommended columns

- `group`: Biological group for downstream comparison and composition plots, for example `Tumor` and `Normal`.
- `batch`: Technical batch for Harmony correction, for example `Run1`, `Run2`, or `PatientA`.

If `group` or `batch` is missing, the pipeline uses `sample_id`.

## Example

```tsv
sample_id	group	batch	data_dir
P01	Primary	Run1	/data/project/P01/filtered_feature_bc_matrix
P02	Primary	Run2	/data/project/P02/filtered_feature_bc_matrix
M01	Metastatic	Run1	/data/project/M01/filtered_feature_bc_matrix
M02	Metastatic	Run2	/data/project/M02/filtered_feature_bc_matrix
```

## Notes

- Point `data_dir` directly at the folder containing `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`.
- For human data, mitochondrial genes are detected with `^MT-`; for mouse data, the script uses `^mt-`.
- Keep `group` and `batch` conceptually separate. Composition plots use `group`; Harmony uses `batch`.

