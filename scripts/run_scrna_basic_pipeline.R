#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(harmony)
  library(SingleR)
  library(celldex)
})

parse_args <- function(args) {
  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      parsed[[key]] <- TRUE
      i <- i + 1
    } else {
      parsed[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  parsed
}

get_arg <- function(args, key, default = NULL, required = FALSE) {
  if (!is.null(args[[key]])) {
    return(args[[key]])
  }
  if (required) {
    stop("Missing required argument --", key, call. = FALSE)
  }
  default
}

parse_numeric <- function(x) {
  if (is.null(x)) {
    return(NA_real_)
  }
  if (x %in% c("Inf", "inf", "INF")) {
    return(Inf)
  }
  as.numeric(x)
}

parse_dims <- function(x) {
  if (grepl(":", x, fixed = TRUE)) {
    parts <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    if (length(parts) != 2 || any(is.na(parts))) {
      stop("Invalid --dims value: ", x, call. = FALSE)
    }
    return(seq.int(parts[1], parts[2]))
  }
  dims <- as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
  if (any(is.na(dims))) {
    stop("Invalid --dims value: ", x, call. = FALSE)
  }
  dims
}

dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

save_plot_pdf <- function(filepath, plot_obj, width = 8, height = 6) {
  pdf(filepath, width = width, height = height)
  print(plot_obj)
  dev.off()
}

load_reference <- function(reference_name, species) {
  ref_name <- gsub("-", "_", tolower(reference_name), fixed = TRUE)
  if (species == "human") {
    ref <- switch(
      ref_name,
      hpca = HumanPrimaryCellAtlasData(),
      blueprint = BlueprintEncodeData(),
      monaco = MonacoImmuneData(),
      dice = DatabaseImmuneCellExpressionData(),
      immgen = stop("immgen is mouse-only. Choose hpca, blueprint, monaco, or dice.", call. = FALSE),
      mouse_rna = stop("mouse_rna is mouse-only. Choose a human reference.", call. = FALSE),
      stop("Unsupported human reference: ", reference_name, call. = FALSE)
    )
  } else if (species == "mouse") {
    ref <- switch(
      ref_name,
      immgen = ImmGenData(),
      mouse_rna = MouseRNAseqData(),
      hpca = stop("hpca is human-only. Choose immgen or mouse_rna.", call. = FALSE),
      blueprint = stop("blueprint is human-only. Choose immgen or mouse_rna.", call. = FALSE),
      monaco = stop("monaco is human-only. Choose immgen or mouse_rna.", call. = FALSE),
      dice = stop("dice is human-only. Choose immgen or mouse_rna.", call. = FALSE),
      stop("Unsupported mouse reference: ", reference_name, call. = FALSE)
    )
  } else {
    stop("Unsupported species: ", species, call. = FALSE)
  }

  label_column <- if ("label.main" %in% colnames(colData(ref))) {
    "label.main"
  } else if ("label.fine" %in% colnames(colData(ref))) {
    "label.fine"
  } else {
    stop("Reference object does not contain label.main or label.fine.", call. = FALSE)
  }

  list(ref = ref, labels = colData(ref)[[label_column]])
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

sample_sheet_path <- get_arg(args, "sample-sheet", required = TRUE)
outdir <- normalizePath(get_arg(args, "outdir", required = TRUE), winslash = "/", mustWork = FALSE)
figure_dir <- normalizePath(get_arg(args, "figure-dir", file.path(outdir, "figures")), winslash = "/", mustWork = FALSE)
bundle_dir <- normalizePath(get_arg(args, "bundle-dir", file.path(outdir, "data_and_code")), winslash = "/", mustWork = FALSE)
project_name <- get_arg(args, "project", "scrna")
species <- tolower(get_arg(args, "species", "human"))
reference_name <- get_arg(args, "reference", ifelse(species == "human", "hpca", "immgen"))
dims <- parse_dims(get_arg(args, "dims", "1:30"))
resolution <- parse_numeric(get_arg(args, "resolution", "0.6"))
min_cells <- as.integer(parse_numeric(get_arg(args, "min-cells", "3")))
min_features_create <- as.integer(parse_numeric(get_arg(args, "min-features-create", "200")))
min_features_filter <- parse_numeric(get_arg(args, "min-features", "300"))
max_features_filter <- parse_numeric(get_arg(args, "max-features", "7500"))
max_mt <- parse_numeric(get_arg(args, "max-mt", "20"))
max_count <- parse_numeric(get_arg(args, "max-count", "Inf"))
nfeatures <- as.integer(parse_numeric(get_arg(args, "nfeatures", "3000")))
batch_column <- get_arg(args, "batch-column", "batch")
group_column <- get_arg(args, "group-column", "group")
seed <- as.integer(parse_numeric(get_arg(args, "seed", "1234")))

set.seed(seed)

sample_sheet <- read.delim(sample_sheet_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("sample_id", "data_dir")
missing_columns <- setdiff(required_columns, colnames(sample_sheet))
if (length(missing_columns) > 0) {
  stop("Sample sheet is missing required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

if (!"group" %in% colnames(sample_sheet)) {
  sample_sheet$group <- sample_sheet$sample_id
}
if (!"batch" %in% colnames(sample_sheet)) {
  sample_sheet$batch <- sample_sheet$sample_id
}
if (!group_column %in% colnames(sample_sheet)) {
  sample_sheet[[group_column]] <- sample_sheet$group
}
if (!batch_column %in% colnames(sample_sheet)) {
  sample_sheet[[batch_column]] <- sample_sheet$batch
}

sample_sheet$data_dir <- normalizePath(sample_sheet$data_dir, winslash = "/", mustWork = TRUE)

plot_dir <- figure_dir
marker_dir <- file.path(bundle_dir, "markers")
annotation_dir <- file.path(bundle_dir, "annotation")
composition_dir <- file.path(bundle_dir, "composition")
metadata_dir <- file.path(bundle_dir, "metadata")
rds_dir <- file.path(bundle_dir, "rds")
code_dir <- file.path(bundle_dir, "code")

invisible(lapply(
  list(outdir, plot_dir, bundle_dir, marker_dir, annotation_dir, composition_dir, metadata_dir, rds_dir, code_dir),
  dir_create
))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)
  script_target <- file.path(code_dir, basename(script_path))
  if (file.exists(script_path) && normalizePath(script_path, winslash = "/", mustWork = FALSE) != normalizePath(script_target, winslash = "/", mustWork = FALSE)) {
    file.copy(script_path, file.path(code_dir, basename(script_path)), overwrite = TRUE)
  }
}
if (file.exists(sample_sheet_path)) {
  sample_target <- file.path(code_dir, basename(sample_sheet_path))
  if (normalizePath(sample_sheet_path, winslash = "/", mustWork = FALSE) != normalizePath(sample_target, winslash = "/", mustWork = FALSE)) {
    file.copy(sample_sheet_path, sample_target, overwrite = TRUE)
  }
}

message("Reading samples...")
object_list <- lapply(seq_len(nrow(sample_sheet)), function(i) {
  counts <- Read10X(data.dir = sample_sheet$data_dir[i])
  obj <- CreateSeuratObject(
    counts = counts,
    project = project_name,
    min.cells = min_cells,
    min.features = min_features_create
  )
  metadata_columns <- setdiff(colnames(sample_sheet), "data_dir")
  for (column_name in metadata_columns) {
    obj[[column_name]] <- sample_sheet[[column_name]][i]
  }
  obj
})

if (length(object_list) == 1) {
  seu <- object_list[[1]]
} else {
  seu <- merge(
    x = object_list[[1]],
    y = object_list[-1],
    add.cell.ids = sample_sheet$sample_id,
    project = project_name
  )
}

mt_pattern <- if (species == "mouse") "^mt-" else "^MT-"
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = mt_pattern)

save_plot_pdf(
  file.path(plot_dir, "qc_violin_before_filtering.pdf"),
  VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1),
  width = 12,
  height = 4
)

qc_keep <- (
  seu$nFeature_RNA >= min_features_filter &
    seu$nFeature_RNA <= max_features_filter &
    seu$percent.mt <= max_mt &
    seu$nCount_RNA <= max_count
)
seu <- subset(seu, cells = colnames(seu)[qc_keep])
if (ncol(seu) == 0) {
  stop("All cells were filtered out. Relax the QC thresholds and rerun.", call. = FALSE)
}

write.table(
  as.data.frame(table(seu$sample_id), stringsAsFactors = FALSE),
  file = file.path(metadata_dir, "cells_per_sample_after_filtering.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

save_plot_pdf(
  file.path(plot_dir, "qc_violin_after_filtering.pdf"),
  VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1),
  width = 12,
  height = 4
)

message("Running Seurat preprocessing...")
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = nfeatures)
seu <- ScaleData(seu, vars.to.regress = "percent.mt")
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = max(dims))

save_plot_pdf(
  file.path(plot_dir, "elbow_plot.pdf"),
  ElbowPlot(seu, ndims = max(dims)),
  width = 8,
  height = 5
)

use_harmony <- batch_column %in% colnames(seu@meta.data) &&
  length(unique(seu@meta.data[[batch_column]])) > 1

reduction_name <- "pca"
if (use_harmony) {
  message("Running Harmony batch correction on column: ", batch_column)
  seu <- RunHarmony(
    object = seu,
    group.by.vars = batch_column,
    reduction.use = "pca",
    dims.use = dims,
    verbose = FALSE
  )
  reduction_name <- "harmony"
}

seu <- RunUMAP(seu, reduction = reduction_name, dims = dims)
seu <- FindNeighbors(seu, reduction = reduction_name, dims = dims)
seu <- FindClusters(seu, resolution = resolution)

save_plot_pdf(
  file.path(plot_dir, "umap_by_cluster.pdf"),
  DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + NoLegend(),
  width = 8,
  height = 6
)

if (group_column %in% colnames(seu@meta.data)) {
  save_plot_pdf(
    file.path(plot_dir, "umap_by_group.pdf"),
    DimPlot(seu, reduction = "umap", group.by = group_column),
    width = 8,
    height = 6
  )
}

if (batch_column %in% colnames(seu@meta.data)) {
  save_plot_pdf(
    file.path(plot_dir, "umap_by_batch.pdf"),
    DimPlot(seu, reduction = "umap", group.by = batch_column),
    width = 8,
    height = 6
  )
}

message("Running SingleR...")
reference_bundle <- load_reference(reference_name, species)
expr_mat <- GetAssayData(seu, slot = "data")
cluster_ids <- seu$seurat_clusters

singler_pred <- SingleR(
  test = expr_mat,
  ref = reference_bundle$ref,
  labels = reference_bundle$labels,
  method = "cluster",
  clusters = cluster_ids,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts"
)

annotation_map <- data.frame(
  seurat_cluster = rownames(singler_pred),
  singler_label = singler_pred$labels,
  singler_pruned = singler_pred$pruned.labels,
  stringsAsFactors = FALSE
)
annotation_map$celltype <- ifelse(
  is.na(annotation_map$singler_pruned) | annotation_map$singler_pruned == "",
  annotation_map$singler_label,
  annotation_map$singler_pruned
)

seu$celltype <- annotation_map$celltype[match(seu$seurat_clusters, annotation_map$seurat_cluster)]

write.table(
  annotation_map,
  file = file.path(annotation_dir, "singler_cluster_annotations.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

save_plot_pdf(
  file.path(plot_dir, "umap_by_celltype.pdf"),
  DimPlot(seu, reduction = "umap", group.by = "celltype", label = TRUE),
  width = 10,
  height = 7
)

message("Finding markers...")
markers <- FindAllMarkers(
  object = seu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.table(
  markers,
  file = file.path(marker_dir, "all_markers.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (nrow(markers) > 0) {
  top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
    ungroup()

  write.table(
    top_markers,
    file = file.path(marker_dir, "top10_markers_per_cluster.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  save_plot_pdf(
    file.path(plot_dir, "marker_heatmap_top10.pdf"),
    DoHeatmap(seu, features = unique(top_markers$gene), size = 3),
    width = 12,
    height = 10
  )
}

message("Building composition plots...")
meta_df <- seu@meta.data

if (!group_column %in% colnames(meta_df)) {
  stop("Metadata column for composition plotting does not exist: ", group_column, call. = FALSE)
}

composition_counts <- as.data.frame(
  table(meta_df[[group_column]], meta_df$celltype),
  stringsAsFactors = FALSE
)
colnames(composition_counts) <- c("group", "celltype", "cell_count")
composition_counts <- composition_counts[composition_counts$cell_count > 0, , drop = FALSE]

composition_props <- composition_counts %>%
  group_by(group) %>%
  mutate(proportion = cell_count / sum(cell_count)) %>%
  ungroup()

write.table(
  composition_counts,
  file = file.path(composition_dir, "celltype_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  composition_props,
  file = file.path(composition_dir, "celltype_proportions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

count_plot <- ggplot(composition_counts, aes(x = celltype, y = cell_count, fill = group)) +
  geom_col(position = position_dodge(width = 0.8)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(x = "Cell type", y = "Cell count", fill = group_column)

prop_plot <- ggplot(composition_props, aes(x = group, y = proportion, fill = celltype)) +
  geom_col(position = "fill") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = group_column, y = "Cell type proportion", fill = "Cell type")

save_plot_pdf(
  file.path(plot_dir, "celltype_counts_barplot.pdf"),
  count_plot,
  width = 11,
  height = 6
)
save_plot_pdf(
  file.path(plot_dir, "celltype_proportions_stacked_barplot.pdf"),
  prop_plot,
  width = 9,
  height = 6
)

write.table(
  meta_df,
  file = file.path(metadata_dir, "cell_metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

saveRDS(seu, file = file.path(rds_dir, "seurat_integrated.rds"))
writeLines(capture.output(sessionInfo()), con = file.path(bundle_dir, "sessionInfo.txt"))

message("Analysis completed. Results written to: ", outdir)
