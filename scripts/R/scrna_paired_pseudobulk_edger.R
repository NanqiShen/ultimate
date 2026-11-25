#!/usr/bin/env Rscript

parse_args <- function(argv) {
  out <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key))
    }
    name <- sub("^--", "", key)
    if (i == length(argv) || startsWith(argv[[i + 1]], "--")) {
      out[[name]] <- TRUE
      i <- i + 1
    } else {
      out[[name]] <- argv[[i + 1]]
      i <- i + 2
    }
  }
  out
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) {
    stop(sprintf("Missing required argument --%s", name))
  }
  value
}

optional_arg <- function(args, name, default = NULL) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) {
    default
  } else {
    value
  }
}

safe_mkdir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

write_tsv <- function(x, path) {
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

package_version_or_na <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    as.character(utils::packageVersion(pkg))
  } else {
    NA_character_
  }
}

safe_id <- function(x) {
  gsub("[^A-Za-z0-9_.+-]+", "_", x)
}

validate_counts_matrix <- function(counts) {
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("counts_missing_dimnames")
  }
  if (anyNA(rownames(counts)) || anyNA(colnames(counts))) {
    stop("counts_na_dimnames")
  }
  if (anyDuplicated(rownames(counts))) {
    stop("counts_duplicate_feature_ids")
  }
  if (anyDuplicated(colnames(counts))) {
    stop("counts_duplicate_sample_ids")
  }
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  if (any(!is.finite(values))) {
    stop("counts_nonfinite")
  }
  if (any(values < 0)) {
    stop("counts_negative")
  }
  if (any(abs(values - round(values)) > 1e-8)) {
    stop("counts_noninteger")
  }
  invisible(TRUE)
}

read_sce_inputs <- function(sce_rda, object_name, assay_name, donor_column, condition_column, cell_type_column, multiplet_column, multiplet_keep) {
  suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(Matrix)
  })
  loaded <- load(sce_rda)
  if (!is.null(object_name) && nzchar(object_name)) {
    if (!object_name %in% loaded) {
      stop(sprintf("sce_object_missing:%s", object_name))
    }
    sce <- get(object_name)
  } else {
    if (length(loaded) != 1) {
      stop(sprintf("sce_object_ambiguous:%s", paste(loaded, collapse = ",")))
    }
    object_name <- loaded[[1]]
    sce <- get(object_name)
  }
  if (!is(sce, "SingleCellExperiment")) {
    stop(sprintf("not_singlecellexperiment:%s", class(sce)[[1]]))
  }
  if (!assay_name %in% assayNames(sce)) {
    stop(sprintf("sce_assay_missing:%s", assay_name))
  }
  counts <- assay(sce, assay_name)
  rownames(counts) <- rownames(sce)
  colnames(counts) <- colnames(sce)
  validate_counts_matrix(counts)
  cell_metadata <- as.data.frame(colData(sce), stringsAsFactors = FALSE)
  required_cols <- c(donor_column, condition_column, cell_type_column)
  if (!is.null(multiplet_column) && nzchar(multiplet_column) && !is.null(multiplet_keep) && nzchar(multiplet_keep)) {
    required_cols <- c(required_cols, multiplet_column)
  }
  missing_cols <- setdiff(required_cols, colnames(cell_metadata))
  if (length(missing_cols) > 0) {
    stop(sprintf("sce_coldata_missing:%s", paste(missing_cols, collapse = ",")))
  }
  cell_metadata$cell_id <- colnames(sce)
  cell_metadata$donor <- as.character(cell_metadata[[donor_column]])
  cell_metadata$condition <- as.character(cell_metadata[[condition_column]])
  cell_metadata$cell_type <- as.character(cell_metadata[[cell_type_column]])
  keep <- rep(TRUE, nrow(cell_metadata))
  metadata_filter <- ""
  if (!is.null(multiplet_column) && nzchar(multiplet_column) && !is.null(multiplet_keep) && nzchar(multiplet_keep)) {
    cell_metadata$multiplet_status <- as.character(cell_metadata[[multiplet_column]])
    keep <- keep & cell_metadata$multiplet_status == multiplet_keep
    metadata_filter <- paste(multiplet_column, multiplet_keep, sep = "==")
  }
  keep[is.na(keep)] <- FALSE
  if (!any(keep)) {
    stop("no_cells_after_metadata_filter")
  }
  list(
    object_name = object_name,
    assay_name = assay_name,
    counts = counts[, keep, drop = FALSE],
    cell_metadata = cell_metadata[keep, , drop = FALSE],
    feature_metadata = data.frame(feature_id = rownames(sce), stringsAsFactors = FALSE),
    metadata_filter = metadata_filter
  )
}

read_preaggregated_inputs <- function(counts_path, design_path, condition_column, donor_column, cell_type_column) {
  counts_frame <- utils::read.table(counts_path, sep = "\t", header = TRUE, check.names = FALSE, comment.char = "", quote = "")
  if (!"feature_id" %in% names(counts_frame)) {
    names(counts_frame)[[1]] <- "feature_id"
  }
  if (anyDuplicated(as.character(counts_frame$feature_id))) {
    stop("counts_duplicate_feature_ids")
  }
  counts <- as.matrix(counts_frame[, setdiff(names(counts_frame), "feature_id"), drop = FALSE])
  mode(counts) <- "numeric"
  rownames(counts) <- as.character(counts_frame$feature_id)
  design <- utils::read.table(design_path, sep = "\t", header = TRUE, check.names = FALSE, comment.char = "", quote = "")
  required_cols <- c("pseudobulk_id", "sample_id", donor_column, condition_column, cell_type_column)
  missing_cols <- setdiff(required_cols, names(design))
  if (length(missing_cols) > 0) {
    stop(sprintf("design_missing_columns:%s", paste(missing_cols, collapse = ",")))
  }
  design$pseudobulk_id <- as.character(design$pseudobulk_id)
  design$sample_id <- as.character(design$sample_id)
  design[[donor_column]] <- as.character(design[[donor_column]])
  design[[condition_column]] <- as.character(design[[condition_column]])
  design[[cell_type_column]] <- as.character(design[[cell_type_column]])
  if (!"cell_count" %in% names(design)) {
    design$cell_count <- NA_integer_
  }
  validate_counts_matrix(counts)
  missing_counts <- setdiff(design$pseudobulk_id, colnames(counts))
  if (length(missing_counts) > 0) {
    stop(sprintf("counts_missing_pseudobulk_ids:%s", paste(missing_counts, collapse = ",")))
  }
  list(counts = counts, design = design, feature_metadata = data.frame(feature_id = rownames(counts), stringsAsFactors = FALSE), cell_metadata = NULL)
}

aggregate_counts <- function(counts, cell_metadata, donor_column, condition_column, cell_type_column, min_cells) {
  suppressPackageStartupMessages({
    library(Matrix)
    library(data.table)
  })
  dt <- as.data.table(cell_metadata)
  dt[, donor_value := as.character(get(donor_column))]
  dt[, condition_value := as.character(get(condition_column))]
  dt[, cell_type_value := as.character(get(cell_type_column))]
  dt[, pseudobulk_id := safe_id(paste(donor_value, condition_value, cell_type_value, sep = "__"))]
  dt[, sample_id := safe_id(paste(donor_value, condition_value, sep = "__"))]
  group_rows <- dt[, .(cell_indices = list(.I), cell_count = .N), by = .(pseudobulk_id, sample_id, donor_value, condition_value, cell_type_value)]
  group_rows <- group_rows[cell_count >= min_cells]
  if (!nrow(group_rows)) {
    stop("no_pseudobulk_samples_after_min_cells")
  }
  groups <- group_rows$cell_indices
  names(groups) <- group_rows$pseudobulk_id
  aggregation <- lapply(groups, function(idx) Matrix::rowSums(counts[, idx, drop = FALSE]))
  pseudobulk_counts <- as.matrix(do.call(cbind, aggregation))
  rownames(pseudobulk_counts) <- rownames(counts)
  colnames(pseudobulk_counts) <- group_rows$pseudobulk_id
  design <- data.frame(
    pseudobulk_id = group_rows$pseudobulk_id,
    sample_id = group_rows$sample_id,
    donor = group_rows$donor_value,
    condition = group_rows$condition_value,
    cell_type = group_rows$cell_type_value,
    cell_count = group_rows$cell_count,
    stringsAsFactors = FALSE
  )
  list(counts = pseudobulk_counts, design = design)
}

condition_contrast <- function(model_formula, donor_values, condition_column, control, case) {
  ctrl_data <- data.frame(donor = factor(donor_values, levels = donor_values), condition = factor(rep(control, length(donor_values)), levels = c(control, case)), check.names = FALSE)
  case_data <- data.frame(donor = factor(donor_values, levels = donor_values), condition = factor(rep(case, length(donor_values)), levels = c(control, case)), check.names = FALSE)
  names(ctrl_data)[names(ctrl_data) == "condition"] <- condition_column
  names(case_data)[names(case_data) == "condition"] <- condition_column
  colMeans(model.matrix(model_formula, case_data) - model.matrix(model_formula, ctrl_data))
}

plot_empty <- function(path, title, message) {
  png(path, width = 1100, height = 780, res = 160)
  on.exit(dev.off(), add = TRUE)
  plot.new()
  title(main = title)
  text(0.5, 0.5, message)
}

plot_volcano <- function(results, path) {
  if (!nrow(results)) {
    plot_empty(path, "Pseudobulk edgeR volcano", "No testable genes")
    return(invisible(NULL))
  }
  png(path, width = 1100, height = 780, res = 160)
  on.exit(dev.off(), add = TRUE)
  padj <- results$padj
  padj[is.na(padj)] <- 1
  y <- -log10(pmax(padj, 1e-300))
  x <- results$log2FoldChange
  cls <- ifelse(padj < 0.05 & x > 0, "up", ifelse(padj < 0.05 & x < 0, "down", "ns"))
  colors <- c(up = "#B65F63", down = "#4E79A7", ns = "#96A0AA")
  plot(x, y, pch = 16, col = adjustcolor(colors[cls], alpha.f = 0.72), cex = 0.55,
       xlab = "log2 fold change (case - control)", ylab = "-log10 FDR", main = "Pseudobulk edgeR QL")
  abline(v = 0, lty = 2, col = "#7A8793")
  grid(col = "#E6E9ED")
}

plot_paired_genes <- function(log_cpm, results, design, donor_column, condition_column, control, case, path, n_genes = 6) {
  top <- head(results$feature_id[order(results$padj, results$pvalue, na.last = TRUE)], n_genes)
  if (!length(top)) {
    plot_empty(path, "Paired expression diagnostics", "No genes available")
    return(invisible(NULL))
  }
  png(path, width = 1200, height = 850, res = 160)
  on.exit(dev.off(), add = TRUE)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(2, ceiling(length(top) / 2)), mar = c(4, 4, 3, 1))
  for (gene in top) {
    values <- log_cpm[gene, design$pseudobulk_id]
    plot(c(1, 2), range(values), type = "n", xaxt = "n", xlab = "", ylab = "logCPM", main = gene)
    axis(1, at = c(1, 2), labels = c(control, case))
    for (donor_id in unique(design[[donor_column]])) {
      rows <- design[design[[donor_column]] == donor_id, , drop = FALSE]
      ctrl <- rows$pseudobulk_id[rows[[condition_column]] == control]
      stim <- rows$pseudobulk_id[rows[[condition_column]] == case]
      if (length(ctrl) == 1 && length(stim) == 1) {
        lines(c(1, 2), c(values[[ctrl]], values[[stim]]), col = "#778899")
        points(c(1, 2), c(values[[ctrl]], values[[stim]]), pch = 16, col = "#334E68")
      }
    }
  }
}

run_edgeR_ql <- function(counts, design, condition_column, donor_column, cell_type_column, control, case, target_cell_type, min_pairs, backend_id, analysis_level, tables_dir, figures_dir, robust) {
  suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
  })
  if (!is.null(target_cell_type) && nzchar(target_cell_type)) {
    design <- design[as.character(design[[cell_type_column]]) == target_cell_type, , drop = FALSE]
  }
  if (!nrow(design)) {
    stop(sprintf("target_cell_type_absent:%s", target_cell_type))
  }
  if (anyDuplicated(design$pseudobulk_id)) {
    stop("design_duplicate_pseudobulk_id")
  }
  missing_counts <- setdiff(design$pseudobulk_id, colnames(counts))
  if (length(missing_counts) > 0) {
    stop(sprintf("counts_missing_pseudobulk_ids:%s", paste(missing_counts, collapse = ",")))
  }
  counts <- counts[, design$pseudobulk_id, drop = FALSE]
  validate_counts_matrix(counts)
  design[[condition_column]] <- as.character(design[[condition_column]])
  design[[donor_column]] <- as.character(design[[donor_column]])
  if (!all(c(control, case) %in% unique(design[[condition_column]]))) {
    stop("need_control_and_case_conditions")
  }
  donor_condition <- as.data.frame(table(design[[donor_column]], design[[condition_column]]), stringsAsFactors = FALSE)
  names(donor_condition) <- c("donor", "condition", "n")
  complete_donors <- unique(donor_condition$donor[donor_condition$condition == control & donor_condition$n == 1])
  complete_donors <- intersect(complete_donors, unique(donor_condition$donor[donor_condition$condition == case & donor_condition$n == 1]))
  if (length(complete_donors) < min_pairs) {
    stop(sprintf("insufficient_complete_paired_donors:%s", length(complete_donors)))
  }
  design <- design[design[[donor_column]] %in% complete_donors, , drop = FALSE]
  design <- design[order(design[[donor_column]], design[[condition_column]]), , drop = FALSE]
  counts <- counts[, design$pseudobulk_id, drop = FALSE]
  design[[condition_column]] <- factor(design[[condition_column]], levels = c(control, case))
  design[[donor_column]] <- factor(design[[donor_column]])
  formula <- stats::as.formula(sprintf("~ %s + %s", donor_column, condition_column))
  design_matrix <- model.matrix(formula, data = design)
  if (qr(design_matrix)$rank < ncol(design_matrix)) {
    stop("design_matrix_rank_deficient")
  }
  residual_df <- nrow(design_matrix) - qr(design_matrix)$rank
  if (residual_df <= 0) {
    stop(sprintf("nonpositive_residual_df:%s", residual_df))
  }
  contrast <- condition_contrast(formula, levels(design[[donor_column]]), condition_column, control, case)
  if (any(!is.finite(contrast)) || length(contrast) != ncol(design_matrix) || all(abs(contrast) < 1e-12)) {
    stop("contrast_not_estimable")
  }
  y <- edgeR::DGEList(counts = round(counts), samples = design)
  keep <- edgeR::filterByExpr(y, design = design_matrix)
  filtered <- data.frame(feature_id = rownames(y), retained = keep, stringsAsFactors = FALSE)
  write_tsv(filtered, file.path(tables_dir, "pseudobulk_feature_filter.tsv"))
  y <- y[keep, , keep.lib.sizes = FALSE]
  if (nrow(y) < 2) {
    stop("too_few_features_after_filter")
  }
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design_matrix)
  fit <- edgeR::glmQLFit(y, design_matrix, robust = robust)
  qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
  tab <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table
  tab$feature_id <- rownames(tab)
  names(tab)[names(tab) == "logFC"] <- "log2FoldChange"
  names(tab)[names(tab) == "PValue"] <- "pvalue"
  names(tab)[names(tab) == "FDR"] <- "padj"
  tab$cell_type <- target_cell_type
  tab$backend_id <- backend_id
  tab$backend_method <- "edgeR"
  tab$comparison <- paste(case, "vs", control, sep = "_")
  tab$design <- sprintf("~ %s + %s", donor_column, condition_column)
  tab$analysis_level <- analysis_level
  tab$interpretation_warning <- "Pseudobulk DE is sample-level statistical evidence, not cell-level independent testing."
  tab <- tab[, c("feature_id", setdiff(colnames(tab), "feature_id"))]
  status <- data.frame(
    backend_id = backend_id,
    cell_type = target_cell_type,
    status = "ready",
    analysis_level = analysis_level,
    backend_method = "edgeR",
    comparison = paste(case, "vs", control, sep = "_"),
    n_features_input = nrow(counts),
    n_features_tested = nrow(tab),
    n_samples = ncol(counts),
    control_replicates = sum(design[[condition_column]] == control),
    case_replicates = sum(design[[condition_column]] == case),
    biological_replicates = length(unique(design[[donor_column]])),
    residual_df = residual_df,
    design_rank = qr(design_matrix)$rank,
    design_columns = ncol(design_matrix),
    reason = "",
    interpretation_warning = "edgeR QL results require raw counts, biological replicates, and design review.",
    stringsAsFactors = FALSE
  )
  design_out <- as.data.frame(design_matrix)
  design_out$pseudobulk_id <- rownames(design_matrix)
  design_out <- design_out[, c("pseudobulk_id", setdiff(colnames(design_out), "pseudobulk_id")), drop = FALSE]
  write_tsv(design_out, file.path(tables_dir, "pseudobulk_design_matrix.tsv"))
  write_tsv(data.frame(term = names(contrast), contrast = as.numeric(contrast)), file.path(tables_dir, "pseudobulk_contrast.tsv"))
  write_tsv(data.frame(pseudobulk_id = colnames(y), library_size = y$samples$lib.size, norm_factor = y$samples$norm.factors), file.path(tables_dir, "pseudobulk_library_norm.tsv"))
  plot_volcano(tab, file.path(figures_dir, "pseudobulk_de_volcano.png"))
  log_cpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)
  png(file.path(figures_dir, "pseudobulk_mds.png"), width = 1100, height = 780, res = 160)
  limma::plotMDS(log_cpm, labels = paste(design[[donor_column]], design[[condition_column]], sep = ":"), main = "Pseudobulk sample MDS")
  dev.off()
  plot_paired_genes(log_cpm, tab, design, donor_column, condition_column, control, case, file.path(figures_dir, "pseudobulk_paired_top_genes.png"))
  list(results = tab, status = status, y = y, fit = fit, qlf = qlf, design = design, design_matrix = design_matrix, contrast = contrast, residual_df = residual_df)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
tables_dir <- required_arg(args, "tables-dir")
backend_id <- optional_arg(args, "backend-id", "scrna.pseudobulk.edger_paired")
analysis_level <- optional_arg(args, "analysis-level", "validated_backend")
condition_column <- optional_arg(args, "condition-column", "condition")
donor_column <- optional_arg(args, "donor-column", "donor")
cell_type_column <- optional_arg(args, "cell-type-column", "cell_type")
control <- optional_arg(args, "control", "ctrl")
case <- optional_arg(args, "case", "stim")
target_cell_type <- optional_arg(args, "cell-type", NULL)
assay_name <- optional_arg(args, "assay", "counts")
min_cells <- as.integer(optional_arg(args, "min-cells", "1"))
min_pairs <- as.integer(optional_arg(args, "min-pairs", "4"))
robust <- tolower(optional_arg(args, "robust", "true")) %in% c("true", "1", "yes", "y")

safe_mkdir(tables_dir)
figures_dir <- normalizePath(file.path(tables_dir, "..", "figures"), mustWork = FALSE)
objects_dir <- normalizePath(file.path(tables_dir, "..", "..", "objects"), mustWork = FALSE)
safe_mkdir(figures_dir)
safe_mkdir(objects_dir)
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("dependency_missing:jsonlite")
}

versions <- data.frame(
  package = c("edgeR", "limma", "SingleCellExperiment", "SummarizedExperiment", "Matrix", "jsonlite"),
  version = vapply(c("edgeR", "limma", "SingleCellExperiment", "SummarizedExperiment", "Matrix", "jsonlite"), package_version_or_na, character(1)),
  stringsAsFactors = FALSE
)
write_tsv(versions, file.path(tables_dir, "pseudobulk_de_backend_versions.tsv"))
if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop("dependency_missing:edgeR")
}
if (!requireNamespace("limma", quietly = TRUE)) {
  stop("dependency_missing:limma")
}

input <- if (!is.null(args[["sce-rda"]])) {
  sce_input <- read_sce_inputs(
    sce_rda = required_arg(args, "sce-rda"),
    object_name = optional_arg(args, "sce-object", NULL),
    assay_name = assay_name,
    donor_column = donor_column,
    condition_column = condition_column,
    cell_type_column = cell_type_column,
    multiplet_column = optional_arg(args, "multiplet-column", NULL),
    multiplet_keep = optional_arg(args, "multiplet-keep", NULL)
  )
  aggregated <- aggregate_counts(sce_input$counts, sce_input$cell_metadata, "donor", "condition", "cell_type", min_cells)
  list(
    counts = aggregated$counts,
    design = aggregated$design,
    feature_metadata = sce_input$feature_metadata,
    cell_metadata = sce_input$cell_metadata,
    source_kind = "SingleCellExperiment_Rda",
    sce_object = sce_input$object_name,
    assay_name = sce_input$assay_name,
    metadata_filter = sce_input$metadata_filter
  )
} else {
  pre <- read_preaggregated_inputs(
    counts_path = required_arg(args, "counts"),
    design_path = required_arg(args, "design"),
    condition_column = condition_column,
    donor_column = donor_column,
    cell_type_column = cell_type_column
  )
  list(counts = pre$counts, design = pre$design, feature_metadata = pre$feature_metadata, cell_metadata = NULL, source_kind = "preaggregated_tsv", sce_object = "", assay_name = "", metadata_filter = "")
}

write_tsv(data.frame(feature_id = rownames(input$counts), input$counts, check.names = FALSE), file.path(tables_dir, "pseudobulk_counts.tsv"))
write_tsv(input$design, file.path(tables_dir, "pseudobulk_design.tsv"))
write_tsv(input$feature_metadata, file.path(tables_dir, "pseudobulk_feature_metadata.tsv"))
if (!is.null(input$cell_metadata)) {
  write_tsv(input$cell_metadata, file.path(tables_dir, "pseudobulk_cell_metadata.tsv"))
}
write_tsv(data.frame(pseudobulk_id = colnames(input$counts), library_size = as.numeric(colSums(input$counts))), file.path(tables_dir, "pseudobulk_raw_library_sizes.tsv"))
write_tsv(as.data.frame(table(input$design$donor, input$design$condition, input$design$cell_type), stringsAsFactors = FALSE), file.path(tables_dir, "pseudobulk_sample_cell_counts_long.tsv"))

backend <- run_edgeR_ql(
  counts = input$counts,
  design = input$design,
  condition_column = "condition",
  donor_column = "donor",
  cell_type_column = "cell_type",
  control = control,
  case = case,
  target_cell_type = target_cell_type,
  min_pairs = min_pairs,
  backend_id = backend_id,
  analysis_level = analysis_level,
  tables_dir = tables_dir,
  figures_dir = figures_dir,
  robust = robust
)

write_tsv(backend$results, file.path(tables_dir, "pseudobulk_de_results.tsv"))
write_tsv(backend$status, file.path(tables_dir, "pseudobulk_de_backend_status.tsv"))
saveRDS(list(results = backend$results, status = backend$status, y = backend$y, fit = backend$fit, qlf = backend$qlf, design = backend$design, design_matrix = backend$design_matrix, contrast = backend$contrast), file.path(objects_dir, "scrna_pseudobulk_de_backend.rds"))

manifest <- list(
  backend_id = backend_id,
  status = "ready",
  analysis_level = analysis_level,
  backend_method = "edgeR",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  source_kind = input$source_kind,
  sce_rda = optional_arg(args, "sce-rda", ""),
  sce_object = input$sce_object,
  assay_name = input$assay_name,
  metadata_filter = input$metadata_filter,
  condition_column = condition_column,
  donor_column = donor_column,
  cell_type_column = cell_type_column,
  control = control,
  case = case,
  comparison = paste(case, "vs", control, sep = "_"),
  logfc_direction = paste(case, "-", control),
  design = "~ donor + condition",
  min_cells = min_cells,
  min_pairs = min_pairs,
  robust = robust,
  target_cell_type = target_cell_type,
  n_features_input = nrow(input$counts),
  n_pseudobulk_samples_all_cell_types = ncol(input$counts),
  n_pseudobulk_samples_tested = ncol(backend$y$counts),
  n_features_tested = nrow(backend$results),
  biological_replicates = backend$status$biological_replicates[[1]],
  residual_df = backend$residual_df,
  design_rank = qr(backend$design_matrix)$rank,
  design_columns = ncol(backend$design_matrix),
  artifacts = list(
    pseudobulk_counts = file.path(tables_dir, "pseudobulk_counts.tsv"),
    pseudobulk_design = file.path(tables_dir, "pseudobulk_design.tsv"),
    pseudobulk_feature_metadata = file.path(tables_dir, "pseudobulk_feature_metadata.tsv"),
    pseudobulk_cell_metadata = file.path(tables_dir, "pseudobulk_cell_metadata.tsv"),
    pseudobulk_feature_filter = file.path(tables_dir, "pseudobulk_feature_filter.tsv"),
    pseudobulk_design_matrix = file.path(tables_dir, "pseudobulk_design_matrix.tsv"),
    pseudobulk_contrast = file.path(tables_dir, "pseudobulk_contrast.tsv"),
    pseudobulk_library_norm = file.path(tables_dir, "pseudobulk_library_norm.tsv"),
    pseudobulk_de_results = file.path(tables_dir, "pseudobulk_de_results.tsv"),
    pseudobulk_de_backend_status = file.path(tables_dir, "pseudobulk_de_backend_status.tsv"),
    pseudobulk_de_backend_versions = file.path(tables_dir, "pseudobulk_de_backend_versions.tsv"),
    pseudobulk_de_volcano = file.path(figures_dir, "pseudobulk_de_volcano.png"),
    pseudobulk_mds = file.path(figures_dir, "pseudobulk_mds.png"),
    pseudobulk_paired_top_genes = file.path(figures_dir, "pseudobulk_paired_top_genes.png"),
    rds = file.path(objects_dir, "scrna_pseudobulk_de_backend.rds")
  ),
  versions = as.list(stats::setNames(versions$version, versions$package)),
  interpretation_warning = "Pseudobulk edgeR output is a sample-level statistical result for the tested design, not a cell-level independent test or causal mechanism claim."
)
jsonlite::write_json(manifest, file.path(tables_dir, "pseudobulk_de_backend_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
