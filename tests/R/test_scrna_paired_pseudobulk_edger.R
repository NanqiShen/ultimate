#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(jsonlite)
})

cmd <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", cmd[grep("^--file=", cmd)][1])
repo <- normalizePath(file.path(dirname(file_arg), "..", ".."), mustWork = FALSE)
script <- file.path(repo, "scripts", "R", "scrna_paired_pseudobulk_edger.R")
if (!file.exists(script)) {
  stop(sprintf("backend_script_missing:%s", script))
}

run_cmd <- function(args, expect_success = TRUE) {
  status <- system2("Rscript", shQuote(c(script, args)), stdout = TRUE, stderr = TRUE)
  code <- attr(status, "status")
  if (is.null(code)) code <- 0
  if (expect_success && code != 0) {
    cat(paste(status, collapse = "\n"), "\n")
    stop(sprintf("expected_success_got_%s", code))
  }
  if (!expect_success && code == 0) {
    stop("expected_failure_got_success")
  }
  status
}

make_fixture <- function(path) {
  genes <- paste0("Gene", LETTERS[1:6])
  donors <- paste0("D", 1:4)
  conditions <- c("ctrl", "stim")
  cells <- list()
  metadata <- data.frame()
  mat <- matrix(0, nrow = length(genes), ncol = 0, dimnames = list(genes, character()))
  set.seed(1)
  for (donor in donors) {
    for (condition in conditions) {
      for (i in seq_len(4)) {
        cell_id <- paste(donor, condition, i, sep = "_")
        baseline <- c(20, 20, 50, 50, 80, 80)
        if (condition == "stim") baseline[1] <- 200
        counts <- rpois(length(genes), lambda = baseline)
        mat <- cbind(mat, counts)
        colnames(mat)[ncol(mat)] <- cell_id
        metadata <- rbind(metadata, data.frame(ind = donor, stim = condition, cell = "CD4 T cells", multiplets = "singlet", stringsAsFactors = FALSE))
      }
    }
  }
  sce <- SingleCellExperiment(assays = list(counts = Matrix(mat, sparse = TRUE)), colData = metadata)
  save(sce, file = path)
  invisible(path)
}

tmp <- tempfile("ultimate_pseudobulk_test_")
dir.create(tmp)
sce_path <- file.path(tmp, "fixture.Rda")
make_fixture(sce_path)
out <- file.path(tmp, "run")
run_cmd(c(
  "--sce-rda", sce_path,
  "--tables-dir", file.path(out, "results", "tables"),
  "--donor-column", "ind",
  "--condition-column", "stim",
  "--cell-type-column", "cell",
  "--cell-type", "CD4 T cells",
  "--control", "ctrl",
  "--case", "stim",
  "--multiplet-column", "multiplets",
  "--multiplet-keep", "singlet",
  "--min-cells", "1",
  "--min-pairs", "4"
))
manifest <- fromJSON(file.path(out, "results", "tables", "pseudobulk_de_backend_manifest.json"))
stopifnot(identical(manifest$status, "ready"))
stopifnot(identical(manifest$biological_replicates, 4L))
result <- read.delim(file.path(out, "results", "tables", "pseudobulk_de_results.tsv"), check.names = FALSE)
gene_a <- result[result$feature_id == "GeneA", ]
stopifnot(nrow(gene_a) == 1)
stopifnot(gene_a$log2FoldChange > 0)
pb <- read.delim(file.path(out, "results", "tables", "pseudobulk_counts.tsv"), check.names = FALSE)
stopifnot("D1__stim__CD4_T_cells" %in% colnames(pb))
env <- new.env(parent = emptyenv())
load(sce_path, envir = env)
manual <- rowSums(assay(env$sce, "counts")[, colData(env$sce)$ind == "D1" & colData(env$sce)$stim == "stim", drop = FALSE])
stopifnot(all(pb$D1__stim__CD4_T_cells == as.numeric(manual[match(pb$feature_id, names(manual))])))

bad_counts <- file.path(tmp, "bad_counts.tsv")
bad_design <- file.path(tmp, "bad_design.tsv")
write.table(data.frame(feature_id = c("GeneA", "GeneB"), S1 = c(1.5, 2), S2 = c(1, 2)), bad_counts, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(pseudobulk_id = c("S1", "S2"), sample_id = c("S1", "S2"), donor = c("D1", "D1"), condition = c("ctrl", "stim"), cell_type = c("CD4 T cells", "CD4 T cells"), cell_count = c(10, 10)), bad_design, sep = "\t", row.names = FALSE, quote = FALSE)
run_cmd(c("--counts", bad_counts, "--design", bad_design, "--tables-dir", file.path(tmp, "bad", "results", "tables"), "--cell-type", "CD4 T cells", "--control", "ctrl", "--case", "stim"), expect_success = FALSE)

dup_design <- file.path(tmp, "dup_design.tsv")
write.table(data.frame(pseudobulk_id = c("S1", "S1"), sample_id = c("S1", "S1dup"), donor = c("D1", "D1"), condition = c("ctrl", "stim"), cell_type = c("CD4 T cells", "CD4 T cells"), cell_count = c(10, 10)), dup_design, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(feature_id = c("GeneA", "GeneB"), S1 = c(1, 2)), bad_counts, sep = "\t", row.names = FALSE, quote = FALSE)
run_cmd(c("--counts", bad_counts, "--design", dup_design, "--tables-dir", file.path(tmp, "dup", "results", "tables"), "--cell-type", "CD4 T cells", "--control", "ctrl", "--case", "stim"), expect_success = FALSE)

single_condition <- file.path(tmp, "single_condition.tsv")
write.table(data.frame(pseudobulk_id = c("S1", "S2"), sample_id = c("S1", "S2"), donor = c("D1", "D2"), condition = c("ctrl", "ctrl"), cell_type = c("CD4 T cells", "CD4 T cells"), cell_count = c(10, 10)), single_condition, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(feature_id = c("GeneA", "GeneB"), S1 = c(1, 2), S2 = c(1, 2)), bad_counts, sep = "\t", row.names = FALSE, quote = FALSE)
run_cmd(c("--counts", bad_counts, "--design", single_condition, "--tables-dir", file.path(tmp, "single", "results", "tables"), "--cell-type", "CD4 T cells", "--control", "ctrl", "--case", "stim"), expect_success = FALSE)

cat("scrna_paired_pseudobulk_edger_regression_pass\n")
