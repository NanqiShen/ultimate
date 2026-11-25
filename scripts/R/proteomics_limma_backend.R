#!/usr/bin/env Rscript
# Real limma moderated t-test for proteomics abundance tables.
# Called by ultimate.bulk when a limma-capable R environment is available.
# Writes a results table and a saveRDS() object (never JSON).

parse_args <- function(argv) {
  out <- list(); i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key))
    name <- sub("^--", "", key)
    if (i == length(argv) || startsWith(argv[[i + 1]], "--")) {
      out[[name]] <- TRUE; i <- i + 1
    } else {
      out[[name]] <- argv[[i + 1]]; i <- i + 2
    }
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!requireNamespace("limma", quietly = TRUE)) {
  cat('{"status":"skipped","reason":"dependency_missing:limma"}\n'); quit(status = 0)
}

matrix_path <- args[["matrix"]]
samples_path <- args[["samples"]]
output_path <- args[["output"]]
object_path <- args[["object"]]
condition <- if (is.null(args[["condition-column"]])) "condition" else args[["condition-column"]]
control <- if (is.null(args[["control"]])) "control" else args[["control"]]
case <- if (is.null(args[["case"]])) "treated" else args[["case"]]

for (path in c(matrix_path, samples_path)) {
  if (is.null(path) || !file.exists(path)) {
    cat(sprintf('{"status":"skipped","reason":"missing_input:%s"}\n', path)); quit(status = 0)
  }
}
mat <- as.matrix(read.delim(matrix_path, row.names = 1, check.names = FALSE))
samples <- read.delim(samples_path, stringsAsFactors = FALSE)
samples$sample_id <- as.character(samples$sample_id)
if (!condition %in% colnames(samples)) {
  cat(sprintf('{"status":"skipped","reason":"missing_condition_column:%s"}\n', condition)); quit(status = 0)
}
ids <- colnames(mat)
keep <- ids %in% samples$sample_id
mat <- mat[, keep, drop = FALSE]
grp <- samples$condition[match(colnames(mat), samples$sample_id)]
control_cols <- colnames(mat)[grp == control]
case_cols <- colnames(mat)[grp == case]
if (length(control_cols) < 2 || length(case_cols) < 2) {
  cat('{"status":"skipped","reason":"insufficient_biological_replicates"}\n'); quit(status = 0)
}

library(limma)
design <- model.matrix(~ 0 + factor(grp, levels = c(control, case)))
colnames(design) <- c(control, case)
fit <- lmFit(mat, design)
contrast <- makeContrasts(paste0(case, "-", control), levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast))
res <- topTable(fit2, coef = 1, number = nrow(mat), sort.by = "none")
res$feature_id <- rownames(mat)
res$log2FC <- res$logFC
res$padj <- res$adj.P.Val
res$method <- "limma::lmFit+eBayes"
out <- res[, c("feature_id", "logFC", "log2FC", "AveExpr", "t", "P.Value", "adj.P.Val", "padj", "method")]
utils::write.table(out, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(results = out, fit = fit2, design = design), object_path)
cat(sprintf('{"status":"ready","engine":"R_limma","n_features":%d}\n', nrow(out)))
