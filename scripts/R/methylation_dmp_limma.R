#!/usr/bin/env Rscript
# Real limma M-value differential methylation (DMP) backend.
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
beta <- as.matrix(read.delim(matrix_path, row.names = 1, check.names = FALSE))
samples <- read.delim(samples_path, stringsAsFactors = FALSE)
samples$sample_id <- as.character(samples$sample_id)
if (!condition %in% colnames(samples)) {
  cat(sprintf('{"status":"skipped","reason":"missing_condition_column:%s"}\n', condition)); quit(status = 0)
}
keep <- colnames(beta) %in% samples$sample_id
beta <- beta[, keep, drop = FALSE]
grp <- samples$condition[match(colnames(beta), samples$sample_id)]
control_cols <- colnames(beta)[grp == control]
case_cols <- colnames(beta)[grp == case]
if (length(control_cols) < 2 || length(case_cols) < 2) {
  cat('{"status":"skipped","reason":"insufficient_biological_replicates"}\n'); quit(status = 0)
}
# M-values are the recommended scale for limma on methylation arrays.
beta <- pmin(pmax(beta, 1e-4), 1 - 1e-4)
mvalues <- log2(beta / (1 - beta))

library(limma)
design <- model.matrix(~ 0 + factor(grp, levels = c(control, case)))
colnames(design) <- c(control, case)
fit <- lmFit(mvalues, design)
contrast <- makeContrasts(paste0(case, "-", control), levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast))
res <- topTable(fit2, coef = 1, number = nrow(mvalues), sort.by = "none")
res$region_id <- rownames(mvalues)
res$m_value_delta <- res$logFC
res$beta_delta <- rowMeans(beta[, case_cols, drop = FALSE]) - rowMeans(beta[, control_cols, drop = FALSE])
res$padj <- res$adj.P.Val
res$method <- "limma::Mvalue_lmFit+eBayes"
res$annotation_status <- "annotation_handoff_required_for_genomic_context"
out <- res[, c("region_id", "logFC", "m_value_delta", "beta_delta", "t", "P.Value", "adj.P.Val", "padj", "method", "annotation_status")]
utils::write.table(out, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(results = out, fit = fit2, design = design, beta = beta), object_path)
cat(sprintf('{"status":"ready","engine":"R_limma_mvalue","n_regions":%d}\n', nrow(out)))
