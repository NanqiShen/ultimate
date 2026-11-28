#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(fgsea)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) {
    return(args[idx + 1L])
  }
  default
}

ranked_path <- get_arg("--ranked")
gmt_path <- get_arg("--gmt")
output_path <- get_arg("--output")
seed_val <- as.integer(get_arg("--seed", "123"))

if (is.null(ranked_path) || is.null(gmt_path) || is.null(output_path)) {
  stop("Usage: enrichment_fgsea_backend.R --ranked <tsv> --gmt <gmt> --output <tsv> [--seed 123]")
}

ranked_df <- utils::read.delim(ranked_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("feature_id", "score") %in% colnames(ranked_df))) {
  stop("Ranked table must contain feature_id and score columns")
}
ranked_df$score <- suppressWarnings(as.numeric(ranked_df$score))
ranked_df <- ranked_df[!is.na(ranked_df$feature_id) & nzchar(ranked_df$feature_id) & is.finite(ranked_df$score), , drop = FALSE]
ranked_df <- ranked_df[!duplicated(ranked_df$feature_id), , drop = FALSE]
stats <- stats::setNames(ranked_df$score, as.character(ranked_df$feature_id))
stats <- sort(stats, decreasing = TRUE)

pathways <- fgsea::gmtPathways(gmt_path)
if (length(pathways) == 0L || length(stats) < 5L) {
  empty_df <- data.frame(
    Term = character(0),
    ES = numeric(0),
    NES = numeric(0),
    `NOM p-val` = numeric(0),
    `FDR q-val` = numeric(0),
    Lead_genes = character(0),
    size = integer(0),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  utils::write.table(empty_df, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
  quit(save = "no", status = 0L)
}

set.seed(seed_val)
res <- suppressWarnings(fgsea::fgsea(pathways = pathways, stats = stats, minSize = 5L, maxSize = 500L))
if (nrow(res) > 0L) {
  res <- res[order(res$padj, res$pval, -abs(res$NES)), , drop = FALSE]
  lead <- vapply(res$leadingEdge, function(x) paste(as.character(x), collapse = ";"), character(1))
  terms <- gsub("[\t\r\n]+", " ", as.character(res$pathway))
  out_df <- data.frame(
    Term = terms,
    ES = as.numeric(res$ES),
    NES = as.numeric(res$NES),
    `NOM p-val` = as.numeric(res$pval),
    `FDR q-val` = as.numeric(res$padj),
    Lead_genes = lead,
    size = as.integer(res$size),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
} else {
  out_df <- data.frame(
    Term = character(0),
    ES = numeric(0),
    NES = numeric(0),
    `NOM p-val` = numeric(0),
    `FDR q-val` = numeric(0),
    Lead_genes = character(0),
    size = integer(0),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
utils::write.table(out_df, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)
