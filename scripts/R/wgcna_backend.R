#!/usr/bin/env Rscript
# Real WGCNA co-expression network: soft-threshold selection, blockwiseModules,
# module-trait correlations and hub genes. Never a quartile binning proxy.
# Called by ultimate.bulk when a WGCNA-capable R environment is available.

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
MIN_SAMPLES <- 15
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  cat('{"status":"skipped","reason":"dependency_missing:WGCNA"}\n'); quit(status = 0)
}

matrix_path <- args[["matrix"]]
samples_path <- args[["samples"]]
modules_path <- args[["modules"]]
trait_path <- args[["trait"]]
object_path <- args[["object"]]
trait_column <- if (is.null(args[["trait-column"]])) "" else args[["trait-column"]]

for (path in c(matrix_path, samples_path)) {
  if (is.null(path) || !file.exists(path)) {
    cat(sprintf('{"status":"skipped","reason":"missing_input:%s"}\n', path)); quit(status = 0)
  }
}
dat <- as.matrix(read.delim(matrix_path, row.names = 1, check.names = FALSE))
samples <- read.delim(samples_path, stringsAsFactors = FALSE)
samples$sample_id <- as.character(samples$sample_id)
keep <- colnames(dat) %in% samples$sample_id
dat <- dat[, keep, drop = FALSE]
if (ncol(dat) < MIN_SAMPLES) {
  cat(sprintf('{"status":"skipped","reason":"insufficient_samples_for_wgcna:<%d","n_samples":%d}\n', MIN_SAMPLES, ncol(dat)))
  quit(status = 0)
}

library(WGCNA)
enableWGCNAThreads()
datExpr <- t(dat)
datExpr <- datExpr[, apply(datExpr, 2, function(x) sum(x != 0) > 0.1 * nrow(datExpr)), drop = FALSE]
if (ncol(datExpr) < 2) {
  cat('{"status":"skipped","reason":"too_few_variable_features"}\n'); quit(status = 0)
}
sft <- pickSoftThreshold(datExpr, powerVector = c(1:20), verbose = FALSE)
softPower <- sft$powerEstimate
if (is.null(softPower) || is.na(softPower)) softPower <- 6
net <- blockwiseModules(datExpr, power = softPower, TOMType = "unsigned",
                        minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE, verbose = FALSE)
moduleColors <- labels2colors(net$colors)
modules <- data.frame(feature_id = colnames(datExpr), module = moduleColors)
utils::write.table(modules, modules_path, sep = "\t", quote = FALSE, row.names = FALSE)

# Module-trait correlations when a numeric trait column is supplied.
trait_rows <- list()
if (nzchar(trait_column) && trait_column %in% colnames(samples)) {
  trait <- suppressWarnings(as.numeric(samples[[trait_column]][match(rownames(datExpr), samples$sample_id)]))
  if (!all(is.na(trait))) {
    for (mod in unique(moduleColors)) {
      if (mod == "grey") next
      eigengene <- as.numeric(moduleEigengenes(datExpr, moduleColors)$eigengenes[[paste0("ME", mod)]])
      ct <- cor.test(eigengene, trait, use = "complete.obs")
      hub <- names(sort(abs(cor(datExpr[, moduleColors == mod], eigengene, use = "complete.obs")), decreasing = TRUE))[1]
      trait_rows[[length(trait_rows) + 1]] <- data.frame(
        module = mod, n_features = sum(moduleColors == mod),
        trait_correlation = unname(ct$estimate), trait_pvalue = ct$p.value,
        hub_gene = hub
      )
    }
  }
}
if (length(trait_rows) == 0) {
  for (mod in unique(moduleColors)) {
    trait_rows[[length(trait_rows) + 1]] <- data.frame(
      module = mod, n_features = sum(moduleColors == mod),
      trait_correlation = NA_real_, trait_pvalue = NA_real_, hub_gene = NA_character_
    )
  }
}
trait_df <- do.call(rbind, trait_rows)
utils::write.table(trait_df, trait_path, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(modules = modules, net = net, softPower = softPower, datExpr = datExpr), object_path)
cat(sprintf('{"status":"ready","engine":"R_WGCNA","softPower":%d,"n_features":%d,"n_samples":%d}\n',
            softPower, ncol(datExpr), nrow(datExpr)))
