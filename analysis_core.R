# analysis_core.R: rnaRAP - Reproducible RNA-seq Analysis Pipeline

# Load aligned Pasilla counts + metadata
load_example_data <- function() {
  pasCts <- system.file("extdata/pasilla_gene_counts.tsv", package = "pasilla", mustWork = TRUE)
  pasAnno <- system.file("extdata/pasilla_sample_annotation.csv", package = "pasilla", mustWork = TRUE)
  
  cts <- as.matrix(read.csv(pasCts, sep = "\t", row.names = "gene_id"))
  coldata <- read.csv(pasAnno, row.names = 1)
  rownames(coldata) <- sub("fb", "", rownames(coldata))  # Strip 'fb' suffix
  
  common <- intersect(colnames(cts), rownames(coldata))
  cts <- cts[, common, drop = FALSE]
  coldata <- coldata[common, , drop = FALSE]
  
  stopifnot(identical(colnames(cts), rownames(coldata)))
  list(cts = cts, coldata = coldata)
}

#' Run DESeq2 on condition (enhanced: size factors, cooks filter)
run_deseq2_condition <- function(cts, coldata) {
  stopifnot("condition" %in% colnames(coldata))
  coldata$condition <- factor(coldata$condition)
  
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
  dds <- DESeq(dds)
  results(dds)
}

#' Run full multi-factor DESeq2 (~ condition + type)
run_deseq2_full <- function(cts, coldata) {
  stopifnot(c("condition", "type") %in% colnames(coldata))
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition + type)
  dds <- DESeq(dds)
  
  # Get actual coefficient names from the model
  coef_names <- resultsNames(dds)
  
  # Find condition and type coefficients (usually 2nd and 3rd positions)
  cond_coef <- grep("^condition", coef_names, value = TRUE)[1]
  type_coef <- grep("^type", coef_names, value = TRUE)[1]
  
  # LFC shrinkage for better visualization
  res_cond <- lfcShrink(dds, coef = cond_coef, type = "apeglm")
  res_type <- lfcShrink(dds, coef = type_coef, type = "apeglm")
  list(condition = res_cond, type = res_type, dds = dds)
}

#' Enhanced MA plot (with LFC shrinkage)
plot_ma_enhanced <- function(res, outpath = "reports/ma_plot.png", alpha = 0.05) {
  res_df <- as.data.frame(res)
  res_df$log10baseMean <- log10(res_df$baseMean + 1)
  sig <- res_df$padj < alpha & !is.na(res_df$padj)
  
  dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)
  png(outpath, width = 1000, height = 800)
  smoothScatter(res_df$log10baseMean, res_df$log2FoldChange,
                xlab = "log10(baseMean + 1)", ylab = "log2FoldChange",
                main = "Enhanced MA Plot (LFC Shrinkage)")
  abline(h = 0, col = "red", lwd = 2)
  points(res_df$log10baseMean[sig], res_df$log2FoldChange[sig], pch = 20, col = "blue", cex = 0.8)
  dev.off()
  invisible(outpath)
}

#' Volcano plot (interactive-ready)
plot_volcano <- function(res, outpath = "reports/volcano.png", lfc_threshold = 1, padj_threshold = 0.05) {
  res_df <- as.data.frame(res)
  res_df$significant <- res_df$padj < padj_threshold & abs(res_df$log2FoldChange) > lfc_threshold
  
  dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)
  png(outpath, width = 1000, height = 800)
  with(res_df,
       plot(log2FoldChange, -log10(padj), pch = 20, main = "Volcano Plot",
            xlab = "log2FoldChange", ylab = "-log10(padj)"))
  abline(v = c(-lfc_threshold, lfc_threshold), col = "red", lwd = 2)
  abline(h = -log10(padj_threshold), col = "red", lwd = 2)
  with(res_df, points(log2FoldChange[significant], -log10(padj)[significant],
                      pch = 20, col = "blue", cex = 1.2))
  dev.off()
  invisible(outpath)
}

#' PCA plot (samples)
plot_pca <- function(dds, outpath = "reports/pca.png", intgroup = "condition") {
  vsd <- vst(dds, blind = FALSE)
  pcaData <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)
  ggplot(pcaData, aes(PC1, PC2, color = .data[[intgroup]])) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA: Sample Clustering by Condition") +
    theme_minimal()
  ggsave(outpath, width = 10, height = 8)
  invisible(outpath)
}

#' Top significant genes table (Markdown + CSV)
top_genes_table <- function(res, n = 20, outdir = "reports") {
  res_df <- as.data.frame(res) |>
    tibble::rownames_to_column("gene") |>
    dplyr::filter(padj < 0.05 & !is.na(padj)) |>
    dplyr::arrange(padj) |>
    dplyr::select(gene, baseMean, log2FoldChange, padj) |>
    dplyr::slice_head(n = n)
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  write.csv(res_df, file.path(outdir, "top_genes.csv"), row.names = FALSE)
  
  # Markdown table for Quarto
  cat("## Top 20 Significant Genes\n\n", file = file.path(outdir, "top_genes.md"))
  cat(knitr::kable(res_df, digits = 3), file = file.path(outdir, "top_genes.md"), append = TRUE)
  invisible(res_df)
}

#' Full session info for reproducibility
save_session_info <- function(outpath = "reports/session_info.txt") {
  dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)
  sink(outpath)
  sessionInfo()
  sink()
  invisible(outpath)
}

#' Plot dispersion estimates
plot_dispersion <- function(dds, outpath = "reports/dispersion.png") {
  dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)
  png(outpath, width = 1000, height = 800)
  plotDispEsts(dds)
  dev.off()
  invisible(outpath)
}

#' P-value distribution histogram
plot_pvalue_dist <- function(res, outpath = "reports/pvalue_dist.png") {
  res_df <- as.data.frame(res) |> dplyr::filter(!is.na(pvalue))
  
  dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)
  png(outpath, width = 1000, height = 800)
  hist(res_df$pvalue, breaks = 50, col = "#3498DB", border = "white",
       main = "Distribution of P-values", xlab = "P-value", ylab = "Frequency")
  dev.off()
  invisible(outpath)
}

#' Heatmap of top genes
plot_heatmap_top_genes <- function(dds, res, n = 30, outpath = "reports/heatmap_top_genes.png") {
  library(pheatmap)
  
  # Get top N genes
  res_df <- as.data.frame(res)
  top_genes <- rownames(res_df)[order(res_df$padj)][1:n]
  top_genes <- top_genes[!is.na(top_genes)]
  
  # Get variance stabilized data
  vsd <- vst(dds, blind = FALSE)
  mat <- assay(vsd)[top_genes, ]
  mat <- t(scale(t(mat)))  # Z-score normalization
  
  # Annotation
  annotation_col <- data.frame(
    Condition = dds$condition,
    Type = dds$type
  )
  rownames(annotation_col) <- colnames(mat)
  
  dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)
  png(outpath, width = 1200, height = 1400)
  pheatmap(mat,
           annotation_col = annotation_col,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           color = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
           main = paste("Top", n, "Differentially Expressed Genes"),
           fontsize_row = 8,
           fontsize_col = 9)
  dev.off()
  invisible(outpath)
}

#' Sample distance heatmap
plot_sample_distances <- function(dds, outpath = "reports/sample_distances.png") {
  library(pheatmap)
  library(RColorBrewer)
  
  vsd <- vst(dds, blind = FALSE)
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep = " - ")
  colnames(sampleDistMatrix) <- NULL
  
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)
  png(outpath, width = 1000, height = 900)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors,
           main = "Sample-to-Sample Distances")
  dev.off()
  invisible(outpath)
}

#' Expression boxplots for top genes
plot_top_gene_boxplots <- function(dds, res, n = 6, outpath = "reports/top_genes_boxplots.png") {
  library(ggplot2)
  library(tidyr)
  
  # Get top N genes
  res_df <- as.data.frame(res)
  top_genes <- rownames(res_df)[order(res_df$padj)][1:n]
  top_genes <- top_genes[!is.na(top_genes)]
  
  # Get variance stabilized data
  vsd <- vst(dds, blind = FALSE)
  counts_long <- assay(vsd)[top_genes, ] |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "sample", values_to = "expression") |>
    dplyr::left_join(
      data.frame(
        sample = colnames(dds),
        condition = dds$condition
      ),
      by = "sample"
    )
  
  p <- ggplot(counts_long, aes(x = condition, y = expression, fill = condition)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    facet_wrap(~gene, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("#3498DB", "#E74C3C")) +
    labs(
      title = paste("Expression of Top", n, "Differentially Expressed Genes"),
      subtitle = "Variance-stabilized expression values",
      x = "Condition",
      y = "Normalized Expression"
    ) +
    theme_minimal() +
    theme(legend.position = "top",
          strip.background = element_rect(fill = "#ECF0F1"))
  
  dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)
  ggsave(outpath, p, width = 12, height = 8)
  invisible(outpath)
}

#' Save analysis results as RDS for later use
save_analysis_results <- function(full_results, data_list, outpath = "reports/analysis_results.rds") {
  dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)
  results <- list(
    full_results = full_results,
    data_list = data_list
  )
  saveRDS(results, outpath)
  invisible(outpath)
}