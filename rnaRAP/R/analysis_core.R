# R/analysis_core.R

#' Load pasilla counts and metadata from DESeq2
#' @return list(cts, coldata)
#' @export

load_example_data <- function() {
  pasCts <- system.file(
    "extdata",
    "pasilla_gene_counts.tsv",
    package = "pasilla",
    mustWork = TRUE
  )
  pasAnno <- system.file(
    "extdata",
    "pasilla_sample_annotation.csv",
    package = "pasilla",
    mustWork = TRUE
  )
  
  cts <- as.matrix(
    read.csv(pasCts, sep = "\t", row.names = "gene_id")
  )
  coldata <- read.csv(pasAnno, row.names = 1)
  
  # Strip 'fb' suffix from rownames of coldata to match counts
  rownames(coldata) <- sub("fb$", "", rownames(coldata))
  
  # Align by common sample names
  common <- intersect(colnames(cts), rownames(coldata))
  cts <- cts[, common, drop = FALSE]
  coldata <- coldata[common, , drop = FALSE]
  
  # Final check
  stopifnot(identical(colnames(cts), rownames(coldata)))
  
  list(cts = cts, coldata = coldata)
}

#' Run simple DESeq2 analysis on pasilla condition
#' @param cts counts matrix
#' @param coldata sample metadata with 'condition' column
#' @return DESeq2 results object
#' @export
run_deseq2_condition <- function(cts, coldata) {
  stopifnot("condition" %in% colnames(coldata))
  
  coldata$condition <- factor(coldata$condition)
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~ condition
  )
  dds <- DESeq2::DESeq(dds)
  DESeq2::results(dds)
}

#' Save MA plot to file
#' @param res DESeq2 results
#' @param out_path where to save PNG
#' @export
plot_ma_simple <- function(res, out_path = "reports/ma_plot.png") {
  res_df <- as.data.frame(res)
  res_df$log10baseMean <- log10(res_df$baseMean + 1)
  
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  
  png(out_path, width = 1000, height = 800)
  smoothScatter(
    res_df$log10baseMean,
    res_df$log2FoldChange,
    xlab = "log10(baseMean)",
    ylab = "log2FC",
    main = "MA plot"
  )
  abline(h = 0, col = "red")
  dev.off()
  
  invisible(out_path)
}
