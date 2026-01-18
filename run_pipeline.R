#!/usr/bin/env Rscript
# run_pipeline.R: rnaRAP - Reproducible RNA-seq Analysis Pipeline

# ===== SETUP =====
cat("=== RNA-seq RAP: Reproducible Analysis Pipeline ===\n")
cat("Starting analysis:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(pasilla)
  library(ggplot2)
  library(dplyr)
  library(knitr)
  library(pheatmap)
  library(RColorBrewer)
  library(tidyr)
})

# Source analysis functions
source("analysis_core.R")

# ===== 1. LOAD DATA =====
cat("1. Loading Pasilla example data...\n")
data_list <- load_example_data()
cat("   Loaded", nrow(data_list$cts), "genes across", ncol(data_list$cts), "samples\n\n")

# ===== 2. RUN DESEQ2 ANALYSES =====
cat("2. Running DESeq2 analyses...\n")

# Simple condition analysis
cat("   - Running condition-only analysis...\n")
res_condition <- run_deseq2_condition(data_list$cts, data_list$coldata)

# Full multi-factor analysis
cat("   - Running full multi-factor analysis (condition + type)...\n")
full_results <- run_deseq2_full(data_list$cts, data_list$coldata)
cat("   Analysis complete!\n\n")

# ===== 3. GENERATE SUMMARIES =====
cat("3. Generating result summaries...\n")
cat("\n=== Condition Analysis Summary ===\n")
print(summary(full_results$condition))

cat("\n=== Type Analysis Summary ===\n")
print(summary(full_results$type))

# Count significant genes
n_sig_condition <- sum(full_results$condition$padj < 0.05, na.rm = TRUE)
n_sig_type <- sum(full_results$type$padj < 0.05, na.rm = TRUE)
cat("\nSignificant genes (padj < 0.05):\n")
cat("  - Condition effect:", n_sig_condition, "\n")
cat("  - Type effect:", n_sig_type, "\n\n")

# ===== 4. CREATE VISUALIZATIONS =====
cat("4. Creating visualizations...\n")

# Core plots
cat("   - Generating MA plot...\n")
plot_ma_enhanced(full_results$condition, outpath = "reports/ma_plot_condition.png")

cat("   - Generating volcano plot...\n")
plot_volcano(full_results$condition, outpath = "reports/volcano_condition.png")

cat("   - Generating PCA plot...\n")
plot_pca(full_results$dds, outpath = "reports/pca_condition.png", intgroup = "condition")

# Additional quality control and exploratory plots
cat("   - Generating dispersion plot...\n")
plot_dispersion(full_results$dds, outpath = "reports/dispersion.png")

cat("   - Generating p-value distribution...\n")
plot_pvalue_dist(full_results$condition, outpath = "reports/pvalue_dist.png")

cat("   - Generating heatmap of top genes...\n")
plot_heatmap_top_genes(full_results$dds, full_results$condition, n = 30, 
                       outpath = "reports/heatmap_top_genes.png")

cat("   - Generating sample distance heatmap...\n")
plot_sample_distances(full_results$dds, outpath = "reports/sample_distances.png")

cat("   - Generating top gene boxplots...\n")
plot_top_gene_boxplots(full_results$dds, full_results$condition, n = 6,
                       outpath = "reports/top_genes_boxplots.png")

cat("   All plots saved to reports/ directory\n\n")

# ===== 5. EXPORT TOP GENES =====
cat("5. Exporting top differentially expressed genes...\n")
top_genes_condition <- top_genes_table(full_results$condition, n = 20, outdir = "reports")
cat("   Top genes table saved to reports/top_genes.csv\n\n")

# ===== 6. SAVE ANALYSIS RESULTS =====
cat("6. Saving analysis results for report generation...\n")
save_analysis_results(full_results, data_list, outpath = "reports/analysis_results.rds")
cat("   Analysis results saved to reports/analysis_results.rds\n\n")

# ===== 7. SAVE SESSION INFO =====
cat("7. Saving session information for reproducibility...\n")
save_session_info(outpath = "reports/session_info.txt")
cat("   Session info saved to reports/session_info.txt\n\n")

# ===== 8. RENDER QUARTO REPORT =====
cat("8. Rendering Quarto report as landing page...\n")
if (file.exists("rnaseq_report.qmd")) {
  cat("   Found rnaseq_report.qmd, rendering to index.html...\n")
  
  # Render to index.html directly
  quarto::quarto_render("rnaseq_report.qmd", output_file = "index.html")
  
  # Check if HTML was created in root
  if (file.exists("index.html")) {
    cat("   Successfully rendered to index.html\n")
    # Move to reports folder
    file.rename("index.html", "reports/index.html")
    cat("   Moved to reports/index.html\n\n")
  } else {
    cat("   WARNING: index.html was not created\n")
    # Fallback: render to default name and rename
    if (file.exists("rnaseq_report.html")) {
      file.rename("rnaseq_report.html", "reports/index.html")
      cat("   Renamed rnaseq_report.html to index.html\n\n")
    }
  }
}

# ===== COMPLETION =====
cat("=== PIPELINE COMPLETE ===\n")
cat("Analysis finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Outputs:\n")
cat("  - MA plot: reports/ma_plot_condition.png\n")
cat("  - Volcano plot: reports/volcano_condition.png\n")
cat("  - PCA plot: reports/pca_condition.png\n")
cat("  - Dispersion plot: reports/dispersion.png\n")
cat("  - P-value distribution: reports/pvalue_dist.png\n")
cat("  - Heatmap (top 30 genes): reports/heatmap_top_genes.png\n")
cat("  - Sample distances: reports/sample_distances.png\n")
cat("  - Top gene boxplots: reports/top_genes_boxplots.png\n")
cat("  - Top genes table: reports/top_genes.csv\n")
cat("  - Top genes markdown: reports/top_genes.md\n")
cat("  - Session info: reports/session_info.txt\n")
cat("  - Analysis results (RDS): reports/analysis_results.rds\n")
cat("  - HTML Report: reports/index.html\n")
